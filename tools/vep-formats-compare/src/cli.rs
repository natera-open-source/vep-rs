// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The command line.

pub use crate::argparse::ArgumentParser;
use crate::argparse::{
    self, path_name, path_text, standard_error, standard_output, Action, Exit, Kind, Taken,
};
use crate::cache::THREADS_VARIABLE;
use crate::cache_classes::{self, CacheTable};
use crate::classes::{
    self, DEFAULT_REFERENCE_RELEASE, REFERENCE_RELEASES, REGULATORY_FEATURE_ROWS,
};
use crate::coding_ends::terms::CodingEndClass;
use crate::compare_files::compare_files;
use crate::consts::{EXAMPLES_PER_GROUP, REFERENCE_PRESENT};
use crate::error::{Error, Result};
use crate::fasta::ReferenceSequence;
use crate::field_class;
use crate::fields::OrganelleCodonTables;
use crate::input::{input_records, read_text};
use crate::manifest::Documented;
use crate::nearest::NearestSites;
use crate::parse::parquet::parquet_as_tab;
use crate::pyfmt;
use crate::pyjson::{self, Value};
use crate::remainder::FrontEnd;
use crate::report::write_report;
use crate::rules::downstream::DownstreamSequences;
use crate::rules::gwas::GwasCatalogPairs;
use crate::rules::intron_ends::IntronEnds;
use crate::rules::revel::RevelTable;
use crate::rules::sequence_edges::SequenceEdges;
use crate::rules::splice_windows::SpliceWindows;
use crate::rules::spliceai::SpliceAiTable;
use crate::rules::structural::StructuralVariants;
use crate::rules::transcript_edges::TranscriptEdges;
use crate::rules::utrannotator::UtrannotatorSeedSweep;
use crate::rules::Side;
use crate::self_parity::SelfParity;
use crate::warnings::compare_warnings;
use crate::COMPONENT;
use indexmap::IndexMap;
use std::collections::{BTreeMap, BTreeSet};
use std::ffi::OsString;
use std::io::Write;
use std::path::Path;
use std::path::PathBuf;
use std::sync::Arc;

/// The name the usage and the usage errors give the program.
pub const PROG: &str = "vep-formats-compare";

/// The parser's usage block and its whole help at 80 columns.
const USAGE: &str = include_str!("../data/usage.txt");
const HELP: &str = include_str!("../data/help.txt");

/// Runs the comparison `argv` names and returns the exit code: 0 for a green verdict, 1
/// for a red one, 2 for a refused command line, an unreadable input and a report that
/// cannot be written. Nothing is written, to the report or to a standard stream, before
/// the whole report exists; the report file is written before the standard output.
///
/// A command line whose first argument is `field-class` is that command's
/// (`field_class::main`).
pub fn main(argv: &[OsString]) -> i32 {
    if argv
        .first()
        .is_some_and(|first| first == field_class::COMMAND)
    {
        return field_class::main(&argv[1..]);
    }
    let parser = build_parser();
    let args = match parse_args(&parser, argv) {
        Ok(args) => args,
        Err(Exit::Help) => {
            let _ = std::io::stdout().write_all(parser.format_help().as_bytes());
            return 0;
        }
        Err(Exit::Error(message)) => {
            let text = standard_error(&parser.error_text(&message));
            let _ = std::io::stderr().write_all(text.as_bytes());
            return 2;
        }
    };
    let written = thread_count()
        .and_then(|_| run(&args, None))
        .and_then(|report| write_report(&report, args.report.as_deref()).map(|()| report));
    let report = match written {
        Ok(report) => report,
        Err(error) => {
            let text = standard_error(&format!("ERROR: [{COMPONENT}] {error}\n"));
            let _ = std::io::stderr().write_all(text.as_bytes());
            return 2;
        }
    };
    let text_of = |name: &str| report.get(name).map(pyjson::py_str).unwrap_or_default();
    let label = match report.get("format") {
        Some(format) if pyjson::truthy(format) => pyjson::py_str(format),
        _ => text_of("mode"),
    };
    let verdict = text_of("verdict");
    let mut lines = format!("{label}: {verdict}\n");
    for reason in report
        .get("reasons")
        .and_then(Value::as_list)
        .unwrap_or(&[])
    {
        lines.push_str(&format!("  {}\n", pyjson::py_str(reason)));
    }
    let _ = std::io::stdout().write_all(&standard_output(&lines));
    if verdict == "green" {
        0
    } else {
        1
    }
}

/// The threads a comparison works on: the value of `VEP_FORMATS_COMPARE_THREADS`, which
/// is refused unless it is a whole number of 1 or more, or the cores of the host where it
/// is unset or empty.
pub fn thread_count() -> Result<usize> {
    let available = || std::thread::available_parallelism().map_or(1, |cores| cores.get());
    match std::env::var_os(THREADS_VARIABLE) {
        None => Ok(available()),
        Some(value) if value.is_empty() => Ok(available()),
        Some(value) => value
            .to_str()
            .and_then(|text| text.trim().parse::<usize>().ok())
            .filter(|threads| *threads >= 1)
            .ok_or_else(|| {
                Error::comparator(format!(
                    "{THREADS_VARIABLE}={}: expected a whole number of threads, 1 or more",
                    argparse::text_of(&value)
                ))
            }),
    }
}

/// The parsed command line: one field an option, as `build_parser` declares them.
#[derive(Debug)]
pub struct Namespace {
    pub self_parity: bool,
    pub warnings: bool,
    pub format: Option<String>,
    pub expected: Option<PathBuf>,
    pub actual: Option<PathBuf>,
    pub input: Option<PathBuf>,
    pub manifest: Option<PathBuf>,
    pub run: Option<String>,
    pub reference_sequence: Option<String>,
    pub masks: Option<String>,
    /// The Ensembl VEP release that wrote the reference output; `None` reads as
    /// `DEFAULT_REFERENCE_RELEASE`.
    pub reference_release: Option<String>,
    pub keep_only: Option<String>,
    pub floors: Option<String>,
    pub field_pairs: Option<Vec<String>>,
    pub nearest_sites: Option<PathBuf>,
    pub json_cache: Option<PathBuf>,
    pub reference_fasta: Option<PathBuf>,
    pub gwas_catalog: Option<PathBuf>,
    pub revel_table: Option<PathBuf>,
    pub revel_no_match: bool,
    pub spliceai_table: Option<Vec<PathBuf>>,
    pub utrannotator_sweep: Option<Vec<PathBuf>>,
    pub default: Option<PathBuf>,
    pub tab: Option<PathBuf>,
    pub vcf: Option<PathBuf>,
    pub json: Option<PathBuf>,
    pub parquet: Option<PathBuf>,
    pub expected_skipped: Option<PathBuf>,
    pub actual_skipped: Option<PathBuf>,
    pub unordered: bool,
    pub report: Option<PathBuf>,
    pub examples: i64,
}

impl Default for Namespace {
    /// The command line without an option: no switch, no value, `EXAMPLES_PER_GROUP`
    /// examples.
    fn default() -> Self {
        Namespace {
            self_parity: false,
            warnings: false,
            format: None,
            expected: None,
            actual: None,
            input: None,
            manifest: None,
            run: None,
            reference_sequence: None,
            masks: None,
            reference_release: None,
            keep_only: None,
            floors: None,
            field_pairs: None,
            nearest_sites: None,
            json_cache: None,
            reference_fasta: None,
            gwas_catalog: None,
            revel_table: None,
            revel_no_match: false,
            spliceai_table: None,
            utrannotator_sweep: None,
            default: None,
            tab: None,
            vcf: None,
            json: None,
            parquet: None,
            expected_skipped: None,
            actual_skipped: None,
            unordered: false,
            report: None,
            examples: EXAMPLES_PER_GROUP,
        }
    }
}

/// The annotation formats, and with the recoder's JSON every format the comparison mode takes.
pub const VEP_FORMATS: &[&str] = &["default", "tab", "vcf", "json", "parquet"];

pub const FORMATS: &[&str] = &["default", "tab", "vcf", "json", "parquet", "recoder"];

/// The value of `--reference-sequence` for a candidate run that had no reference sequence.
pub const REFERENCE_ABSENT: &str = "absent";

/// `COL=MIN,...` into a column-to-minimum map. An item is stripped of the blanks around
/// it, not of those around its `=`; a column named twice keeps its place and takes the
/// later minimum.
///
/// A minimum of 2**63 or more, which no count reaches, ends the run once the whole list
/// was read, unless a later item gives its column another minimum.
pub fn parse_floors(spec: Option<&str>) -> Result<IndexMap<String, i64>> {
    let mut floors = IndexMap::new();
    let mut beyond: IndexMap<&str, &str> = IndexMap::new();
    for item in spec.unwrap_or("").split(',') {
        let item = pyfmt::strip(item);
        if item.is_empty() {
            continue;
        }
        let (name, sep, value) = pyfmt::partition(item, "=");
        if sep.is_empty() || name.is_empty() || !pyfmt::is_digit(value) {
            return Err(Error::comparator(format!(
                "--floors expects COL=MIN with a non-negative integer, got {}",
                argparse::repr(item)
            )));
        }
        // `str.isdigit` passes digits `int` does not read, and more of them than it converts.
        let minimum = match pyfmt::int_checked(value) {
            Ok(decimal) => decimal,
            Err(error) => crate::uncaught!("{error}"),
        };
        match minimum.parse::<i64>() {
            Ok(minimum) => {
                floors.insert(name.to_string(), minimum);
                beyond.shift_remove(name);
            }
            Err(_) => {
                floors.insert(name.to_string(), i64::MAX);
                beyond.insert(name, item);
            }
        }
    }
    if let Some(item) = beyond.values().next() {
        crate::uncaught!("--floors {item}: a minimum of 2**63 or more is not read");
    }
    Ok(floors)
}

/// `COL=EXPECTED=>ACTUAL` items, several per argument separated by `;`. A blank item is
/// skipped; an item keeps the blanks around its column and its values.
pub fn parse_field_pairs(
    specs: Option<&[String]>,
) -> Result<IndexMap<String, BTreeSet<(String, String)>>> {
    let mut pairs: IndexMap<String, BTreeSet<(String, String)>> = IndexMap::new();
    for spec in specs.unwrap_or_default() {
        for item in spec.split(';') {
            if pyfmt::strip(item).is_empty() {
                continue;
            }
            let (name, sep, rest) = pyfmt::partition(item, "=");
            let (expected, arrow, actual) = pyfmt::partition(rest, "=>");
            if sep.is_empty() || name.is_empty() || arrow.is_empty() {
                return Err(Error::comparator(format!(
                    "--field-pairs expects COL=EXPECTED=>ACTUAL, got {}",
                    argparse::repr(item)
                )));
            }
            pairs
                .entry(name.to_string())
                .or_default()
                .insert((expected.to_string(), actual.to_string()));
        }
    }
    Ok(pairs)
}

/// The classes of a comma-separated list, each stripped of the blanks around it, the empty
/// ones left out.
pub fn parse_masks(spec: Option<&str>) -> BTreeSet<String> {
    spec.unwrap_or("")
        .split(',')
        .map(pyfmt::strip)
        .filter(|class| !class.is_empty())
        .map(str::to_string)
        .collect()
}

/// The classes of `named` a comparison applies where Ensembl VEP `release` wrote its
/// reference output, and the classes of `named` that release corrects
/// (`classes::corrected_under`), each with the first release that corrects its defect.
pub fn masks_under(
    named: BTreeSet<String>,
    release: &str,
) -> (BTreeSet<String>, BTreeMap<String, String>) {
    let mut corrected = BTreeMap::new();
    let applied = named
        .into_iter()
        .filter(|klass| match classes::corrected_under(klass, release) {
            Some(first) => {
                corrected.insert(klass.clone(), first.to_string());
                false
            }
            None => true,
        })
        .collect();
    (applied, corrected)
}

/// The manifest at `path`, which must be a JSON object; `None` without a path.
pub fn load_manifest(path: Option<&Path>) -> Result<Option<Value>> {
    let Some(path) = path else {
        return Ok(None);
    };
    let manifest = match pyjson::loads(&read_text(path)?) {
        Ok(manifest) => manifest,
        Err(error) if error.is_value_error() => {
            return Err(Error::comparator(format!(
                "{} is not JSON: {error}",
                path_text(path)
            )))
        }
        Err(error) => crate::uncaught!("{}: {error}", path_text(path)),
    };
    if manifest.as_dict().is_none() {
        return Err(Error::comparator(format!(
            "{} is not a manifest object",
            path_text(path)
        )));
    }
    Ok(Some(manifest))
}

/// The comparator's options, in the order the help lists them.
pub fn build_parser() -> ArgumentParser {
    let mut ap = ArgumentParser::new(PROG, USAGE, HELP);
    let switch = |ap: &mut ArgumentParser, option| {
        ap.add_argument(option, Action::StoreTrue, Kind::Text, None)
    };
    let text = |ap: &mut ArgumentParser, option| {
        ap.add_argument(option, Action::Store, Kind::Text, None);
    };
    let choice = |ap: &mut ArgumentParser, option, choices| {
        ap.add_argument(option, Action::Store, Kind::Text, Some(choices));
    };
    let path = |ap: &mut ArgumentParser, option| {
        ap.add_argument(option, Action::Store, Kind::Path, None);
    };
    let self_parity = switch(&mut ap, "--self-parity");
    let warnings = switch(&mut ap, "--warnings");
    ap.add_mutually_exclusive_group(&[self_parity, warnings]);
    choice(&mut ap, "--format", FORMATS);
    path(&mut ap, "--expected");
    path(&mut ap, "--actual");
    path(&mut ap, "--input");
    path(&mut ap, "--manifest");
    text(&mut ap, "--run");
    choice(
        &mut ap,
        "--reference-sequence",
        &[REFERENCE_PRESENT, REFERENCE_ABSENT],
    );
    text(&mut ap, "--masks");
    choice(&mut ap, "--reference-release", REFERENCE_RELEASES);
    choice(&mut ap, "--keep-only", &[REGULATORY_FEATURE_ROWS]);
    text(&mut ap, "--floors");
    ap.add_argument("--field-pairs", Action::Append, Kind::Text, None);
    path(&mut ap, "--nearest-sites");
    path(&mut ap, "--json-cache");
    path(&mut ap, "--reference-fasta");
    path(&mut ap, "--gwas-catalog");
    path(&mut ap, "--revel-table");
    switch(&mut ap, "--revel-no-match");
    ap.add_argument("--spliceai-table", Action::Append, Kind::Path, None);
    ap.add_argument("--utrannotator-sweep", Action::Append, Kind::Path, None);
    path(&mut ap, "--default");
    path(&mut ap, "--tab");
    path(&mut ap, "--vcf");
    path(&mut ap, "--json");
    path(&mut ap, "--parquet");
    path(&mut ap, "--expected-skipped");
    path(&mut ap, "--actual-skipped");
    switch(&mut ap, "--unordered");
    path(&mut ap, "--report");
    ap.add_argument("--examples", Action::Store, Kind::Int, None);
    ap
}

/// The command line `argv` as `parser` reads it, or how the parser ends the run: an
/// option given twice keeps its last value, and `--field-pairs`, `--spliceai-table` and
/// `--utrannotator-sweep` keep every value in order.
pub fn parse_args(
    parser: &ArgumentParser,
    argv: &[OsString],
) -> std::result::Result<Namespace, Exit> {
    use argparse::Value as Given;
    let mut args = Namespace::default();
    for Taken { argument, value } in parser.parse_args(argv)? {
        match (parser.dest(argument), value) {
            ("self_parity", _) => args.self_parity = true,
            ("warnings", _) => args.warnings = true,
            ("format", Given::Text(text)) => args.format = Some(text),
            ("expected", Given::Path(path)) => args.expected = Some(path),
            ("actual", Given::Path(path)) => args.actual = Some(path),
            ("input", Given::Path(path)) => args.input = Some(path),
            ("manifest", Given::Path(path)) => args.manifest = Some(path),
            ("run", Given::Text(text)) => args.run = Some(text),
            ("reference_sequence", Given::Text(text)) => args.reference_sequence = Some(text),
            ("masks", Given::Text(text)) => args.masks = Some(text),
            ("reference_release", Given::Text(text)) => args.reference_release = Some(text),
            ("keep_only", Given::Text(text)) => args.keep_only = Some(text),
            ("floors", Given::Text(text)) => args.floors = Some(text),
            ("field_pairs", Given::Text(text)) => {
                args.field_pairs.get_or_insert_default().push(text)
            }
            ("nearest_sites", Given::Path(path)) => args.nearest_sites = Some(path),
            ("json_cache", Given::Path(path)) => args.json_cache = Some(path),
            ("reference_fasta", Given::Path(path)) => args.reference_fasta = Some(path),
            ("gwas_catalog", Given::Path(path)) => args.gwas_catalog = Some(path),
            ("revel_table", Given::Path(path)) => args.revel_table = Some(path),
            ("revel_no_match", _) => args.revel_no_match = true,
            ("spliceai_table", Given::Path(path)) => {
                args.spliceai_table.get_or_insert_default().push(path);
            }
            ("utrannotator_sweep", Given::Path(path)) => {
                args.utrannotator_sweep.get_or_insert_default().push(path);
            }
            ("default", Given::Path(path)) => args.default = Some(path),
            ("tab", Given::Path(path)) => args.tab = Some(path),
            ("vcf", Given::Path(path)) => args.vcf = Some(path),
            ("json", Given::Path(path)) => args.json = Some(path),
            ("parquet", Given::Path(path)) => args.parquet = Some(path),
            ("expected_skipped", Given::Path(path)) => args.expected_skipped = Some(path),
            ("actual_skipped", Given::Path(path)) => args.actual_skipped = Some(path),
            ("unordered", _) => args.unordered = true,
            ("report", Given::Path(path)) => args.report = Some(path),
            ("examples", Given::Int(examples)) => args.examples = examples,
            _ => {}
        }
    }
    Ok(args)
}

/// The report of the mode `args` names. The inputs and the side tables are read in one
/// fixed order, so that of several things wrong with a command line the one reported is
/// the first in that order.
///
/// A class `--masks` names whose defect the reference release corrects (`masks_under`) is
/// taken out of the named set before the manifest, a side table or a rule reads the set,
/// in the comparison mode and in the warnings mode.
///
/// `_front_end` is not read: `compare_files` opens the stream of record groups itself.
pub fn run(args: &Namespace, _front_end: Option<&mut FrontEnd>) -> Result<Value> {
    if args.keep_only.is_some() && (args.warnings || args.self_parity) {
        return Err(Error::comparator(
            "--keep-only is available only in --format comparison mode",
        ));
    }
    let release = args
        .reference_release
        .as_deref()
        .unwrap_or(DEFAULT_REFERENCE_RELEASE);
    if !REFERENCE_RELEASES.contains(&release) {
        return Err(Error::comparator(format!(
            "--reference-release expects one of {}, got {}",
            REFERENCE_RELEASES.join(", "),
            argparse::repr(release)
        )));
    }
    if args.warnings {
        let (Some(expected), Some(actual)) = (&args.expected, &args.actual) else {
            return Err(Error::comparator(
                "--warnings needs --expected and --actual",
            ));
        };
        let (masks, _) = masks_under(parse_masks(args.masks.as_deref()), release);
        let masks: Vec<&str> = masks.iter().map(String::as_str).collect();
        return compare_warnings(
            expected,
            actual,
            args.expected_skipped.as_deref(),
            args.actual_skipped.as_deref(),
            &masks,
        );
    }
    if args.self_parity {
        return self_parity(args);
    }
    let (Some(format), Some(expected), Some(actual)) =
        (args.format.as_deref(), &args.expected, &args.actual)
    else {
        return Err(Error::comparator(
            "--format, --expected and --actual are required",
        ));
    };
    let (masks, corrected) = masks_under(parse_masks(args.masks.as_deref()), release);
    let manifest = load_manifest(args.manifest.as_deref())?;
    let input_name = args.input.as_deref().map(path_name);
    let mut documented = Documented::from_manifest(
        manifest.as_ref(),
        &masks,
        input_name.as_deref(),
        args.run.as_deref(),
        args.reference_sequence.as_deref(),
    )?;
    documented.corrected = corrected;
    let tables = classes::tables();
    let names_one_of = |classes: &[&str]| classes.iter().any(|class| masks.contains(*class));
    let recoder = format == "recoder";
    let cache = args.json_cache.as_deref();
    let structural_variants_requested =
        !recoder && names_one_of(tables.structural_variant_classes());
    // The records are read where a parser or class consumes them; a default or tab
    // comparison names its input for the manifest's flag set alone and reads none of its
    // lines, but for the records a class of `RECORD_REFERENCE_CLASSES` reads by position.
    let reads_records = matches!(format, "vcf" | "json")
        || structural_variants_requested
        || names_one_of(tables.allele_reading_classes());
    // The side tables that place an entry by its input line share the lines.
    let inputs: Option<Arc<[String]>> = match &args.input {
        Some(input) if reads_records => Some(input_records(input)?.into()),
        _ => None,
    };
    let floors = parse_floors(args.floors.as_deref())?;
    let field_pairs = parse_field_pairs(args.field_pairs.as_deref())?;
    let nearest_sites = args
        .nearest_sites
        .as_deref()
        .map(NearestSites::from_cache)
        .transpose()?;
    let mut structural_variants = None;
    if structural_variants_requested {
        let classes: Vec<&str> = masks.iter().map(String::as_str).collect();
        let table = StructuralVariants::new(&classes, args.input.as_deref(), cache)?;
        if let (Some(inputs), Some(input)) = (&inputs, &args.input) {
            if table.input_variants.len() != inputs.len() {
                return Err(Error::comparator(format!(
                    "{}: {} data lines, {} records the structural-variant comparator reads; \
                     a record needs CHROM, POS, ID, REF and ALT",
                    path_text(input),
                    inputs.len(),
                    table.input_variants.len()
                )));
            }
        }
        structural_variants = Some(table);
    }
    // Each side table below is built from the candidate's cache where the masks name a
    // class that `data/json_cache_classes.tsv` lists under it.
    let named_of = |table: CacheTable| -> Vec<&'static str> {
        cache_classes::classes_reading(table)
            .iter()
            .copied()
            .filter(|class| masks.contains(*class))
            .collect()
    };
    let reads = |table: CacheTable| !named_of(table).is_empty();
    let mut sequence_edges = None;
    if let Some(fasta) = &args.reference_fasta {
        let edge_classes = named_of(CacheTable::SequenceEdges);
        if !edge_classes.is_empty() {
            let Some(cache) = cache else {
                return Err(Error::comparator(format!(
                    "--masks {} with --reference-fasta needs \
                     --json-cache, the candidate JSON cache",
                    edge_classes.join(",")
                )));
            };
            sequence_edges = Some(SequenceEdges::new(fasta, cache)?);
        }
    }
    // The five classes of a coding end take a pair on the reading of its transcript's
    // cache record alone.
    let coding_end_classes: Vec<&str> = CodingEndClass::IN_ORDER
        .iter()
        .map(|class| class.name())
        .filter(|class| masks.contains(*class))
        .collect();
    if !coding_end_classes.is_empty() && cache.is_none() {
        return Err(Error::comparator(format!(
            "--masks {} needs --json-cache, the candidate JSON cache",
            coding_end_classes.join(",")
        )));
    }
    // The classes of `TRANSCRIPT_READING_CLASSES` take a row whose two consequence sets
    // differ as a class of a coding end takes them, so their comparison reads the coding
    // ends beside the table the file lists them under.
    let reads_the_coding_ends =
        reads(CacheTable::PartialCodonModels) || names_one_of(tables.transcript_reading_classes());
    let mut partial_codon_models = None;
    if reads_the_coding_ends && !recoder {
        partial_codon_models = cache
            .map(crate::rules::partial_codon_models::PartialCodonModels::new)
            .transpose()?;
    }
    let mut length_edits = None;
    if reads(CacheTable::LengthEdits) && !recoder {
        length_edits = cache
            .map(crate::rules::length_edits::LengthEdits::new)
            .transpose()?;
    }
    let mut reference_sequence = None;
    let completed = CodingEndClass::StopCodonCompleted.name();
    if !recoder && partial_codon_models.is_some() && masks.contains(completed) {
        reference_sequence = args
            .reference_fasta
            .as_deref()
            .map(ReferenceSequence::new)
            .transpose()?;
    }
    let mut codon_tables = None;
    if reads(CacheTable::CodonTables) && !recoder {
        codon_tables = cache.map(OrganelleCodonTables::new).transpose()?;
    }
    let mut downstream_sequences = None;
    if !recoder && reads(CacheTable::DownstreamSequences) {
        downstream_sequences = cache.map(DownstreamSequences::new).transpose()?;
    }
    let mut transcript_edges = None;
    if cache_classes::classes_reading(CacheTable::TranscriptEdges)
        .iter()
        .any(|class| documented.masks.contains(*class))
    {
        transcript_edges = cache
            .map(|cache| TranscriptEdges::new(cache, inputs.as_deref()))
            .transpose()?;
    }
    let mut intron_ends = None;
    let intron_gap_classes = named_of(CacheTable::IntronEnds);
    if !intron_gap_classes.is_empty() {
        let Some(cache) = cache else {
            return Err(Error::comparator(format!(
                "--masks {} needs --json-cache, the candidate JSON cache",
                intron_gap_classes.join(",")
            )));
        };
        intron_ends = Some(IntronEnds::sharing(
            cache,
            inputs.clone().unwrap_or_default(),
        )?);
    }
    let mut splice_windows = None;
    if reads(CacheTable::SpliceWindows) && !recoder {
        splice_windows = cache
            .map(|cache| SpliceWindows::sharing(cache, inputs.clone().unwrap_or_default()))
            .transpose()?;
    }
    let gwas_pairs = args
        .gwas_catalog
        .as_deref()
        .map(GwasCatalogPairs::from_catalog)
        .transpose()?;
    let utrannotator_sweep = paths_of(&args.utrannotator_sweep)
        .map(|outputs| UtrannotatorSeedSweep::from_outputs(&outputs))
        .transpose()?;
    let input_path = match &args.input {
        Some(input) if names_one_of(tables.record_reference_classes()) => Some(input.as_path()),
        _ => None,
    };
    let revel_table = args
        .revel_table
        .as_deref()
        .map(|table| RevelTable::new(table, args.revel_no_match))
        .transpose()?;
    let spliceai_table = paths_of(&args.spliceai_table)
        .map(|tables| SpliceAiTable::new(&tables, cache))
        .transpose()?;
    let mut side = Side {
        gwas_pairs,
        utrannotator_sweep,
        codon_tables,
        partial_codon_models,
        length_edits,
        reference_sequence,
        downstream_sequences,
        transcript_edges,
        intron_ends,
        revel_table,
        spliceai_table,
        splice_windows,
        sequence_edges,
        nearest_sites,
        structural_variants,
    };
    compare_files(
        format,
        expected,
        actual,
        inputs.as_deref(),
        &mut documented,
        &masks,
        &floors,
        &field_pairs,
        args.examples,
        args.unordered,
        args.keep_only.as_deref(),
        input_path,
        args.input.as_deref(),
        &mut side,
    )
}

/// The paths of an option that may be given several times, `None` where it was not given.
fn paths_of(paths: &Option<Vec<PathBuf>>) -> Option<Vec<&Path>> {
    paths
        .as_ref()
        .filter(|paths| !paths.is_empty())
        .map(|paths| paths.iter().map(PathBuf::as_path).collect())
}

/// The self-parity report with the files it was made from, in the order default, tab,
/// vcf, json, parquet, input, the ones not given left out.
fn self_parity(args: &Namespace) -> Result<Value> {
    let files = [
        ("default", &args.default),
        ("tab", &args.tab),
        ("vcf", &args.vcf),
        ("json", &args.json),
        ("parquet", &args.parquet),
        ("input", &args.input),
    ];
    let (Some(default), Some(tab), Some(vcf), Some(json)) =
        (&args.default, &args.tab, &args.vcf, &args.json)
    else {
        let missing: Vec<&str> = files[..4]
            .iter()
            .filter(|(_, path)| path.is_none())
            .map(|(name, _)| *name)
            .collect();
        return Err(Error::comparator(format!(
            "--self-parity needs --{}",
            missing.join(", --")
        )));
    };
    let inputs = args.input.as_deref().map(input_records).transpose()?;
    let parquet_text = args
        .parquet
        .as_deref()
        .map(|parquet| parquet_as_tab(parquet, false))
        .transpose()?;
    let default_text = read_text(default)?;
    let tab_text = read_text(tab)?;
    let vcf_text = read_text(vcf)?;
    let json_text = read_text(json)?;
    let mut report = SelfParity::new(
        &default_text,
        &tab_text,
        &vcf_text,
        &json_text,
        parquet_text.as_deref(),
        inputs.as_deref(),
        args.examples,
    )?
    .run();
    let files: IndexMap<String, Value> = files
        .iter()
        .filter_map(|(name, path)| Some((name.to_string(), Value::Str(path_text(path.as_ref()?)))))
        .collect();
    if let Value::Dict(members) = &mut report {
        members.insert("files".to_string(), Value::Dict(files));
    }
    Ok(report)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parsed(args: &[&str]) -> std::result::Result<Namespace, Exit> {
        let argv: Vec<OsString> = args.iter().map(OsString::from).collect();
        parse_args(&build_parser(), &argv)
    }

    fn refusal(args: &[&str]) -> String {
        let args = parsed(args).expect("the command line parses");
        match run(&args, None) {
            Err(error) => error.to_string(),
            Ok(_) => panic!("the run is not refused"),
        }
    }

    #[test]
    fn floors_masks_and_field_pairs_parse() {
        let floors = parse_floors(Some("HGVSc=30, IMPACT=977")).expect("the floors parse");
        assert_eq!(
            floors.into_iter().collect::<Vec<_>>(),
            [("HGVSc".to_string(), 30), ("IMPACT".to_string(), 977)]
        );
        assert_eq!(parse_floors(None), Ok(IndexMap::new()));
        assert!(parse_floors(Some("HGVSc=many")).is_err());
        assert!(parse_floors(Some("HGVSc")).is_err());
        let masks: Vec<String> = parse_masks(Some("a, b,,c")).into_iter().collect();
        assert_eq!(masks, ["a", "b", "c"]);
        assert!(parse_masks(None).is_empty());
        let pairs = parse_field_pairs(Some(&["HGVSp=p.Arg97%3D=>p.Arg97=".to_string()]))
            .expect("the pair parses");
        let pair = ("p.Arg97%3D".to_string(), "p.Arg97=".to_string());
        assert_eq!(
            pairs.into_iter().collect::<Vec<_>>(),
            [("HGVSp".to_string(), BTreeSet::from([pair]))]
        );
        assert!(parse_field_pairs(Some(&["HGVSp=a:b".to_string()])).is_err());
    }

    #[test]
    fn field_pairs_take_several_items_and_an_equals_sign_inside_a_value() {
        let spec = "IMPACT=MODIFIER=>LOW;Extra=IMPACT=MODIFIER=>IMPACT=LOW".to_string();
        let pairs = parse_field_pairs(Some(&[spec])).expect("the pairs parse");
        let pair = |expected: &str, actual: &str| (expected.to_string(), actual.to_string());
        assert_eq!(
            pairs.into_iter().collect::<Vec<_>>(),
            [
                (
                    "IMPACT".to_string(),
                    BTreeSet::from([pair("MODIFIER", "LOW")])
                ),
                (
                    "Extra".to_string(),
                    BTreeSet::from([pair("IMPACT=MODIFIER", "IMPACT=LOW")])
                ),
            ]
        );
    }

    #[test]
    fn floors_and_field_pairs_keep_the_blanks_and_repeats_their_readers_state() {
        let floors = parse_floors(Some("IMPACT=1,,Gene=0007, SYMBOL=\u{663} ,ZZZ=5,IMPACT=2"))
            .expect("the floors parse");
        assert_eq!(
            floors.into_iter().collect::<Vec<_>>(),
            [
                ("IMPACT".to_string(), 2),
                ("Gene".to_string(), 7),
                ("SYMBOL".to_string(), 3),
                ("ZZZ".to_string(), 5)
            ]
        );
        let message = |spec: &str| parse_floors(Some(spec)).expect_err("refused").to_string();
        let expects = "--floors expects COL=MIN with a non-negative integer, got";
        assert_eq!(message("IMPACT = 1"), format!("{expects} 'IMPACT = 1'"));
        assert_eq!(message(" =3 "), format!("{expects} '=3'"));
        assert_eq!(message("A=1=2"), format!("{expects} 'A=1=2'"));
        assert_eq!(message("A=-1"), format!("{expects} 'A=-1'"));
        assert_eq!(parse_floors(Some(",")), Ok(IndexMap::new()));
        let specs = [
            "IMPACT=HIGH=>LOW;;IMPACT=A=>B; B==>".to_string(),
            "SYMBOL=a=>b".to_string(),
        ];
        let pairs = parse_field_pairs(Some(&specs)).expect("the pairs parse");
        let columns: Vec<&str> = pairs.keys().map(String::as_str).collect();
        assert_eq!(columns, ["IMPACT", " B", "SYMBOL"]);
        assert_eq!(pairs["IMPACT"].len(), 2);
        assert!(pairs[" B"].contains(&(String::new(), String::new())));
        let message = |spec: &str| {
            parse_field_pairs(Some(&[spec.to_string()]))
                .expect_err("refused")
                .to_string()
        };
        let expects = "--field-pairs expects COL=EXPECTED=>ACTUAL, got";
        assert_eq!(message("A=>y"), format!("{expects} 'A=>y'"));
        assert_eq!(message("=x=>y"), format!("{expects} '=x=>y'"));
        assert_eq!(message("ok=a=>b; it's"), format!("{expects} \" it's\""));
    }

    #[test]
    fn a_minimum_the_interpreter_cannot_convert_ends_the_run() {
        let uncaught = |spec: &'static str| {
            let ended = std::panic::catch_unwind(|| parse_floors(Some(spec)));
            let payload = ended.expect_err("the run ends");
            let crate::error::Uncaught(what) = payload
                .downcast_ref::<crate::error::Uncaught>()
                .expect("a run that ends on an input no reading expects");
            what.clone()
        };
        assert_eq!(
            uncaught("IMPACT=\u{b2}"),
            "invalid literal for int() with base 10: '\u{b2}'"
        );
        assert_eq!(
            uncaught("A=9223372036854775808"),
            "--floors A=9223372036854775808: a minimum of 2**63 or more is not read"
        );
        // A later item that is refused is refused before the minimum ends the run, and
        // a later minimum of the same column takes its place.
        assert!(parse_floors(Some("A=9223372036854775808,B")).is_err());
        let floors = parse_floors(Some("A=9223372036854775808,B=1,A=2")).expect("the floors parse");
        assert_eq!(
            floors.into_iter().collect::<Vec<_>>(),
            [("A".to_string(), 2), ("B".to_string(), 1)]
        );
        assert_eq!(
            parse_floors(Some("A=9223372036854775807")).map(|floors| floors["A"]),
            Ok(i64::MAX)
        );
    }

    #[test]
    fn every_option_reaches_its_field() {
        let args = parsed(&[
            "--self-parity",
            "--format",
            "tab",
            "--expected",
            "e",
            "--actual",
            "a",
            "--input",
            "i",
            "--manifest",
            "m",
            "--run",
            "r",
            "--reference-sequence",
            "absent",
            "--masks",
            "x,y",
            "--reference-release",
            "116.2",
            "--keep-only",
            "regulatory_feature_rows",
            "--floors",
            "A=1",
            "--field-pairs",
            "p",
            "--field-pairs",
            "q",
            "--nearest-sites",
            "n",
            "--json-cache",
            "c",
            "--reference-fasta",
            "f",
            "--gwas-catalog",
            "g",
            "--revel-table",
            "t",
            "--revel-no-match",
            "--spliceai-table",
            "s1",
            "--spliceai-table",
            "s2",
            "--utrannotator-sweep",
            "u",
            "--default",
            "d",
            "--tab",
            "b",
            "--vcf",
            "v",
            "--json",
            "j",
            "--parquet",
            "q",
            "--expected-skipped",
            "es",
            "--actual-skipped",
            "as",
            "--unordered",
            "--report",
            "o",
            "--examples",
            "7",
        ])
        .expect("the command line parses");
        let path = |text: &str| Some(PathBuf::from(text));
        let text = |text: &str| Some(text.to_string());
        assert!(args.self_parity && !args.warnings && args.revel_no_match && args.unordered);
        assert_eq!(args.format, text("tab"));
        assert_eq!(
            (args.expected, args.actual, args.input),
            (path("e"), path("a"), path("i"))
        );
        assert_eq!((args.manifest, args.run), (path("m"), text("r")));
        assert_eq!(args.reference_sequence, text("absent"));
        assert_eq!(
            (args.masks, args.keep_only),
            (text("x,y"), text("regulatory_feature_rows"))
        );
        assert_eq!(args.reference_release, text("116.2"));
        assert_eq!(args.floors, text("A=1"));
        assert_eq!(
            args.field_pairs,
            Some(vec!["p".to_string(), "q".to_string()])
        );
        assert_eq!(
            (args.nearest_sites, args.json_cache),
            (path("n"), path("c"))
        );
        assert_eq!(
            (args.reference_fasta, args.gwas_catalog),
            (path("f"), path("g"))
        );
        assert_eq!(args.revel_table, path("t"));
        assert_eq!(
            args.spliceai_table,
            Some(vec![PathBuf::from("s1"), PathBuf::from("s2")])
        );
        assert_eq!(args.utrannotator_sweep, Some(vec![PathBuf::from("u")]));
        assert_eq!(
            (args.default, args.tab, args.vcf),
            (path("d"), path("b"), path("v"))
        );
        assert_eq!((args.json, args.parquet), (path("j"), path("q")));
        assert_eq!(
            (args.expected_skipped, args.actual_skipped),
            (path("es"), path("as"))
        );
        assert_eq!((args.report, args.examples), (path("o"), 7));
        let args = parsed(&["--warnings"]).expect("the command line parses");
        assert!(args.warnings && !args.self_parity);
        assert_eq!(
            (args.examples, args.format, args.field_pairs),
            (EXAMPLES_PER_GROUP, None, None)
        );
        assert_eq!(args.reference_release, None);
        let parser = build_parser();
        assert_eq!(parser.arguments().len(), 33);
        assert!(parser.format_help().starts_with(parser.format_usage()));
        assert!(parser
            .format_usage()
            .starts_with("usage: vep-formats-compare [-h] "));
    }

    #[test]
    fn the_choices_and_the_class_lists_are_the_ones_the_command_line_reads() {
        assert_eq!(VEP_FORMATS, ["default", "tab", "vcf", "json", "parquet"]);
        assert_eq!(FORMATS[..5], *VEP_FORMATS);
        assert_eq!(FORMATS[5], "recoder");
        for format in FORMATS {
            assert_eq!(
                parsed(&["--format", format]).map(|args| args.format),
                Ok(Some(format.to_string()))
            );
        }
        for setting in [REFERENCE_PRESENT, REFERENCE_ABSENT] {
            assert!(parsed(&["--reference-sequence", setting]).is_ok());
        }
        assert_eq!((REFERENCE_PRESENT, REFERENCE_ABSENT), ("present", "absent"));
        assert_eq!(REFERENCE_RELEASES, ["115.2", "116.2", "116.3"]);
        assert_eq!(REFERENCE_RELEASES.last(), Some(&DEFAULT_REFERENCE_RELEASE));
        for release in REFERENCE_RELEASES {
            for spelt in [
                vec!["--reference-release", release],
                vec!["--reference-r", release],
            ] {
                assert_eq!(
                    parsed(&spelt).map(|args| args.reference_release),
                    Ok(Some(release.to_string()))
                );
            }
            assert_eq!(
                parsed(&[&format!("--reference-release={release}")])
                    .map(|args| args.reference_release),
                Ok(Some(release.to_string()))
            );
        }
        assert_eq!(
            parsed(&["--reference-release", "116.4"]).map(|_| ()),
            Err(Exit::Error(
                "argument --reference-release: invalid choice: '116.4' (choose from '115.2', '116.2', '116.3')"
                    .to_string()
            ))
        );
        assert_eq!(
            parsed(&["--reference", "116.3"]).map(|_| ()),
            Err(Exit::Error(
                "ambiguous option: --reference could match --reference-sequence, --reference-release, --reference-fasta"
                    .to_string()
            ))
        );
        assert_eq!(
            parsed(&["--keep-only", "rows"]).map(|_| ()),
            Err(Exit::Error(
                "argument --keep-only: invalid choice: 'rows' (choose from 'regulatory_feature_rows')"
                    .to_string()
            ))
        );
        let tables = classes::tables();
        assert_eq!(
            tables.intron_gap_classes(),
            ["hgvsp_intron_end_insertion", "hgvsp_intron_start_insertion"]
        );
        assert_eq!(
            tables.record_reference_classes(),
            [
                "splice_term_past_the_record_span",
                "plugin_value_follows_splice_term_past_the_record_span",
                "stop_codon_overlap_insertion_range",
                "plugin_value_follows_stop_codon_term_pair",
                "stop_codon_completed_by_polyadenylation",
                "stop_codon_overlap_without_stop_codon",
                "stop_codon_behind_last_residue_without_stop_codon",
                "start_codon_overlap_without_start_codon",
                "refseq_edited_model_equal_span_insertion_term",
                "hgvsp_frameshift_count_at_the_first_residue"
            ]
        );
        assert_eq!(tables.allele_reading_classes().len(), 12);
    }

    #[test]
    fn a_class_the_reference_release_corrects_leaves_the_named_set() {
        let corrected_by_116_3 = [
            classes::FREQUENCY_FILTER_ALLELE_KEY_ORDER,
            classes::FREQUENCY_FILTER_ZERO_VALUE,
        ];
        let named = || {
            parse_masks(Some(
                "order_only,frequency_filter_zero_value,frequency_filter_allele_key_order",
            ))
        };
        for release in ["115.2", "116.2"] {
            assert_eq!(masks_under(named(), release), (named(), BTreeMap::new()));
        }
        let (applied, corrected) = masks_under(named(), "116.3");
        assert_eq!(applied.into_iter().collect::<Vec<_>>(), ["order_only"]);
        assert_eq!(
            corrected.into_iter().collect::<Vec<_>>(),
            corrected_by_116_3.map(|klass| (klass.to_string(), "116.3".to_string()))
        );
        let others = parse_masks(Some("order_only,regulatory_feature_rows"));
        assert_eq!(
            masks_under(others.clone(), "116.3"),
            (others, BTreeMap::new())
        );
        // A release is read by its two numbers, and a text that is no release corrects
        // nothing.
        for (release, corrects) in [
            ("116.10", true),
            ("117.0", true),
            ("116.03", true),
            ("115.9", false),
            ("116", false),
            ("116.3.1", false),
            ("+116.3", false),
            ("", false),
        ] {
            assert_eq!(
                classes::corrected_under(classes::FREQUENCY_FILTER_ZERO_VALUE, release),
                corrects.then_some("116.3"),
                "{release}"
            );
        }
        assert_eq!(classes::corrected_under("order_only", "116.3"), None);
        assert_eq!(
            refusal(&["--masks", "order_only"]).as_str(),
            "--format, --expected and --actual are required"
        );
        let mut args = parsed(&["--warnings"]).expect("the command line parses");
        args.reference_release = Some("116.4".to_string());
        assert_eq!(
            run(&args, None).map_err(|error| error.to_string()),
            Err("--reference-release expects one of 115.2, 116.2, 116.3, got '116.4'".to_string())
        );
    }

    #[test]
    fn a_mode_refuses_what_it_lacks_before_it_reads_anything() {
        let keep_only = "--keep-only is available only in --format comparison mode";
        assert_eq!(
            refusal(&["--keep-only", "regulatory_feature_rows", "--warnings"]),
            keep_only
        );
        assert_eq!(
            refusal(&["--keep-only", "regulatory_feature_rows", "--self-parity"]),
            keep_only
        );
        assert_eq!(
            refusal(&["--warnings", "--expected", "w.txt"]),
            "--warnings needs --expected and --actual"
        );
        assert_eq!(
            refusal(&["--self-parity"]),
            "--self-parity needs --default, --tab, --vcf, --json"
        );
        assert_eq!(
            refusal(&["--self-parity", "--default", "d", "--vcf", "v"]),
            "--self-parity needs --tab, --json"
        );
        let required = "--format, --expected and --actual are required";
        assert_eq!(refusal(&[]), required);
        assert_eq!(refusal(&["--expected", "e", "--actual", "a"]), required);
        assert_eq!(
            refusal(&["--format", "default", "--expected", "e"]),
            required
        );
    }

    /// The class of a stop named on an incomplete coding sequence reads the cache
    /// `--json-cache` names: a row whose FLAGS and terms show nothing stays a difference
    /// without the option and is set aside once with it, and a cache without a transcripts
    /// directory is refused.
    #[test]
    fn a_stop_named_on_an_incomplete_coding_sequence_reads_the_cache_the_command_line_names() {
        use crate::cache::testing::{plant, scratch};
        let klass = classes::HGVSP_STOP_NAMED_ON_INCOMPLETE_CODING_SEQUENCE;
        let root = scratch("partial-codon-command-line");
        plant(
            &root,
            "cache/transcripts/3/1-5.json",
            r#"[{"stable_id": "ENST1", "variation_effect_feature_cache": {"translateable_seq": "NNNN"}}]"#,
        );
        let columns = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\t\
                       Consequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\t\
                       Codons\tExisting_variation\tExtra";
        let row = "r1\t3:7-8\t-\tENSG1\tENST1\tTranscript\tframeshift_variant\t-\t-\t-\t-\t-\t-";
        let output = |extra: &str| format!("{columns}\n{row}\t{extra}\n");
        plant(
            &root,
            "reference.txt",
            &output("IMPACT=HIGH;HGVSp=ENSP1.2:p.Ter612%3D"),
        );
        plant(&root, "candidate.txt", &output("IMPACT=HIGH"));
        let text = |name: &str| path_text(&root.join(name));
        let (reference, candidate) = (text("reference.txt"), text("candidate.txt"));
        // The verdict and the rows set aside under the class, or the refusal.
        let compared = |cache: &[&str]| {
            let mut line = vec!["--format", "default", "--expected", reference.as_str()];
            line.extend(["--actual", candidate.as_str(), "--masks", klass]);
            line.extend(cache);
            let args = parsed(&line).expect("the command line parses");
            run(&args, None)
                .map(|report| {
                    let applied = report.get("masks").and_then(|masks| masks.get("applied"));
                    (
                        report.get("verdict").map(pyjson::py_str),
                        applied.and_then(|applied| applied.get(klass)).cloned(),
                    )
                })
                .map_err(|refused| refused.to_string())
        };
        assert_eq!(compared(&[]), Ok((Some("red".to_string()), None)));
        assert_eq!(
            compared(&["--json-cache", &text("cache")]),
            Ok((Some("green".to_string()), Some(Value::Int(1))))
        );
        assert_eq!(
            compared(&["--json-cache", &text("nowhere")]),
            Err(format!(
                "--json-cache: {} has no transcripts/ directory",
                text("nowhere")
            ))
        );
        std::fs::remove_dir_all(root).unwrap();
    }
}
