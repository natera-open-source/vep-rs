// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `vep recode`: a variant translated between its identifiers and notations,
//! as Ensembl VEP's `variant_recoder` translates it (`VariantRecoder.pm` at
//! release/116.3, whose line numbers the citations below name), from the cache
//! alone.
//!
//! The tool is an annotation run under a fixed configuration (96-150):
//! co-located variants with failed ones included, HGVSg on the chromosome's
//! accession, no statistics, quiet, JSON output and an up- and downstream
//! distance of 0, plus the HGVS, SPDI and GA4GH VRS flags the fields named
//! need; `--mane_select` switches on all three HGVS notations (141-146). The
//! run's JSON objects are then regrouped (265-490): one object per input line,
//! keyed by the input's alleles (`variant_allele`, 308-312); each allele's
//! object holds `input` and, for each field named, the unique values found
//! anywhere in that allele's consequence entries (`find_in_ref`, `Utils.pm`
//! 337-368); an `spdi` names its sequence by the chromosome's `NC_` accession
//! when the cache's synonyms carry one, as the recoder's database slices name
//! it (`VariationFeature.pm:2150-2151`). `vcf_string` is the record's VCF form
//! built by [`to_vcf_record`] for every kind of input (a VCF input's own
//! columns are kept only by a `VCFVariationFeature`, which no input parser
//! makes) and matched to its allele through the minimised and VCF-converted
//! allele strings (318-375, 502-612); `id` is the co-located variants' names
//! per allele, the COSMIC and HGMD ones annexed to every allele (377-421), the
//! co-located variants gathered as a database-backed run gathers them
//! (`Lookups::colocated`); `var_synonyms` is the cache's synonyms of the
//! record's own variant, else of the last co-located one (425-437), one
//! `source: id,id` string per source;
//! `mane_select` is the hgvsg, hgvsc and hgvsp triple of each MANE Select
//! transcript's entry, once per input (439-463); `ga4gh_vrs` is one GA4GH VRS
//! Allele object per `spdi` (474-482, `ga4gh_vrs_from_spdi`). An input whose
//! record has no consequence entry is `null` (its object is never created),
//! and the warnings raised while an input is read or annotated are attached to
//! the next object written, as the recoder attaches the warnings of each
//! buffer fill to its first output line (485-486); a run that writes nothing
//! is one object carrying every warning, an empty list when there was none
//! (`variant_recoder:116`).
//!
//! The variant's own synonyms and its co-located variants' are the cache's
//! `var_synonyms` column; Ensembl VEP reads its database for them (`:427`,
//! `OutputFactory.pm:955-957`).
//!
//! A run that completes writes one `INFO:` line to standard error (`notice`),
//! naming what the recoder's database gives and the cache read does not.
//! Before it stand the lines the annotation run, being quiet, keeps for the
//! tool (`Runner::holding_own_notices`), in the order `vep` prints them: one
//! for each part the cache's record names as left out that the run asks for
//! or reads (`cache_left_out::notices`); one
//! for each FASTA file passed over for want of its index; one for each of
//! `_source_cache` and `dbID` that a `--transcript_filter` names over a cache
//! whose transcripts lack it; one naming each chromosome on which an entry
//! the run writes, as a field or in a `mane_select` triple, is without an
//! `hgvsg` or an `hgvsc` that takes the sequence the FASTA does not hold to
//! write (`Fields::hgvsg_rows`, `Fields::hgvsc_rows`); and,
//! in a run that writes `hgvsc`, one counting the
//! insertion and deletion alleles with a written `hgvsc` on a transcript
//! edited away from the genome (`Fields::hgvsc_rows`). Between those lines
//! and the tool's own stands the tool's line of a cache without variation
//! data, in a run whose fields read known variants
//! (`Lookups::no_variation_data_notice`).

use std::collections::{BTreeSet, HashMap};
use std::ffi::OsString;
use std::fs;
use std::io::{self, BufRead, BufReader, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::time::{SystemTime, UNIX_EPOCH};

use indexmap::IndexMap;
use serde_json::ser::{Formatter, PrettyFormatter};
use serde_json::{Map, Value};
use vep_core::codon::reverse_complement;
use vep_core::coordinate::Strand;
use vep_core::synonyms::ChromosomeSynonyms;
use vep_core::variant::{ColocatedVariant, InputVariant};
use vep_core::variation::CachedVariation;
use vep_io::input::ReferenceSequence;
use vep_io::output::vcf_output::to_vcf_record;

use crate::argv;
use crate::cache_left_out;
use crate::cache_resolve;
use crate::config::{CliError, Config, ParseFailure};
use crate::recode_cli::RecodeArgs;
use crate::runner::{self, NotationRows, Runner};
use crate::variation_matcher::{find_colocated_variants, ColocationRules};
use crate::variation_source::VariationSource;
use crate::warnings::{self, Target};

const COMPONENT: &str = "recode";

/// The default fields (111).
const DEFAULT_FIELDS: &str = "id,hgvsg,hgvsc,hgvsp,spdi";

/// The recoder's own buffer size (96-109 set `buffer_size` to 1): each input
/// line is read, annotated and written on its own, so the warnings raised
/// while a line is read attach to the next object written.
const DEFAULT_BUFFER_SIZE: usize = 1;

/// How far before a record the cache's variation rows are read when its
/// co-located variants and its own row are gathered: a row overlapping the
/// record starts at most its own span before it, and no span in the human
/// variation data reaches 5,000 bases.
const COLOCATION_LOOKBACK: u64 = 5_000;

/// The span of the regions a database-backed run reads the known variants of,
/// Ensembl VEP's `cache_region_size`, which the recoder takes no flag for.
const REGION_SIZE: u64 = 1_000_000;

fn error(message: impl Into<String>) -> CliError {
    CliError::new(COMPONENT, message)
}

/// The refusal of a `--warning_file` that cannot be created or written.
fn warnings_file_unwritable(path: &str) -> CliError {
    error(format!("Could not write to warnings file {path}"))
}

/// Runs the tool for `args`, writing the JSON array to `out`.
pub fn run(args: RecodeArgs, out: &mut dyn Write) -> Result<(), CliError> {
    if let Some(refused) = args.refused_flag() {
        return Err(refused);
    }
    // variant_recoder:104-106: the check reads the species as given, so only
    // a species named on the command line can fail it.
    if args.mane_select
        && args
            .species
            .as_deref()
            .is_some_and(|species| !species.contains("homo_sapiens") && !species.contains("human"))
    {
        return Err(error("MANE Select is only available for human"));
    }
    // VariantRecoder.pm:98 forces `merged` on, so `--refseq` meets Config.pm:964's
    // check of the `refseq => [gencode_basic merged]` pair (:621) on every run,
    // naming `gencode_basic` first where it is given.
    if args.refseq {
        let partner = if args.gencode_basic {
            "gencode_basic"
        } else {
            "merged"
        };
        return Err(error(format!(
            "Can't use --refseq and --{partner} together"
        )));
    }
    let fields = Fields::settle(&args);
    let input = StagedInput::stage(&args)?;
    let output = TempPath::new("json")?;
    warnings::global().configure(Target::Memory, true, input.input_data.is_some());
    let mut config = inner_config(&args, &fields, &input, &output.path, args.merged, None)?;
    // VariantRecoder.pm:98 presets `merged`, so the recoder reads the merged set
    // without the flag (CacheDir.pm:312 suffixes the species directory on every
    // run); a cache that holds the set is read the same way here, and a cache
    // of the Ensembl set alone is read as it is.
    if !args.merged {
        if let Some(cache) = named_merged_cache(&config) {
            let assembly = cache.assembly.as_deref();
            config = inner_config(&args, &fields, &input, &output.path, true, assembly)?;
        }
    }
    let lookups = Lookups::open(&config, &fields, args.failed.unwrap_or(1))?;
    let merged_set = config.merged;
    let run_notices = Runner::new(config)
        .map(|runner| {
            runner
                .holding_own_notices()
                .writing_hgvsc_of(fields.hgvsc_rows())
                .writing_hgvsg_of(fields.hgvsg_rows())
        })
        // The file is written once the annotation has ended, so the run
        // tests the name with its own files, before its first record.
        .map(|runner| match &args.warning_file {
            Some(path) => runner
                .testing_a_file_of_the_caller(PathBuf::from(path), warnings_file_unwritable(path)),
            None => runner,
        })
        .and_then(|runner| runner.run().map(|()| runner.take_own_notices()))
        .map_err(|e| match e.downcast::<CliError>() {
            Ok(refusal) => refusal,
            Err(other) => error(format!("{other:#}")),
        })?;
    let lines = input.numbered_lines()?;
    let records = read_records(&output.path, &lines)?;
    let collected = warnings::global().take_collected();
    // The file is written as `warning_fh` writes it, once there is a warning.
    if let Some(path) = &args.warning_file {
        if !collected.is_empty() {
            let text: String = collected.iter().map(|(_, text)| text.as_str()).collect();
            fs::write(path, text).map_err(|_| warnings_file_unwritable(path))?;
        }
    }
    let buffer_size = args.buffer_size.unwrap_or(DEFAULT_BUFFER_SIZE).max(1);
    let results = regroup(&records, &lines, &fields, &lookups, &collected, buffer_size)?;
    // The annotation run is quiet, as the recoder's is (107), and prints no
    // notice of its own, so the lines it owes for what it gives less of are
    // written here, before the tool's; a run that was refused or could not
    // write its warnings file has ended above.
    for line in &run_notices {
        crate::stderr::write(format_args!("{line}\n"));
    }
    if let Some(line) = lookups.no_variation_data_notice(&fields) {
        crate::stderr::write(format_args!("{line}\n"));
    }
    crate::stderr::write(format_args!(
        "{}\n",
        notice(merged_set, fields.var_synonyms)
    ));
    write_json(out, &Value::Array(results), args.pretty).map_err(|e| error(e.to_string()))
}

/// The line a completed run writes to standard error, naming what the recoder
/// reads from the Ensembl database and the cache does not give. The preset
/// `lrg` (99) maps a variant between its sequence and an LRG, either way
/// (`Parser.pm:901-944`), through the database alone (`Runner.pm:766`). Where
/// `var_synonyms` is written, the `Archive dbSNP` source is named: the dump
/// leaves it out of the cache's column
/// (`Pipeline/DumpVEP/InitDump.pm:263-266`). A cache read as the Ensembl
/// transcript set (`merged_set` false) is named as that, with the merged cache
/// as the one carrying the RefSeq models the preset `merged` (98) reads.
fn notice(merged_set: bool, var_synonyms: bool) -> String {
    let set = if merged_set {
        ""
    } else {
        "this cache is read as the Ensembl transcript set, and RefSeq notation is written from a \
         merged cache; "
    };
    let database = if var_synonyms {
        "no variant is mapped to an LRG and no Archive dbSNP synonym is written, which need the \
         Ensembl database"
    } else {
        "no variant is mapped to an LRG, which needs the Ensembl database"
    };
    format!("INFO: {set}{database}")
}

/// The keys the fields name, settled as `new` settles them (111-150): the fields
/// list with the appended keys, and which flags the annotation run needs.
struct Fields {
    /// The keys read from each allele's consequence entries (271, read at 468).
    want: BTreeSet<String>,
    /// The keys the entries do not carry, attached per allele (275-287).
    id: bool,
    vcf_string: bool,
    var_synonyms: bool,
    mane_select: bool,
    /// `ga4gh_vrs` is wanted, which wants `spdi` (297-301).
    ga4gh_vrs: bool,
    hgvs: bool,
    hgvsg: bool,
    spdi: bool,
}

impl Fields {
    fn settle(args: &RecodeArgs) -> Fields {
        let mut fields: Vec<String> = args
            .fields
            .as_deref()
            .unwrap_or(DEFAULT_FIELDS)
            .split(',')
            .map(str::to_string)
            .collect();
        // 120-122: the flags switched on by the fields named alone. A field
        // starting with `hgvs` sets the configuration key of its own name, and
        // the output reads `hgvsc` and `hgvsp` (`OutputFactory.pm:1759`, `:1767`),
        // which `hgvs` sets both of (`Config.pm:383-386`); `hgvs_offset` sets a
        // key nothing reads, so it switches no notation on.
        let mut hgvsc = fields
            .iter()
            .any(|f| matches!(f.as_str(), "hgvs" | "hgvsc"));
        let mut hgvs = hgvsc || fields.iter().any(|f| f == "hgvsp");
        let mut hgvsg = fields.iter().any(|f| f == "hgvsg");
        let mut spdi = fields.iter().any(|f| f.starts_with("spdi"));
        // 131-150: the appended keys.
        if args.vcf_string {
            fields.push("vcf_string".to_string());
        }
        if args.var_synonyms && !fields.iter().any(|f| f.contains("var_synonyms")) {
            fields.push("var_synonyms".to_string());
        }
        if args.mane_select {
            fields.push("mane_select".to_string());
            hgvsc = true;
            hgvs = true;
            hgvsg = true;
        }
        if args.ga4gh_vrs {
            fields.push("ga4gh_vrs".to_string());
        }
        let mut want: BTreeSet<String> = fields.into_iter().collect();
        // The offset is written beside the transcript notation alone: the
        // protein branch (`OutputFactory.pm:1767-1776`) reads an offset only
        // `hgvs_transcript` sets (`TranscriptVariationAllele.pm:1447`), since
        // `hgvs_protein` drops the shift it made (`:1796`), so without `hgvsc`
        // switched on no entry carries the key.
        if !hgvsc {
            want.remove("hgvs_offset");
        }
        // 275-301: the keys with no entry of their own leave the wanted set;
        // ga4gh_vrs stays in it and wants spdi beside it.
        let id = want.remove("id");
        let vcf_string = want.remove("vcf_string");
        let var_synonyms = want.remove("var_synonyms");
        let mane_select = want.remove("mane_select");
        let ga4gh_vrs = want.contains("ga4gh_vrs");
        if ga4gh_vrs {
            want.insert("spdi".to_string());
            spdi = true;
        }
        Fields {
            want,
            id,
            vcf_string,
            var_synonyms,
            mane_select,
            ga4gh_vrs,
            hgvs,
            hgvsg,
            spdi,
        }
    }

    /// The consequence entries whose `hgvsc` the output writes: every entry
    /// where the keys read hold `hgvsc` (271, read at 468), the MANE Select
    /// entries alone where `mane_select` writes their triples without it
    /// (439-463), and none otherwise. A field named `hgvs` switches the
    /// notation on and reads a key no entry carries.
    fn hgvsc_rows(&self) -> NotationRows {
        self.rows_of("hgvsc")
    }

    /// The consequence entries whose `hgvsg` the output writes, settled as
    /// [`Self::hgvsc_rows`] settles those of `hgvsc`: every entry where the
    /// keys read hold `hgvsg`, the MANE Select entries alone where
    /// `mane_select` writes their triples without it, and none otherwise.
    fn hgvsg_rows(&self) -> NotationRows {
        self.rows_of("hgvsg")
    }

    /// The consequence entries whose `key` the output writes, `key` being a
    /// notation of the `mane_select` triple.
    fn rows_of(&self, key: &str) -> NotationRows {
        if self.want.contains(key) {
            NotationRows::Every
        } else if self.mane_select {
            NotationRows::ManeSelect
        } else {
            NotationRows::None
        }
    }
}

/// The input the annotation run reads: the file or string as given, or
/// standard input copied to a file so the run and the line count read the
/// same bytes.
struct StagedInput {
    input_file: Option<String>,
    input_data: Option<String>,
    _stdin_copy: Option<TempPath>,
}

impl StagedInput {
    fn stage(args: &RecodeArgs) -> Result<StagedInput, CliError> {
        if let Some(data) = args.input_data.as_deref().filter(|data| argv::truthy(data)) {
            return Ok(StagedInput {
                input_file: None,
                input_data: Some(data.to_string()),
                _stdin_copy: None,
            });
        }
        match args.input_file.as_deref() {
            Some(file) if !file.eq_ignore_ascii_case("STDIN") => Ok(StagedInput {
                input_file: Some(file.to_string()),
                input_data: None,
                _stdin_copy: None,
            }),
            _ => {
                let copy = TempPath::new("input")?;
                let mut bytes = Vec::new();
                io::stdin()
                    .read_to_end(&mut bytes)
                    .map_err(|e| error(format!("Could not read standard input: {e}")))?;
                fs::write(&copy.path, bytes).map_err(|e| error(e.to_string()))?;
                Ok(StagedInput {
                    input_file: Some(copy.path.to_string_lossy().into_owned()),
                    input_data: None,
                    _stdin_copy: Some(copy),
                })
            }
        }
    }

    /// The input's lines, numbered as the parsers number them: every line but
    /// a header line counts, an empty one included (`Parser.pm:273-276`).
    fn numbered_lines(&self) -> Result<Vec<(u64, String)>, CliError> {
        let text = match (&self.input_data, &self.input_file) {
            (Some(data), _) => data.clone(),
            (None, Some(file)) => read_text(file)?,
            (None, None) => String::new(),
        };
        let mut number = 0;
        Ok(text
            .lines()
            .filter(|line| !line.starts_with('#'))
            .map(|line| {
                number += 1;
                (number, line.to_string())
            })
            .collect())
    }
}

/// The text of an input file, a gzip or BGZF file decompressed.
fn read_text(path: &str) -> Result<String, CliError> {
    let file = fs::File::open(path)
        .map_err(|_| error(format!("Could not read from input file {path}")))?;
    let mut reader = BufReader::new(file);
    let gzip = reader
        .fill_buf()
        .map_err(|e| error(e.to_string()))?
        .starts_with(&[0x1f, 0x8b]);
    let mut text = String::new();
    let read = if gzip {
        flate2::read::MultiGzDecoder::new(reader).read_to_string(&mut text)
    } else {
        reader.read_to_string(&mut text)
    };
    read.map_err(|e| error(format!("Could not read from input file {path}: {e}")))?;
    Ok(text)
}

/// A file under the temporary directory, removed when dropped.
struct TempPath {
    path: PathBuf,
}

impl TempPath {
    fn new(kind: &str) -> Result<TempPath, CliError> {
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .map(|d| d.as_nanos())
            .unwrap_or(0);
        let path =
            std::env::temp_dir().join(format!("vep-recode-{}-{nanos}.{kind}", std::process::id()));
        if path.exists() {
            return Err(error(format!(
                "temporary file {} already exists",
                path.display()
            )));
        }
        Ok(TempPath { path })
    }
}

impl Drop for TempPath {
    fn drop(&mut self) {
        let _ = fs::remove_file(&self.path);
    }
}

/// The annotation run's configuration, reading the merged transcript set when
/// `merged`, under `assembly` when the command line names none.
fn inner_config(
    args: &RecodeArgs,
    fields: &Fields,
    input: &StagedInput,
    output: &Path,
    merged: bool,
    assembly: Option<&str>,
) -> Result<Config, CliError> {
    match Config::from_argv(inner_argv(args, fields, input, output, merged, assembly)) {
        Ok(config) => Ok(Config {
            spdi_on_the_cache_name: true,
            ..config
        }),
        Err(ParseFailure::Refused(refusal)) => Err(refusal),
        // The tool's own parser holds a value to the list the annotation run's
        // parser declares, and each value is handed on joined to its option, so
        // this arm is for arguments built in code, which can carry a value
        // outside a list: the refusal is the first line of the parser's text,
        // which names the value.
        Err(ParseFailure::Usage(usage)) => {
            let text = usage.to_string();
            let named = text.lines().next().unwrap_or_default();
            Err(error(format!(
                "the annotation run refused its arguments: {}",
                named.strip_prefix("error: ").unwrap_or(named)
            )))
        }
    }
}

/// The merged cache a run names without `--merged`, as [`named_merged_cache`]
/// finds it.
struct MergedCache {
    /// The assembly of the Ensembl set's cache beside it where that cache's
    /// directory has the name the merged run looks up under that assembly,
    /// `<cache_version>_<assembly>` with the assembly's canonical name: the
    /// merged run is given it so that it reads the cache of that name where
    /// the `_merged` directory holds several. Under any other name the merged
    /// run is given none and resolves the directory as `--merged` does.
    assembly: Option<String>,
}

/// What the merged run is given when the cache `config` names holds the merged
/// transcript set, and `None` when it does not. A `--json_cache` directory
/// holds it when its `info.json` declares the set. A cache root holds it when
/// its `<species>_merged` directory, which `CacheDir.pm:312` resolves for the
/// recoder on every run, has an entry named as the cache `config` resolves,
/// with `info.json` or `transcripts/` and declaring no other set or assembly
/// (Ensembl VEP's own Storable cache has neither file), or when the root has
/// that directory and no `<species>` one, which `--merged` is left to resolve.
fn named_merged_cache(config: &Config) -> Option<MergedCache> {
    let info_of = |dir: &Path| {
        let text = fs::read_to_string(dir.join("info.json")).ok()?;
        serde_json::from_str::<Value>(&text).ok()
    };
    if let Some(dir) = &config.json_cache {
        let info = info_of(Path::new(dir))?;
        return (info["source_type"] == "merged").then_some(MergedCache { assembly: None });
    }
    let species_dir = format!("{}/{}", config.dir_cache, config.species);
    let merged_dir = PathBuf::from(format!("{species_dir}_merged"));
    match cache_resolve::resolve(&mut config.clone()) {
        Ok(cache) => {
            let leaf = Path::new(&cache.dir).file_name()?.to_str()?;
            let merged_cache = merged_dir.join(leaf);
            if !merged_cache.join("info.json").is_file()
                && !merged_cache.join("transcripts").is_dir()
            {
                return None;
            }
            // An absent `info.json`, or one that is not JSON, declares nothing
            // here: the merged run refuses such a cache where it resolves it.
            let declares_another = info_of(&merged_cache).is_some_and(|info| {
                let differs = |key: &str, own: Option<&str>| match (info[key].as_str(), own) {
                    (Some(declared), Some(own)) => !declared.is_empty() && declared != own,
                    _ => false,
                };
                differs("source_type", Some("merged"))
                    || differs("assembly", cache.assembly.as_deref())
            });
            if declares_another {
                return None;
            }
            let assembly = cache.assembly.filter(|assembly| {
                let looked_up = cache_resolve::canonical_assembly(assembly);
                leaf == format!("{}_{looked_up}", config.cache_version)
            });
            Some(MergedCache { assembly })
        }
        Err(_) => (merged_dir.is_dir() && !Path::new(&species_dir).is_dir())
            .then_some(MergedCache { assembly: None }),
    }
}

/// The annotation run's command line: the forced configuration (96-129), the
/// flags the fields need, `--merged` when `merged`, `--assembly` as given or
/// else as `assembly`, and the other options forwarded as given. Each value is
/// handed on joined to its option (`--option=value`), so the run's parser
/// reads it as that option's value whatever its first characters are: a
/// transcript filter, an input or a path that begins with a hyphen names no
/// option of the run.
fn inner_argv(
    args: &RecodeArgs,
    fields: &Fields,
    input: &StagedInput,
    output: &Path,
    merged: bool,
    assembly: Option<&str>,
) -> Vec<OsString> {
    let mut argv: Vec<String> = vec![
        "vep".into(),
        "--json".into(),
        "--hgvsg_use_accession".into(),
        "--no_stats".into(),
        "--quiet".into(),
        "--distance=0".into(),
        "--force_overwrite".into(),
        format!("--output_file={}", output.to_string_lossy()),
    ];
    let mut push_switch = |flag: &str, on: bool| {
        if on {
            argv.push(format!("--{flag}"));
        }
    };
    push_switch("hgvs", fields.hgvs);
    push_switch("hgvsg", fields.hgvsg);
    push_switch("spdi", fields.spdi);
    push_switch("mane_select", fields.mane_select);
    push_switch("grch37", args.grch37);
    push_switch("merged", merged);
    push_switch("all_refseq", args.all_refseq);
    push_switch("gencode_basic", args.gencode_basic);
    push_switch("exclude_predicted", args.exclude_predicted);
    push_switch("minimal", args.minimal);
    push_switch("pick", args.pick);
    push_switch("pick_allele", args.pick_allele);
    push_switch("per_gene", args.per_gene);
    push_switch("pick_allele_gene", args.pick_allele_gene);
    push_switch("use_transcript_ref", args.use_transcript_ref);
    push_switch("use_given_ref", args.use_given_ref);
    let valued: [(&str, Option<&str>); 12] = [
        ("input_file", input.input_file.as_deref()),
        ("input_data", input.input_data.as_deref()),
        ("format", args.format.as_deref()),
        ("delimiter", args.delimiter.as_deref()),
        ("species", args.species.as_deref()),
        ("assembly", args.assembly.as_deref().or(assembly)),
        ("pick_order", args.pick_order.as_deref()),
        ("shift_hgvs", args.shift_hgvs.as_deref()),
        ("synonyms", args.synonyms.as_deref()),
        ("dir_cache", args.dir_cache.as_deref()),
        ("json_cache", args.json_cache.as_deref()),
        ("fasta", args.fasta.as_deref()),
    ];
    for (flag, value) in valued {
        if let Some(value) = value {
            argv.push(format!("--{flag}={value}"));
        }
    }
    if let Some(version) = args.cache_version {
        argv.push(format!("--cache_version={version}"));
    }
    if let Some(size) = args.buffer_size {
        argv.push(format!("--buffer_size={size}"));
    }
    for filter in &args.transcript_filter {
        argv.push(format!("--transcript_filter={filter}"));
    }
    argv.into_iter().map(OsString::from).collect()
}

/// What the regrouping reads beside the run's output: the reference the VCF
/// form's anchor base comes from, the cache's variation files the co-located
/// variants and the synonyms are read from, and the chromosome synonyms the
/// SPDI's accession comes from.
struct Lookups {
    reference: Option<Arc<vep_fasta::IndexedFasta>>,
    variations: Option<VariationSource>,
    /// The cache's directory where the fields read known variants and the
    /// cache holds no variation data.
    no_variation_data_under: Option<String>,
    synonyms: Arc<ChromosomeSynonyms>,
    /// The co-location rules of the run: alleles checked, null alleles kept,
    /// `--failed` as given (the recoder's default includes failed variants).
    rules: ColocationRules,
    /// The variation rows a test stands in for the cache with.
    #[cfg(test)]
    rows: Vec<CachedVariation>,
}

impl Lookups {
    fn open(config: &Config, fields: &Fields, failed: i64) -> Result<Lookups, CliError> {
        let rules = ColocationRules {
            failed,
            ..ColocationRules::default()
        };
        if !fields.vcf_string && !fields.var_synonyms && !fields.spdi && !fields.id {
            return Ok(Lookups {
                reference: None,
                variations: None,
                no_variation_data_under: None,
                synonyms: Arc::new(ChromosomeSynonyms::default()),
                rules,
                #[cfg(test)]
                rows: Vec::new(),
            });
        }
        let mut probe = config.clone();
        let cache = cache_resolve::resolve(&mut probe)?;
        let reference = if fields.vcf_string {
            runner::open_fasta(&probe, &cache)?
        } else {
            None
        };
        let variations = if fields.var_synonyms || fields.id {
            VariationSource::open(
                Path::new(&cache.dir),
                &cache.info,
                probe.old_maf,
                Arc::clone(&cache.manifest),
            )
            .map_err(|e| CliError::new("cache", format!("{e:#}")))?
        } else {
            None
        };
        let reads_known_variants = fields.var_synonyms || fields.id;
        let no_variation_data_under =
            (reads_known_variants && variations.is_none()).then(|| cache.dir.clone());
        Ok(Lookups {
            reference,
            variations,
            no_variation_data_under,
            synonyms: Arc::clone(&cache.synonyms),
            rules,
            #[cfg(test)]
            rows: Vec::new(),
        })
    }

    /// The line of a run whose fields read known variants over a cache
    /// without variation data, naming those fields
    /// (`cache_left_out::no_variation_data_recode_notice`): the tool opens
    /// the variation files itself and its annotation run, which never asks
    /// for known variants, raises no warning there. `None` over a cache that
    /// holds variation data and in a run whose fields read none.
    fn no_variation_data_notice(&self, fields: &Fields) -> Option<String> {
        let cache_dir = self.no_variation_data_under.as_deref()?;
        let read: Vec<&str> = [("id", fields.id), ("var_synonyms", fields.var_synonyms)]
            .into_iter()
            .filter_map(|(field, read)| read.then_some(field))
            .collect();
        Some(cache_left_out::no_variation_data_recode_notice(
            &read, cache_dir,
        ))
    }

    /// The cache's variation rows a database-backed Ensembl VEP reads for the
    /// record: the rows starting in the regions of its span widened by one base
    /// (`AnnotationSource.pm:189`, `Database/Variation.pm:167-179`), read from
    /// [`COLOCATION_LOOKBACK`] bases before the record at most, in the cache's
    /// order.
    fn rows_near(&self, head: &RecordHead) -> Result<Vec<CachedVariation>, CliError> {
        let Some(source) = &self.variations else {
            #[cfg(test)]
            return Ok(self.rows.clone());
            #[cfg(not(test))]
            return Ok(Vec::new());
        };
        let (low, high) = (head.start.min(head.end), head.start.max(head.end));
        let (widened_low, widened_high) = (low.saturating_sub(1), high + 1);
        let region_start = widened_low.saturating_sub(1) / REGION_SIZE * REGION_SIZE + 1;
        let first = region_start.max(low.saturating_sub(COLOCATION_LOOKBACK));
        source
            .records_starting_within(&head.chr, first, widened_high)
            .map_err(|e| CliError::new("cache", format!("{e:#}")))
    }

    /// The known variants co-located with the record as a database-backed run
    /// finds them (`AnnotationType/Variation.pm:91-109`): the rows of
    /// [`Lookups::rows_near`] whose span overlaps the record's as the input
    /// buffer tests it (`InputBuffer.pm:286-324`, an insertion overlapping only
    /// a row that spans both of its flanking bases), each name once, kept when
    /// `compare_existing` matches an allele.
    fn colocated(&self, head: &RecordHead, rows: &[CachedVariation]) -> Vec<ColocatedVariant> {
        let candidates: Vec<CachedVariation> = rows
            .iter()
            .filter(|row| {
                let (row_low, row_high) = (row.start.min(row.end), row.start.max(row.end));
                head.end >= row_low && head.start <= row_high
            })
            .cloned()
            .collect();
        find_colocated_variants(&head.variant(), candidates, &self.rules)
    }

    /// An SPDI with its sequence named by the RefSeq name among the name's
    /// synonyms, an `NC_` accession before one of another prefix, as the
    /// recoder's database slices name it (`VariationFeature.pm:2150-2151`);
    /// unchanged when the name has none or is a later form of that accession.
    fn spdi_on_accession(&self, spdi: &str) -> String {
        let Some((name, rest)) = spdi.split_once(':') else {
            return spdi.to_string();
        };
        match vep_effects::hgvsg::refseq_name(
            self.synonyms.synonyms_of(name).iter().map(String::as_str),
        ) {
            Some(accession) if vep_effects::hgvsg::later_form(accession, name) != Some(name) => {
                format!("{accession}:{rest}")
            }
            _ => spdi.to_string(),
        }
    }
}

/// The `var_synonyms` column of the row of `rows` named `id`, as
/// [`written_synonyms`](vep_core::variant::written_synonyms) reads it; none
/// for an empty name. The row is the one spelled as `id` is, else the first
/// named `id` in another case, the cache's spelling of a name the input gave.
/// The recoder fetches the variant by `v.name = ?` (`VariationAdaptor.pm:544`
/// of ensembl-variation, from 427 and `Parser/ID.pm:162`) over a column that
/// holds each name once (`sql/table.sql:77`): a database that compares names
/// without regard to case holds one row for the name, and one that holds two
/// spellings of it returns the row spelled as asked.
fn synonyms_of(rows: &[CachedVariation], id: &str) -> Option<String> {
    if id.is_empty() {
        return None;
    }
    rows.iter()
        .find(|row| row.variation_name == id)
        .or_else(|| {
            rows.iter()
                .find(|row| row.variation_name.eq_ignore_ascii_case(id))
        })
        .and_then(|row| vep_core::variant::written_synonyms(row.var_synonyms.as_deref()))
        .map(str::to_string)
}

/// The run's FASTA as [`to_vcf_record`] reads it.
struct Reference(Arc<vep_fasta::IndexedFasta>);

impl ReferenceSequence for Reference {
    fn sequence(&self, chr: &str, start: u64, end: u64) -> Option<Vec<u8>> {
        self.0.sequence(chr, start, end)
    }

    fn record_name(&self, chr: &str) -> Option<String> {
        self.0.record_name(chr).map(str::to_string)
    }

    fn sequence_length(&self, chr: &str) -> Option<u64> {
        self.0.sequence_length(chr)
    }
}

/// One record of the run's JSON output with the number of the input line it
/// came from, found by walking the input's lines in order.
struct Record {
    line: Option<u64>,
    value: Map<String, Value>,
}

/// The run's JSON objects, one per input record and in input order, each
/// matched to its line of the input (`lines`, numbered).
fn read_records(output: &Path, lines: &[(u64, String)]) -> Result<Vec<Record>, CliError> {
    let text =
        fs::read_to_string(output).map_err(|e| error(format!("reading the run's output: {e}")))?;
    let mut cursor = 0usize;
    let mut records = Vec::new();
    for json in text.lines().filter(|l| !l.trim().is_empty()) {
        let value: Value = serde_json::from_str(json)
            .map_err(|e| error(format!("reading the run's output: {e}")))?;
        let Value::Object(object) = value else {
            return Err(error("the run's output holds a line that is not an object"));
        };
        let line = object
            .get("input")
            .and_then(Value::as_str)
            .and_then(|input| {
                let found = lines[cursor..].iter().position(|(_, text)| text == input)?;
                let (number, _) = &lines[cursor + found];
                cursor += found + 1;
                Some(*number)
            });
        records.push(Record {
            line,
            value: object,
        });
    }
    Ok(records)
}

/// The output of a record's `vcf_string` field and the per-allele keys the
/// consequence entries do not carry (318-437).
#[derive(Default)]
struct AlleleExtras {
    vcf_string: Option<String>,
    id: Vec<String>,
    var_synonyms: Vec<String>,
}

/// `_get_all_results` (265-490) over the run's records; a variation file that
/// cannot be read is the run's error, never an empty `id`.
fn regroup(
    records: &[Record],
    lines: &[(u64, String)],
    fields: &Fields,
    lookups: &Lookups,
    collected: &[(Option<u64>, String)],
    buffer_size: usize,
) -> Result<Vec<Value>, CliError> {
    let mut results: IndexMap<String, IndexMap<String, IndexMap<String, Value>>> = IndexMap::new();
    let mut order: Vec<String> = Vec::new();
    let mut mane_seen: BTreeSet<String> = BTreeSet::new();
    for record in records {
        let line = &record.value;
        let Some(input) = line
            .get("input")
            .and_then(Value::as_str)
            .map(str::to_string)
        else {
            continue;
        };
        let consequences: Vec<&Value> = line
            .get("transcript_consequences")
            .or_else(|| line.get("intergenic_consequences"))
            .and_then(Value::as_array)
            .map(|entries| entries.iter().collect())
            .unwrap_or_default();
        // 308-312: the entries by variant allele, in the order the alleles appear.
        let mut by_allele: IndexMap<String, Vec<&Value>> = IndexMap::new();
        for entry in &consequences {
            let allele = entry
                .get("variant_allele")
                .map(value_text)
                .unwrap_or_default();
            by_allele.entry(allele).or_default().push(entry);
        }
        let mut extras: IndexMap<String, AlleleExtras> = IndexMap::new();
        let head = RecordHead::of(line);

        if fields.vcf_string {
            let minimised = minimised_alleles(&head.allele_string);
            let vcf_alleles = alleles_in_vcf_form(&head.allele_string, head.strand);
            for (vcf_string, reference_allele, alt) in vcf_strings(&head, lookups) {
                // 328-375: the allele the VCF form's alternate spells.
                let allele = vcf_alleles
                    .get(&format!("{reference_allele}/{alt}"))
                    .cloned()
                    .or_else(|| {
                        allele_of_vcf_pair(&reference_allele, &alt, &minimised, &vcf_alleles)
                    })
                    .or_else(|| {
                        allele_of_vcf_pair(&alt, &reference_allele, &minimised, &vcf_alleles)
                    });
                if let Some(allele) = allele {
                    extras.entry(allele).or_default().vcf_string = Some(vcf_string);
                }
            }
        }

        let rows = if fields.id || fields.var_synonyms {
            lookups.rows_near(&head)?
        } else {
            Vec::new()
        };
        // 379-421: the co-located variants' names per allele.
        let mut last_named: Option<String> = None;
        if fields.id {
            let colocated: Vec<(String, String)> = lookups
                .colocated(&head, &rows)
                .into_iter()
                .map(|cv| (cv.id, cv.allele_string.unwrap_or_default()))
                .collect();
            let mut unassigned: Vec<String> = Vec::new();
            for (cv_id, cv_alleles) in &colocated {
                if cv_alleles.contains("COSMIC") || cv_alleles.contains("HGMD") {
                    unassigned.push(cv_id.clone());
                    continue;
                }
                last_named = Some(cv_id.clone());
                let mut parts = cv_alleles.split('/');
                let cv_reference = parts.next().unwrap_or_default();
                for allele in parts {
                    let key = if by_allele.contains_key(allele) {
                        allele.to_string()
                    } else {
                        let reversed =
                            String::from_utf8_lossy(&reverse_complement(allele.as_bytes()))
                                .into_owned();
                        if by_allele.contains_key(&reversed) {
                            reversed
                        } else {
                            minimised_pair(cv_reference, allele)
                                .split_once('/')
                                .map(|(_, alt)| alt.to_string())
                                .unwrap_or_default()
                        }
                    };
                    extras.entry(key).or_default().id.push(cv_id.clone());
                }
            }
            for allele in by_allele.keys() {
                extras
                    .entry(allele.clone())
                    .or_default()
                    .id
                    .extend(unassigned.iter().cloned());
            }
        }

        // 425-437: the synonyms of the record's own variant (the row named as
        // the input named it), else of the last named co-located one, attached
        // to every allele.
        if fields.var_synonyms {
            let synonyms = synonyms_of(&rows, &head.name)
                .or_else(|| last_named.as_deref().and_then(|id| synonyms_of(&rows, id)))
                .map(|column| synonym_strings(&column))
                .unwrap_or_default();
            for allele in by_allele.keys() {
                extras.entry(allele.clone()).or_default().var_synonyms = synonyms.clone();
            }
        }

        // 439-463: the MANE Select triples, unique per input.
        let mut mane_by_allele: IndexMap<String, Vec<Value>> = IndexMap::new();
        if fields.mane_select {
            for (allele, entries) in &by_allele {
                for entry in entries {
                    let Some(object) = entry.as_object() else {
                        continue;
                    };
                    if !object.keys().any(|k| k.contains("mane_select")) {
                        continue;
                    }
                    let text = |key: &str| object.get(key).map(value_text).unwrap_or_default();
                    let (hgvsg, hgvsc, hgvsp) = (text("hgvsg"), text("hgvsc"), text("hgvsp"));
                    if mane_seen.insert(format!("{input}-{hgvsg}-{hgvsc}")) {
                        let triple = serde_json::json!({
                            "hgvsg": dash_when_empty(hgvsg),
                            "hgvsc": dash_when_empty(hgvsc),
                            "hgvsp": dash_when_empty(hgvsp),
                        });
                        mane_by_allele
                            .entry(allele.clone())
                            .or_default()
                            .push(triple);
                    }
                }
            }
        }

        if !order.contains(&input) {
            order.push(input.clone());
        }
        let by_input = results.entry(input.clone()).or_default();
        for (allele, entries) in &by_allele {
            let object = by_input.entry(allele.clone()).or_insert_with(|| {
                IndexMap::from([("input".to_string(), Value::String(input.clone()))])
            });
            // The SPDI is put on its accession before the merge, so a repeated
            // input's SPDI is one text, which the merge keeps once.
            let mut found: IndexMap<String, Value> = IndexMap::new();
            for entry in entries {
                find_in_ref(entry, &fields.want, &mut found, None);
            }
            if let Some(Value::Array(spdis)) = found.get_mut("spdi") {
                for spdi in spdis.iter_mut() {
                    *spdi = Value::String(lookups.spdi_on_accession(&value_text(spdi)));
                }
            }
            for (key, values) in found {
                if let Value::Array(values) = values {
                    for value in values {
                        push_unique(object, &key, value);
                    }
                }
            }
            if let Some(extra) = extras.get(allele) {
                if let Some(vcf_string) = &extra.vcf_string {
                    push_unique(object, "vcf_string", Value::String(vcf_string.clone()));
                }
                for id in &extra.id {
                    push_unique(object, "id", Value::String(id.clone()));
                }
                for synonym in &extra.var_synonyms {
                    push_unique(object, "var_synonyms", Value::String(synonym.clone()));
                }
            }
            if let Some(triples) = mane_by_allele.get(allele) {
                push_all(object, "mane_select", triples.clone());
            }
        }
        // 474-482: one GA4GH VRS allele per spdi on an `NC_` accession. The push
        // of 479 makes the list before it appends what the SPDI gives, so an
        // allele with an spdi carries the key, an empty list when no spdi is on
        // an `NC_` accession, and an allele without one (476) has no key.
        if fields.ga4gh_vrs {
            for object in by_input.values_mut() {
                let spdis: Vec<String> = object
                    .get("spdi")
                    .and_then(Value::as_array)
                    .map(|spdis| spdis.iter().map(value_text).collect())
                    .unwrap_or_default();
                if spdis.is_empty() {
                    continue;
                }
                let vrs_alleles = spdis
                    .iter()
                    .filter_map(|spdi| ga4gh_vrs_from_spdi(spdi))
                    .collect();
                push_all(object, "ga4gh_vrs", vrs_alleles);
            }
        }
    }

    // 485-486 and variant_recoder:116: the warnings of each buffer fill go to
    // the fill's first output line; a run with no output line carries them all.
    let warnings = attach_warnings(records, lines, collected, buffer_size);
    let mut out: Vec<Value> = order
        .iter()
        .enumerate()
        .map(|(i, input)| {
            let alleles = &results[input];
            let mut object: Map<String, Value> = alleles
                .iter()
                .map(|(allele, fields)| {
                    (
                        allele.clone(),
                        Value::Object(fields.iter().map(|(k, v)| (k.clone(), v.clone())).collect()),
                    )
                })
                .collect();
            if let Some(texts) = warnings.get(&i) {
                object.insert(
                    "warnings".to_string(),
                    Value::Array(texts.iter().map(|t| Value::String(t.clone())).collect()),
                );
            }
            if object.is_empty() {
                Value::Null
            } else {
                Value::Object(object)
            }
        })
        .collect();
    if out.is_empty() {
        let texts: Vec<Value> = collected
            .iter()
            .flat_map(|(_, text)| warning_messages(text))
            .map(Value::String)
            .collect();
        out.push(Value::Object(Map::from_iter([(
            "warnings".to_string(),
            Value::Array(texts),
        )])));
    }
    Ok(out)
}

/// Which output object (by position in `order`, the records' input order with
/// repeats merged) carries each warning: the first record of the buffer fill
/// that read the warning's line (a fill of `buffer_size` records consumes the
/// lines up to its last record's, the final short fill every line after it). A
/// warning without a line number takes the number of the input line it quotes,
/// else goes to the first record. A warning read after the last fill is lost
/// when that fill is full, as the recoder loses it (485-486).
fn attach_warnings(
    records: &[Record],
    lines: &[(u64, String)],
    collected: &[(Option<u64>, String)],
    buffer_size: usize,
) -> HashMap<usize, Vec<String>> {
    let mut inputs: Vec<String> = Vec::new();
    let mut position_of_record: Vec<usize> = Vec::with_capacity(records.len());
    for record in records {
        let input = record
            .value
            .get("input")
            .map(value_text)
            .unwrap_or_default();
        let position = match inputs.iter().position(|i| *i == input) {
            Some(position) => position,
            None => {
                inputs.push(input);
                inputs.len() - 1
            }
        };
        position_of_record.push(position);
    }
    let mut fills: Vec<(usize, Option<u64>)> = Vec::new();
    for (i, chunk) in records.chunks(buffer_size).enumerate() {
        let full = chunk.len() == buffer_size;
        let last_line = chunk.iter().filter_map(|r| r.line).max();
        fills.push((i * buffer_size, if full { last_line } else { None }));
    }
    let mut by_position: HashMap<usize, Vec<String>> = HashMap::new();
    for (line, text) in collected {
        let messages = warning_messages(text);
        if messages.is_empty() {
            continue;
        }
        let line = line.or_else(|| quoted_line_number(text, lines));
        let fill = match line {
            None => fills.first(),
            Some(n) => fills
                .iter()
                .find(|(_, last)| last.is_none_or(|last| n <= last)),
        };
        let Some((first_record, _)) = fill else {
            continue;
        };
        let Some(position) = position_of_record.get(*first_record) else {
            continue;
        };
        by_position.entry(*position).or_default().extend(messages);
    }
    by_position
}

/// The number of the first input line a warning quotes, in single or double
/// quotes or in parentheses (`No variant found with ID 'rs1'`).
fn quoted_line_number(text: &str, lines: &[(u64, String)]) -> Option<u64> {
    lines
        .iter()
        .find(|(_, line)| {
            !line.is_empty()
                && [
                    format!("'{line}'"),
                    format!("\"{line}\""),
                    format!("({line})"),
                ]
                .iter()
                .any(|quoted| text.contains(quoted.as_str()))
        })
        .map(|(number, _)| *number)
}

/// The `msg` texts of a warning as `Runner::warnings` reads them (`Runner.pm:321-357`):
/// a `WARNING:` or `ERROR:` line opens a message, a `MSG:` line extends the
/// open one after `: `, every other line is dropped.
fn warning_messages(text: &str) -> Vec<String> {
    let mut messages: Vec<String> = Vec::new();
    for line in text.lines().filter(|l| !l.is_empty()) {
        let opened = ["WARNING", "ERROR"].iter().find_map(|kind| {
            let rest = line.strip_prefix(kind)?.trim_start();
            let rest = rest.strip_prefix(':')?.trim_start();
            (!rest.is_empty()).then(|| rest.to_string())
        });
        if let Some(message) = opened {
            messages.push(message);
        } else if let Some(rest) = line
            .strip_prefix("MSG")
            .and_then(|r| r.trim_start().strip_prefix(':'))
        {
            let rest = rest.trim_start();
            if rest.is_empty() {
                continue;
            }
            if let Some(last) = messages.last_mut() {
                last.push_str(": ");
                last.push_str(rest);
            }
        }
    }
    messages
}

/// A record's identity and coordinates as its JSON object states them.
struct RecordHead {
    chr: String,
    start: u64,
    end: u64,
    strand: Strand,
    allele_string: String,
    name: String,
}

impl RecordHead {
    fn of(line: &Map<String, Value>) -> RecordHead {
        let start = line.get("start").and_then(Value::as_u64).unwrap_or(0);
        RecordHead {
            chr: line
                .get("seq_region_name")
                .map(value_text)
                .unwrap_or_default(),
            start,
            end: line.get("end").and_then(Value::as_u64).unwrap_or(start),
            strand: line
                .get("strand")
                .and_then(Value::as_i64)
                .and_then(|s| Strand::from_i8(s as i8))
                .unwrap_or(Strand::Forward),
            allele_string: line
                .get("allele_string")
                .map(value_text)
                .unwrap_or_default(),
            name: line.get("id").map(value_text).unwrap_or_default(),
        }
    }

    /// A variant carrying the record's coordinates, strand and whole allele
    /// string, unnamed, for [`find_colocated_variants`] and [`to_vcf_record`].
    fn variant(&self) -> InputVariant {
        variant_for(
            &self.chr,
            self.start,
            self.end,
            self.strand,
            &self.allele_string,
        )
    }
}

/// The record's VCF forms, one per alternate allele, as `(chromosome-position-
/// reference-alternate, reference, alternate)` (`OutputFactory.pm:938-949`):
/// the line [`to_vcf_record`] builds from the record's alleles with the
/// reference's anchor base, whatever the input's format.
fn vcf_strings(head: &RecordHead, lookups: &Lookups) -> Vec<(String, String, String)> {
    let mut variant = head.variant();
    variant.id = Some(head.name.clone());
    let reference = lookups
        .reference
        .as_ref()
        .map(|fasta| Reference(Arc::clone(fasta)));
    let built = to_vcf_record(
        &variant,
        reference.as_ref().map(|r| r as &dyn ReferenceSequence),
    );
    let built = String::from_utf8_lossy(&built).into_owned();
    let columns: Vec<&str> = built.split('\t').collect();
    if columns.len() < 5 {
        return Vec::new();
    }
    let (position, reference, alts) = (columns[1], columns[3], columns[4]);
    alts.split(',')
        .map(|alt| {
            (
                format!("{}-{position}-{reference}-{alt}", head.chr),
                reference.to_string(),
                alt.to_string(),
            )
        })
        .collect()
}

/// A variant carrying the record's coordinates, strand and whole allele
/// string, for [`to_vcf_record`].
fn variant_for(
    chr: &str,
    start: u64,
    end: u64,
    strand: Strand,
    allele_string: &str,
) -> InputVariant {
    let mut parts = allele_string.split('/');
    let reference = parts.next().unwrap_or_default().as_bytes().to_vec();
    let alt = parts.next().unwrap_or("-").as_bytes().to_vec();
    let mut variant = InputVariant::new(chr.to_string(), start, end, reference, alt);
    variant.strand = strand;
    variant.allele_string = allele_string.to_string();
    variant
}

/// `_minimise_allele` with its flag set (502-531): each alternate allele's
/// minimised `reference/alternate` pair mapped to the allele as written.
fn minimised_alleles(allele_string: &str) -> HashMap<String, String> {
    let mut parts = allele_string.split('/');
    let reference = parts.next().unwrap_or_default();
    parts
        .map(|alt| (minimised_pair(reference, alt), alt.to_string()))
        .collect()
}

/// `minimise_alleles` (`Variation/Utils/VEP.pm:976-1036`) on one pair: the
/// shared leading bases are removed, then the shared trailing ones, and an
/// allele left empty is `-`.
fn minimised_pair(reference: &str, alt: &str) -> String {
    let mut reference = reference.to_string();
    let mut alt = alt.to_string();
    while !reference.is_empty() && !alt.is_empty() && reference.chars().next() == alt.chars().next()
    {
        reference.remove(0);
        alt.remove(0);
    }
    while !reference.is_empty()
        && !alt.is_empty()
        && reference.chars().next_back() == alt.chars().next_back()
    {
        reference.pop();
        alt.pop();
    }
    format!("{}/{}", dash_when_empty(reference), dash_when_empty(alt))
}

/// `-` for an empty text: an allele string's spelling of an emptied allele, and
/// a MANE Select triple's of a notation its transcript lacks.
fn dash_when_empty(text: String) -> String {
    if text.is_empty() {
        "-".to_string()
    } else {
        text
    }
}

/// `_convert_allele_to_vcf` (578-612): each alternate allele's VCF
/// `reference/alternate` pair, built without a reference so the anchor base
/// is `N` and removed, an allele left empty written `-`, mapped to the allele
/// as written.
fn alleles_in_vcf_form(allele_string: &str, strand: Strand) -> HashMap<String, String> {
    let mut parts = allele_string.split('/');
    let reference = parts.next().unwrap_or_default();
    let mut map = HashMap::new();
    for alt in parts {
        let variant = variant_for("1", 2, 2, strand, &format!("{reference}/{alt}"));
        let built = String::from_utf8_lossy(&to_vcf_record(&variant, None)).into_owned();
        let columns: Vec<&str> = built.split('\t').collect();
        if columns.len() < 5 {
            continue;
        }
        let strip_anchor = |column: &str| dash_when_empty(column.replacen('N', "", 1));
        let vcf_reference = strip_anchor(columns[3]);
        for vcf_alt in strip_anchor(columns[4]).split(',') {
            map.insert(format!("{vcf_reference}/{vcf_alt}"), alt.to_string());
        }
    }
    map
}

/// `_minimise_allele_vcf` (537-564): the allele as written whose minimised or
/// VCF pair the minimised VCF pair `reference/alt` matches.
fn allele_of_vcf_pair(
    reference: &str,
    alt: &str,
    minimised: &HashMap<String, String>,
    vcf_alleles: &HashMap<String, String>,
) -> Option<String> {
    let pair = minimised_pair(reference, alt);
    minimised
        .get(&pair)
        .or_else(|| vcf_alleles.get(&pair))
        .cloned()
}

/// The `source: id,id` strings of a cache `var_synonyms` column
/// (`source::id,id--source::id`), one per source, as
/// `Variation::get_all_synonyms` with its source flag writes them.
fn synonym_strings(column: &str) -> Vec<String> {
    column
        .split("--")
        .filter(|group| !group.is_empty())
        .map(|group| match group.split_once("::") {
            Some((source, ids)) => format!("{source}: {ids}"),
            None => format!("{group}: "),
        })
        .collect()
}

/// `ga4gh_vrs_from_spdi` (`Variation/Utils/Sequence.pm:1430-1472`): the GA4GH
/// VRS Allele of an SPDI on an `NC_` accession, its interval start the SPDI's
/// position as written and its end that position plus the deleted length;
/// nothing for an SPDI on any other sequence or not in SPDI form, the bare
/// `return` of 1433-1434 and 1438, which a caller pushing the result onto a
/// list appends nothing for.
fn ga4gh_vrs_from_spdi(spdi: &str) -> Option<Value> {
    if !spdi.starts_with("NC") {
        return None;
    }
    let parts: Vec<&str> = spdi.split(':').collect();
    let [accession, position, deleted, inserted] = parts.as_slice() else {
        return None;
    };
    if accession.is_empty()
        || position.is_empty()
        || !position.bytes().all(|b| b.is_ascii_digit())
        || !inserted
            .bytes()
            .all(|b| b.is_ascii_alphanumeric() || b == b'_')
    {
        return None;
    }
    let start = position.parse::<u64>().ok()?;
    let deleted_length = if !deleted.is_empty() && deleted.bytes().all(|b| b.is_ascii_digit()) {
        deleted.parse::<u64>().unwrap_or(0)
    } else {
        deleted.chars().count() as u64
    };
    Some(serde_json::json!({
        "location": {
            "type": "SequenceLocation",
            "interval": {"type": "SimpleInterval", "start": position, "end": start + deleted_length},
            "sequence_id": format!("refseq:{accession}"),
        },
        "state": {"type": "SequenceState", "sequence": inserted},
        "type": "Allele",
    }))
}

/// `find_in_ref` (`Utils.pm:337-368`): every scalar found under a wanted key
/// anywhere in `value`, added once to `out` under that key.
fn find_in_ref(
    value: &Value,
    want: &BTreeSet<String>,
    out: &mut IndexMap<String, Value>,
    this_key: Option<&str>,
) {
    match value {
        Value::Object(object) => {
            for (key, child) in object {
                if child.is_object() || child.is_array() {
                    find_in_ref(child, want, out, Some(key));
                } else if want.contains(key) {
                    push_unique(out, key, child.clone());
                }
            }
        }
        Value::Array(items) => {
            for item in items {
                if item.is_object() || item.is_array() {
                    find_in_ref(item, want, out, None);
                } else if let Some(key) = this_key.filter(|key| want.contains(*key)) {
                    push_unique(out, key, item.clone());
                }
            }
        }
        scalar => {
            if let Some(key) = this_key.filter(|key| want.contains(*key)) {
                push_unique(out, key, scalar.clone());
            }
        }
    }
}

/// `merge_arrays` (`Utils.pm:303-310`) of one value into the list under
/// `key`: added unless a value of the same text is there.
fn push_unique(out: &mut IndexMap<String, Value>, key: &str, value: Value) {
    let list = out
        .entry(key.to_string())
        .or_insert_with(|| Value::Array(Vec::new()));
    if let Value::Array(items) = list {
        let text = value_text(&value);
        if !items.iter().any(|item| value_text(item) == text) {
            items.push(value);
        }
    }
}

/// Every value appended to the list under `key`, as the `mane_select` triples
/// and the `ga4gh_vrs` objects are (one per triple or SPDI, no text to merge on).
fn push_all(out: &mut IndexMap<String, Value>, key: &str, values: Vec<Value>) {
    let list = out
        .entry(key.to_string())
        .or_insert_with(|| Value::Array(Vec::new()));
    if let Value::Array(items) = list {
        items.extend(values);
    }
}

/// A JSON scalar as Perl stringifies it: a string as itself, a number as
/// written, a boolean as `1` or `0`, null as the empty string; a list or
/// object as its JSON text.
fn value_text(value: &Value) -> String {
    match value {
        Value::String(s) => s.clone(),
        Value::Number(n) => n.to_string(),
        Value::Bool(true) => "1".to_string(),
        Value::Bool(false) => "0".to_string(),
        Value::Null => String::new(),
        other => other.to_string(),
    }
}

/// `JSON->new->pretty`'s layout: three spaces per level and a space on each
/// side of the colon.
struct PerlPretty<'a>(PrettyFormatter<'a>);

impl Formatter for PerlPretty<'_> {
    fn begin_array<W: ?Sized + io::Write>(&mut self, writer: &mut W) -> io::Result<()> {
        self.0.begin_array(writer)
    }

    fn end_array<W: ?Sized + io::Write>(&mut self, writer: &mut W) -> io::Result<()> {
        self.0.end_array(writer)
    }

    fn begin_array_value<W: ?Sized + io::Write>(
        &mut self,
        writer: &mut W,
        first: bool,
    ) -> io::Result<()> {
        self.0.begin_array_value(writer, first)
    }

    fn end_array_value<W: ?Sized + io::Write>(&mut self, writer: &mut W) -> io::Result<()> {
        self.0.end_array_value(writer)
    }

    fn begin_object<W: ?Sized + io::Write>(&mut self, writer: &mut W) -> io::Result<()> {
        self.0.begin_object(writer)
    }

    fn end_object<W: ?Sized + io::Write>(&mut self, writer: &mut W) -> io::Result<()> {
        self.0.end_object(writer)
    }

    fn begin_object_key<W: ?Sized + io::Write>(
        &mut self,
        writer: &mut W,
        first: bool,
    ) -> io::Result<()> {
        self.0.begin_object_key(writer, first)
    }

    fn begin_object_value<W: ?Sized + io::Write>(&mut self, writer: &mut W) -> io::Result<()> {
        writer.write_all(b" : ")
    }

    fn end_object_value<W: ?Sized + io::Write>(&mut self, writer: &mut W) -> io::Result<()> {
        self.0.end_object_value(writer)
    }
}

/// The JSON text as `variant_recoder:118-121` prints it: compact, or indented
/// with a final newline under `--pretty`.
fn write_json(out: &mut dyn Write, value: &Value, pretty: bool) -> io::Result<()> {
    if pretty {
        let mut serializer = serde_json::Serializer::with_formatter(
            &mut *out,
            PerlPretty(PrettyFormatter::with_indent(b"   ")),
        );
        serde::Serialize::serialize(value, &mut serializer)?;
        out.write_all(b"\n")
    } else {
        serde_json::to_writer(&mut *out, value)?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    fn args(argv: &[&str]) -> RecodeArgs {
        use clap::Parser;
        RecodeArgs::try_parse_from(std::iter::once("vep recode").chain(argv.iter().copied()))
            .unwrap()
    }

    /// The default fields and the appended ones settle as `new` settles them
    /// (111-150): hgvsc and hgvsp want the HGVS flag, hgvsg its own, spdi its
    /// own; `--mane_select` switches all three HGVS notations on; `--ga4gh_vrs`
    /// wants spdi beside it; `--var_synonyms` is not appended twice.
    #[test]
    fn fields_settle_the_flags_the_run_needs() {
        let fields = Fields::settle(&args(&["--id", "rs699"]));
        assert!(fields.id && fields.hgvs && fields.hgvsg && fields.spdi);
        assert!(
            !fields.vcf_string && !fields.var_synonyms && !fields.mane_select && !fields.ga4gh_vrs
        );
        assert_eq!(
            fields.want.iter().map(String::as_str).collect::<Vec<_>>(),
            ["hgvsc", "hgvsg", "hgvsp", "spdi"]
        );
        let fields = Fields::settle(&args(&["--id", "rs699", "--fields", "hgvsg"]));
        assert!(fields.hgvsg && !fields.hgvs && !fields.spdi && !fields.id);
        let fields = Fields::settle(&args(&[
            "--id",
            "rs699",
            "--fields",
            "hgvsc",
            "--mane_select",
        ]));
        assert!(fields.hgvs && fields.hgvsg && fields.mane_select && !fields.spdi);
        let fields = Fields::settle(&args(&["--id", "rs699", "--fields", "id", "--ga4gh_vrs"]));
        assert!(
            fields.ga4gh_vrs
                && fields.spdi
                && fields.want.contains("spdi")
                && fields.want.contains("ga4gh_vrs")
        );
        let fields = Fields::settle(&args(&[
            "--id",
            "rs699",
            "--fields",
            "var_synonyms,hgvsg",
            "--var_synonyms",
        ]));
        assert!(fields.var_synonyms && !fields.want.contains("var_synonyms"));
    }

    /// The notice has four forms: the clause on the transcript set stands for
    /// a cache read as the Ensembl set alone, and the synonyms are named where
    /// `var_synonyms` is written alone.
    #[test]
    fn notice_names_the_set_read_and_the_synonyms_written() {
        let set = "this cache is read as the Ensembl transcript set, and RefSeq notation is \
                   written from a merged cache; ";
        let lrg = "no variant is mapped to an LRG, which needs the Ensembl database";
        let lrg_and_synonyms = "no variant is mapped to an LRG and no Archive dbSNP synonym is \
                                written, which need the Ensembl database";
        assert_eq!(notice(false, false), format!("INFO: {set}{lrg}"));
        assert_eq!(
            notice(false, true),
            format!("INFO: {set}{lrg_and_synonyms}")
        );
        assert_eq!(notice(true, false), format!("INFO: {lrg}"));
        assert_eq!(notice(true, true), format!("INFO: {lrg_and_synonyms}"));
    }

    /// The line of a cache without variation data is drawn where the tool
    /// found none for fields that read known variants, naming `id`,
    /// `var_synonyms` or both as the fields are settled; a run that found the
    /// data, or whose fields read none, draws no line.
    #[test]
    fn the_line_of_a_cache_without_variation_data_names_the_fields_that_read_it() {
        let lookups = |under: Option<&str>| Lookups {
            reference: None,
            variations: None,
            no_variation_data_under: under.map(str::to_string),
            synonyms: Arc::new(ChromosomeSynonyms::default()),
            rules: ColocationRules::default(),
            rows: Vec::new(),
        };
        let line = |fields: &[&str]| {
            Some(cache_left_out::no_variation_data_recode_notice(
                fields, "/c",
            ))
        };
        let id = Fields::settle(&args(&["--id", "rs699", "--fields", "id"]));
        let both = Fields::settle(&args(&["--id", "rs699", "--var_synonyms"]));
        let synonyms = Fields::settle(&args(&["--id", "rs699", "--fields", "var_synonyms"]));
        let without = lookups(Some("/c"));
        assert_eq!(without.no_variation_data_notice(&id), line(&["id"]));
        assert_eq!(
            without.no_variation_data_notice(&both),
            line(&["id", "var_synonyms"])
        );
        assert_eq!(
            without.no_variation_data_notice(&synonyms),
            line(&["var_synonyms"])
        );
        assert_eq!(lookups(None).no_variation_data_notice(&both), None);
    }

    /// `--fields hgvs_offset` switches no HGVS notation on: Ensembl VEP 116.3
    /// sets the configuration key of each field starting with `hgvs`
    /// (`VariantRecoder.pm:120`), and the output reads `hgvsc` and `hgvsp`
    /// alone (`OutputFactory.pm:1759`, `:1767`), so `hgvs_offset` names a key
    /// nothing reads, while `hgvs` sets both (`Config.pm:383-386`). The key is
    /// wanted only beside `hgvsc`: the protein branch reads an offset only
    /// `hgvs_transcript` sets (`TranscriptVariationAllele.pm:1447`, `:1796`).
    #[test]
    fn hgvs_offset_field_switches_no_notation_on() {
        let fields = Fields::settle(&args(&["--id", "rs699", "--fields", "hgvs_offset"]));
        assert!(!fields.hgvs && !fields.hgvsg && !fields.spdi);
        assert!(fields.want.is_empty());
        let fields = Fields::settle(&args(&["--id", "rs699", "--fields", "hgvsc,hgvs_offset"]));
        assert!(fields.hgvs && !fields.hgvsg && fields.want.contains("hgvs_offset"));
        let fields = Fields::settle(&args(&["--id", "rs699", "--fields", "hgvsp,hgvs_offset"]));
        assert!(fields.hgvs);
        assert_eq!(
            fields.want.iter().map(String::as_str).collect::<Vec<_>>(),
            ["hgvsp"]
        );
        let fields = Fields::settle(&args(&["--id", "rs699", "--fields", "hgvs,hgvs_offset"]));
        assert!(fields.hgvs && fields.want.contains("hgvs_offset"));
        let fields = Fields::settle(&args(&[
            "--id",
            "rs699",
            "--fields",
            "hgvsp,hgvs_offset",
            "--mane_select",
        ]));
        assert!(fields.want.contains("hgvs_offset"));
    }

    /// The entries whose `hgvsc` the output writes follow the keys read, not
    /// the HGVS flag: every entry under the default fields and under `hgvsc`,
    /// with `--mane_select` as without; none under `hgvsp`, under `hgvs`,
    /// which is a key no entry carries, and under fields that need no HGVS
    /// flag; and the MANE Select entries alone where `--mane_select` writes
    /// its triples beside fields without `hgvsc`.
    #[test]
    fn fields_settle_the_entries_whose_hgvsc_is_written() {
        let rows = |more: &[&str]| {
            Fields::settle(&args(&[&["--id", "rs699"][..], more].concat())).hgvsc_rows()
        };
        assert_eq!(rows(&[]), NotationRows::Every);
        assert_eq!(rows(&["--fields", "hgvsc"]), NotationRows::Every);
        assert_eq!(rows(&["--fields", "id,hgvsc,hgvsp"]), NotationRows::Every);
        assert_eq!(
            rows(&["--fields", "hgvsc", "--mane_select"]),
            NotationRows::Every
        );
        for fields in [
            "hgvsp",
            "hgvs",
            "hgvs,hgvs_offset",
            "id,hgvsp,spdi",
            "hgvsg",
        ] {
            assert_eq!(rows(&["--fields", fields]), NotationRows::None, "{fields}");
            assert_eq!(
                rows(&["--fields", fields, "--mane_select"]),
                NotationRows::ManeSelect,
                "{fields}"
            );
        }
    }

    /// `find_in_ref` collects a wanted key's scalars from an entry, from a
    /// list under the key and from nested objects, once each, in order.
    #[test]
    fn find_in_ref_collects_unique_scalars_under_wanted_keys() {
        let want: BTreeSet<String> = ["hgvsc", "flags", "type"]
            .iter()
            .map(|s| s.to_string())
            .collect();
        let entries = json!([
            {"hgvsc": "ENST1:c.1A>G", "flags": ["cds_start_NF", "x"], "ga4gh_vrs": {"type": "Allele"}},
            {"hgvsc": "ENST1:c.1A>G", "flags": ["x"], "strand": 1},
            {"hgvsc": "ENST2:c.5A>G"}
        ]);
        let mut out = IndexMap::new();
        for entry in entries.as_array().unwrap() {
            find_in_ref(entry, &want, &mut out, None);
        }
        assert_eq!(out["hgvsc"], json!(["ENST1:c.1A>G", "ENST2:c.5A>G"]));
        assert_eq!(out["flags"], json!(["cds_start_NF", "x"]));
        assert_eq!(out["type"], json!(["Allele"]));
        assert!(!out.contains_key("strand"));
        push_unique(&mut out, "n", json!(1));
        push_unique(&mut out, "n", json!("1"));
        assert_eq!(out["n"], json!([1]));
    }

    /// `minimise_alleles` trims the shared leading bases, then the trailing
    /// ones, and writes an emptied allele `-`.
    #[test]
    fn minimised_pair_trims_both_ends() {
        assert_eq!(minimised_pair("CCCCCCC", "CCC"), "CCCC/-");
        assert_eq!(minimised_pair("GTTA", "G"), "TTA/-");
        assert_eq!(minimised_pair("C", "CC"), "-/C");
        assert_eq!(minimised_pair("C", "T"), "C/T");
        assert_eq!(minimised_pair("-", "A"), "-/A");
        assert_eq!(minimised_pair("ATA", "TAT"), "ATA/TAT");
    }

    /// `_convert_allele_to_vcf` maps each alternate's VCF pair, the anchor `N`
    /// removed and an emptied allele `-`, to the allele as written.
    #[test]
    fn alleles_in_vcf_form_match_perl() {
        let map = alleles_in_vcf_form("C/T", Strand::Forward);
        assert_eq!(map.get("C/T").map(String::as_str), Some("T"));
        let map = alleles_in_vcf_form("TTA/-", Strand::Forward);
        assert_eq!(map.get("TTA/-").map(String::as_str), Some("-"));
        let map = alleles_in_vcf_form("-/C", Strand::Forward);
        assert_eq!(map.get("-/C").map(String::as_str), Some("C"));
        let map = alleles_in_vcf_form("A/C/G", Strand::Reverse);
        assert_eq!(map.get("T/G").map(String::as_str), Some("C"));
        assert_eq!(map.get("T/C").map(String::as_str), Some("G"));
    }

    /// The VCF form of a record built without its own VCF line is matched to
    /// its allele directly, through the minimised pair, or through the pair
    /// reversed.
    #[test]
    fn allele_of_vcf_pair_resolves_minimised_and_reversed_pairs() {
        let minimised = minimised_alleles("GTTA/G");
        let vcf = alleles_in_vcf_form("GTTA/G", Strand::Forward);
        assert_eq!(
            allele_of_vcf_pair("GTTA", "G", &minimised, &vcf).as_deref(),
            Some("G")
        );
        let minimised = minimised_alleles("G/C");
        let vcf = alleles_in_vcf_form("G/C", Strand::Forward);
        assert_eq!(allele_of_vcf_pair("C", "G", &minimised, &vcf), None);
        assert_eq!(
            allele_of_vcf_pair("G", "C", &minimised, &vcf).as_deref(),
            Some("C")
        );
    }

    /// `ga4gh_vrs_from_spdi` writes the interval start as the SPDI's text and
    /// the end as a number, the deleted length from the sequence or the count,
    /// and nothing for an SPDI not on an `NC_` accession or not in SPDI form,
    /// as Ensembl VEP 116.3 returns the empty list for them
    /// (`Variation/Utils/Sequence.pm:1433-1434`, `:1438`).
    #[test]
    fn ga4gh_vrs_follows_the_spdi_and_is_nothing_off_an_nc_accession() {
        assert_eq!(
            ga4gh_vrs_from_spdi("NC_000021.9:25585732:C:T"),
            Some(json!({
                "location": {
                    "type": "SequenceLocation",
                    "interval": {"type": "SimpleInterval", "start": "25585732", "end": 25585733},
                    "sequence_id": "refseq:NC_000021.9"
                },
                "state": {"type": "SequenceState", "sequence": "T"},
                "type": "Allele"
            }))
        );
        let deletion = ga4gh_vrs_from_spdi("NC_000021.9:25585750:TTA:").unwrap();
        assert_eq!(deletion["location"]["interval"]["end"], json!(25585753));
        assert_eq!(deletion["state"]["sequence"], json!(""));
        let counted = ga4gh_vrs_from_spdi("NC_000021.9:100:3:AC").unwrap();
        assert_eq!(counted["location"]["interval"]["end"], json!(103));
        assert_eq!(ga4gh_vrs_from_spdi("21:25585732:C:T"), None);
        assert_eq!(ga4gh_vrs_from_spdi("LRG_18:12276:TG:"), None);
        assert_eq!(ga4gh_vrs_from_spdi("NC_000021.9:g.25585732C>T"), None);
    }

    /// Where the table lists one accession with and without its version, the
    /// SPDI's sequence is named by the versioned form, whichever of the
    /// sequence's names the cache holds it under: its own name, the accession
    /// without its version, or the versioned accession itself, which is kept.
    #[test]
    fn spdi_on_accession_writes_the_versioned_form() {
        let lookups = Lookups {
            reference: None,
            variations: None,
            no_variation_data_under: None,
            synonyms: Arc::new(ChromosomeSynonyms::parse(
                "21\tNC_000021\n21\tNC_000021.8\nNC_000021\tNC_000021.8\n",
            )),
            rules: ColocationRules::default(),
            rows: Vec::new(),
        };
        for name in ["21", "NC_000021", "NC_000021.8"] {
            assert_eq!(
                lookups.spdi_on_accession(&format!("{name}:99:1:2")),
                "NC_000021.8:99:1:2",
                "{name}"
            );
        }
    }

    /// Where the table lists two versions of one accession, the SPDI's sequence
    /// is named by the higher, whichever of the sequence's names the cache
    /// holds it under: its own name, the lower version, or the higher, which
    /// is kept.
    #[test]
    fn spdi_on_accession_writes_the_highest_version() {
        let lookups = Lookups {
            reference: None,
            variations: None,
            no_variation_data_under: None,
            synonyms: Arc::new(ChromosomeSynonyms::parse(
                "21\tNC_000021.7\n21\tNC_000021.8\nNC_000021.7\tNC_000021.8\n",
            )),
            rules: ColocationRules::default(),
            rows: Vec::new(),
        };
        for name in ["21", "NC_000021.7", "NC_000021.8"] {
            assert_eq!(
                lookups.spdi_on_accession(&format!("{name}:99:1:2")),
                "NC_000021.8:99:1:2",
                "{name}"
            );
        }
    }

    /// An allele whose SPDI stays on a chromosome name with no `NC_` synonym
    /// carries `ga4gh_vrs` as an empty list, and an allele with no SPDI carries
    /// no key: Ensembl VEP 116.3 skips an allele without `spdi` and pushes
    /// the empty list `ga4gh_vrs_from_spdi` returns for the SPDI onto the list
    /// the push itself makes (`VariantRecoder.pm:476-479`).
    #[test]
    fn ga4gh_vrs_is_an_empty_list_for_an_spdi_off_an_nc_accession() {
        let record = json!({
            "input": "rs1", "id": "rs1", "seq_region_name": "21", "start": 5, "end": 5, "strand": 1,
            "allele_string": "A/G/T",
            "intergenic_consequences": [
                {"variant_allele": "G", "spdi": "21:4:A:G"},
                {"variant_allele": "T", "hgvsg": "21:g.5A>T"}
            ]
        });
        let fields = Fields {
            want: ["spdi", "ga4gh_vrs"]
                .iter()
                .map(|key| key.to_string())
                .collect(),
            id: false,
            vcf_string: false,
            var_synonyms: false,
            mane_select: false,
            ga4gh_vrs: true,
            hgvs: false,
            hgvsg: false,
            spdi: true,
        };
        let out = regrouped_on(
            std::slice::from_ref(&record),
            &fields,
            ChromosomeSynonyms::parse("21 chr21"),
            Vec::new(),
        );
        assert_eq!(
            out[0]["G"],
            json!({"ga4gh_vrs": [], "input": "rs1", "spdi": ["21:4:A:G"]})
        );
        assert_eq!(out[0]["T"], json!({"input": "rs1"}));
        let out = regrouped(&[record], &fields);
        assert_eq!(out[0]["G"]["ga4gh_vrs"], json!([]));
    }

    /// A regrouping over hand-built records with no cache behind it: the
    /// keys named alone, no reference and no variation rows.
    fn regrouped(records: &[Value], fields: &Fields) -> Vec<Value> {
        regrouped_on(records, fields, ChromosomeSynonyms::default(), Vec::new())
    }

    /// [`regrouped`] with the chromosome synonyms the SPDI's accession comes
    /// from and the variation rows the cache would hold near the records.
    fn regrouped_on(
        records: &[Value],
        fields: &Fields,
        synonyms: ChromosomeSynonyms,
        rows: Vec<CachedVariation>,
    ) -> Vec<Value> {
        let records: Vec<Record> = records
            .iter()
            .enumerate()
            .map(|(i, value)| Record {
                line: Some(i as u64 + 1),
                value: value.as_object().unwrap().clone(),
            })
            .collect();
        let lookups = Lookups {
            reference: None,
            variations: None,
            no_variation_data_under: None,
            synonyms: Arc::new(synonyms),
            rules: ColocationRules {
                failed: 1,
                ..ColocationRules::default()
            },
            rows,
        };
        regroup(&records, &[], fields, &lookups, &[], 1).unwrap()
    }

    /// A cache variation row of chromosome 21.
    fn row(name: &str, start: u64, end: u64, allele_string: &str) -> CachedVariation {
        let mut variation = CachedVariation::default();
        variation.chr = Some("21".to_string());
        variation.variation_name = name.to_string();
        variation.start = start;
        variation.end = end;
        variation.allele_string = allele_string.to_string();
        variation.strand = 1;
        variation
    }

    fn id_fields() -> Fields {
        Fields {
            want: ["hgvsg".to_string()].into_iter().collect(),
            id: true,
            vcf_string: false,
            var_synonyms: false,
            mane_select: false,
            ga4gh_vrs: false,
            hgvs: false,
            hgvsg: true,
            spdi: false,
        }
    }

    /// The co-located variants' names go to the allele their alternate allele
    /// spells (379-421); a COSMIC or HGMD identifier, which has no alleles, is
    /// annexed to every allele (377-378, 413-420).
    #[test]
    fn cosmic_and_hgmd_ids_are_annexed_to_every_allele() {
        let record = json!({
            "input": "21\t100\t.\tC\tA,G\t.\t.\t.",
            "id": "21_100_C/A/G",
            "seq_region_name": "21", "start": 100, "end": 100, "strand": 1,
            "allele_string": "C/A/G",
            "transcript_consequences": [
                {"variant_allele": "A", "hgvsg": "21:g.100C>A"},
                {"variant_allele": "G", "hgvsg": "21:g.100C>G"}
            ]
        });
        let rows = vec![
            row("rs1", 100, 100, "C/A"),
            row("COSV1", 100, 100, "COSMIC_MUTATION"),
            row("CM001", 100, 100, "HGMD_MUTATION"),
            row("rs2", 100, 100, "C/A/G"),
        ];
        let out = regrouped_on(&[record], &id_fields(), ChromosomeSynonyms::default(), rows);
        assert_eq!(out[0]["A"]["id"], json!(["rs1", "rs2", "COSV1", "CM001"]));
        assert_eq!(out[0]["G"]["id"], json!(["rs2", "COSV1", "CM001"]));
        assert_eq!(out[0]["A"]["hgvsg"], json!(["21:g.100C>A"]));
    }

    /// A co-located variant written on the other strand reaches its allele by
    /// reverse complement (397-401), one written unminimised through the
    /// minimised allele (402-404), and one spelling no allele of the record
    /// reaches none.
    #[test]
    fn colocated_ids_match_alleles_by_reverse_complement_and_minimised_form() {
        let record = json!({
            "input": "rs9",
            "id": "rs9",
            "seq_region_name": "21", "start": 100, "end": 100, "strand": 1,
            "allele_string": "C/A",
            "transcript_consequences": [
                {"variant_allele": "A", "hgvsg": "21:g.100C>A"}
            ]
        });
        let mut reversed = row("rs_rev", 100, 100, "G/T");
        reversed.strand = -1;
        let rows = vec![
            reversed,
            row("rs_long", 98, 100, "CCC/CCA"),
            row("rs_other", 100, 100, "C/G"),
        ];
        let out = regrouped_on(&[record], &id_fields(), ChromosomeSynonyms::default(), rows);
        assert_eq!(
            out[0]["A"]["id"],
            json!(["rs_rev", "rs_long"]),
            "{}",
            out[0]
        );
        assert_eq!(out[0].as_object().unwrap().len(), 1);
    }

    /// The co-located variants are the cache rows a database-backed run finds
    /// (`AnnotationType/Variation.pm:91-109`): a row is kept when its span
    /// overlaps the record's and an allele matches, however the two are
    /// spelled, so a deletion minimised to its last base finds the row dbSNP
    /// spells from the run's first; an insertion overlaps only a row spanning
    /// both flanking bases, so its own row, spelled from the repeat's first
    /// base, is not co-located and gives no `id`, while its synonyms are still
    /// the row's, taken by the input's name (`OutputFactory.pm:955-957`).
    #[test]
    fn colocated_rows_are_found_by_span_overlap_and_allele_match() {
        let deletion = json!({
            "input": "NC_000021.9:34886861:G:", "id": "NC_000021.9:34886861:G:",
            "seq_region_name": "21", "start": 34886862, "end": 34886862, "strand": 1,
            "allele_string": "G/-",
            "transcript_consequences": [{"variant_allele": "-", "hgvsg": "21:g.34886862del"}]
        });
        let insertion = json!({
            "input": "rs1601245742", "id": "rs1601245742",
            "seq_region_name": "21", "start": 33543120, "end": 33543119, "strand": 1,
            "allele_string": "-/T",
            "transcript_consequences": [{"variant_allele": "T", "hgvsg": "21:g.33543119dup"}]
        });
        let mut deletion_row = row("rs2057997486", 34886860, 34886862, "GGG/GG");
        deletion_row.var_synonyms = Some("ClinVar::VCV000869209".to_string());
        let mut insertion_row = row("rs1601245742", 33543115, 33543119, "TTTTT/TTTTTT");
        insertion_row.var_synonyms = Some("ClinVar::VCV000817560".to_string());
        let rows = vec![
            insertion_row,
            deletion_row,
            row("rs_before", 34886858, 34886859, "AC/A"),
        ];
        let mut fields = id_fields();
        fields.var_synonyms = true;
        let out = regrouped_on(
            &[deletion, insertion],
            &fields,
            ChromosomeSynonyms::default(),
            rows,
        );
        assert_eq!(out[0]["-"]["id"], json!(["rs2057997486"]));
        assert_eq!(
            out[0]["-"]["var_synonyms"],
            json!(["ClinVar: VCV000869209"])
        );
        assert!(out[1]["T"].get("id").is_none(), "{}", out[1]);
        assert_eq!(
            out[1]["T"]["var_synonyms"],
            json!(["ClinVar: VCV000817560"])
        );
    }

    /// The synonyms of a row whose `var_synonyms` column the cache cut inside
    /// an id, at 1,024 bytes, are the ids before the cut one: the stored value
    /// is a `UniProt` id, 76 `ClinVar` ids and the first six characters of
    /// another, every id invented.
    #[test]
    fn synonyms_the_cache_cut_inside_an_id_end_before_it() {
        let deletion = json!({
            "input": "rs9920000001", "id": "rs9920000001",
            "seq_region_name": "21", "start": 34886862, "end": 34886862, "strand": 1,
            "allele_string": "G/-",
            "transcript_consequences": [{"variant_allele": "-", "hgvsg": "21:g.34886862del"}]
        });
        let ids: Vec<String> = (1..=76)
            .map(|n| format!("RCV{}", 999_000_000 + n))
            .collect();
        let whole = ids.join(",");
        let stored = format!("UniProt::VAR_999001--ClinVar::{whole},RCV999");
        assert_eq!(stored.len(), 1024);
        let mut cut_row = row("rs9920000001", 34886860, 34886862, "GGG/GG");
        cut_row.var_synonyms = Some(stored);
        let mut fields = id_fields();
        fields.var_synonyms = true;
        let out = regrouped_on(
            &[deletion],
            &fields,
            ChromosomeSynonyms::default(),
            vec![cut_row],
        );
        assert_eq!(
            out[0]["-"]["var_synonyms"],
            json!(["UniProt: VAR_999001", format!("ClinVar: {whole}")])
        );
    }

    /// The rows the co-located variants are gathered from are the ones a
    /// database-backed run reads (`Database/Variation.pm:167-179`): the rows
    /// starting in the region of the record's span widened by one base. Of two
    /// rows that both spell the record's deletion, the one starting in the
    /// region before the record's is not read, however far its span reaches
    /// into it, while the one starting in the record's region is.
    #[test]
    fn rows_near_read_the_regions_of_the_record_by_start() {
        let cache = tempfile::tempdir().unwrap();
        let write_region = |start: u64, end: u64, line: &str| {
            let dir = cache.path().join("variations").join("21");
            fs::create_dir_all(&dir).unwrap();
            let mut encoder = flate2::write::GzEncoder::new(
                fs::File::create(dir.join(format!("{start}-{end}_var.gz"))).unwrap(),
                flate2::Compression::default(),
            );
            writeln!(encoder, "{line}").unwrap();
            encoder.finish().unwrap();
        };
        write_region(1, 1_000_000, "rs_before   1000000 1000002 GAC/GA");
        write_region(1_000_001, 2_000_000, "rs_here   1000001 1000002 AC/A");
        let mut info = vep_core::cache_info::CacheInfo::default();
        info.variation_cols = "variation_name,failed,somatic,start,end,allele_string,strand"
            .split(',')
            .map(str::to_string)
            .collect();
        let manifest = vep_core::cache_manifest::CacheManifest::absent(cache.path());
        let source = VariationSource::open(cache.path(), &info, false, Arc::new(manifest))
            .unwrap()
            .unwrap();
        let record = json!({
            "input": "21:g.1000002del", "id": "21:g.1000002del",
            "seq_region_name": "21", "start": 1000002, "end": 1000002, "strand": 1,
            "allele_string": "C/-",
            "transcript_consequences": [{"variant_allele": "-", "hgvsg": "21:g.1000002del"}]
        });
        // Read together, both rows co-locate with the record.
        let rows = vec![
            row("rs_before", 1_000_000, 1_000_002, "GAC/GA"),
            row("rs_here", 1_000_001, 1_000_002, "AC/A"),
        ];
        let out = regrouped_on(
            std::slice::from_ref(&record),
            &id_fields(),
            ChromosomeSynonyms::default(),
            rows,
        );
        assert_eq!(out[0]["-"]["id"], json!(["rs_before", "rs_here"]));
        let lookups = Lookups {
            reference: None,
            variations: Some(source),
            no_variation_data_under: None,
            synonyms: Arc::new(ChromosomeSynonyms::default()),
            rules: ColocationRules {
                failed: 1,
                ..ColocationRules::default()
            },
            rows: Vec::new(),
        };
        let records = [Record {
            line: Some(1),
            value: record.as_object().unwrap().clone(),
        }];
        let out = regroup(&records, &[], &id_fields(), &lookups, &[], 1).unwrap();
        assert_eq!(out[0]["-"]["id"], json!(["rs_here"]), "{}", out[0]);
    }

    /// An input whose record has no consequence entry is `null`, and the same
    /// input given twice is one object (`merge_arrays` over the order, 465).
    #[test]
    fn null_without_consequences_and_one_object_per_input_text() {
        let empty = json!({"input": "rs0", "id": "rs0", "allele_string": "DBSNP_NOVARIATION"});
        let record = json!({
            "input": "rs1", "id": "rs1", "seq_region_name": "21", "start": 5, "end": 5, "strand": 1,
            "allele_string": "A/G",
            "intergenic_consequences": [{"variant_allele": "G", "hgvsg": "21:g.5A>G"}]
        });
        let out = regrouped(&[empty, record.clone(), record], &id_fields());
        assert_eq!(out.len(), 2);
        assert_eq!(out[0], Value::Null);
        assert_eq!(out[1]["G"]["hgvsg"], json!(["21:g.5A>G"]));
    }

    /// The same input given twice carries its SPDI once, on the chromosome's
    /// accession, as `merge_arrays` keeps one text (465), and one GA4GH VRS
    /// allele per record of it (474-482).
    #[test]
    fn repeated_input_keeps_one_spdi_on_the_accession() {
        let record = json!({
            "input": "rs1", "id": "rs1", "seq_region_name": "21", "start": 5, "end": 5, "strand": 1,
            "allele_string": "A/G",
            "intergenic_consequences": [{"variant_allele": "G", "spdi": "21:4:A:G"}]
        });
        let fields = Fields {
            want: ["spdi", "ga4gh_vrs"]
                .iter()
                .map(|key| key.to_string())
                .collect(),
            id: false,
            vcf_string: false,
            var_synonyms: false,
            mane_select: false,
            ga4gh_vrs: true,
            hgvs: false,
            hgvsg: false,
            spdi: true,
        };
        let out = regrouped_on(
            &[record.clone(), record],
            &fields,
            ChromosomeSynonyms::parse("21 NC_000021.9"),
            Vec::new(),
        );
        assert_eq!(out.len(), 1);
        assert_eq!(out[0]["G"]["spdi"], json!(["NC_000021.9:4:A:G"]));
        assert_eq!(out[0]["G"]["ga4gh_vrs"].as_array().unwrap().len(), 2);
    }

    /// A run with no object is one object of its warnings, an empty list when
    /// there was none (`variant_recoder:116`).
    #[test]
    fn a_run_without_objects_is_one_object_of_its_warnings() {
        assert_eq!(regrouped(&[], &id_fields()), [json!({"warnings": []})]);
        let lookups = Lookups {
            reference: None,
            variations: None,
            no_variation_data_under: None,
            synonyms: Arc::new(ChromosomeSynonyms::default()),
            rules: ColocationRules::default(),
            rows: Vec::new(),
        };
        let collected = vec![(
            None,
            "WARNING: No variant found with ID 'rs0'\n".to_string(),
        )];
        assert_eq!(
            regroup(&[], &[], &id_fields(), &lookups, &collected, 1).unwrap(),
            [json!({"warnings": ["No variant found with ID 'rs0'"]})]
        );
    }

    /// The VCF form is built from the record's alleles whatever the input's
    /// format, as `to_VCF_record` builds it for every variation feature an
    /// input parser makes: a VCF record's columns are not copied, so an indel
    /// the record spells from a repeat's first base is written from its
    /// minimised position with the base before it, and a balanced
    /// substitution's form is its own line.
    #[test]
    fn vcf_strings_are_rebuilt_for_every_input_format() {
        let lookups = Lookups {
            reference: None,
            variations: None,
            no_variation_data_under: None,
            synonyms: Arc::new(ChromosomeSynonyms::default()),
            rules: ColocationRules::default(),
            rows: Vec::new(),
        };
        let head = |start: u64, end: u64, allele_string: &str| RecordHead {
            chr: "21".to_string(),
            start,
            end,
            strand: Strand::Forward,
            allele_string: allele_string.to_string(),
            name: "21\t44293150\trs759281545\tTGTGTGTG\tTGTGTG\t.\t.\t.".to_string(),
        };
        let forms = vcf_strings(&head(44293156, 44293157, "TG/-"), &lookups);
        assert_eq!(
            forms,
            [(
                "21-44293155-NTG-N".to_string(),
                "NTG".to_string(),
                "N".to_string()
            )]
        );
        let forms = vcf_strings(&head(100, 100, "C/A/G"), &lookups);
        assert_eq!(
            forms
                .iter()
                .map(|(text, _, _)| text.as_str())
                .collect::<Vec<_>>(),
            ["21-100-C-A", "21-100-C-G"]
        );
    }

    /// A scalar's text is Perl's: a string as itself, a number as written,
    /// `1` and `0` for the booleans, the empty string for null.
    #[test]
    fn value_text_stringifies_as_perl() {
        assert_eq!(value_text(&json!("rs1")), "rs1");
        assert_eq!(value_text(&json!(25702257)), "25702257");
        assert_eq!(value_text(&json!(true)), "1");
        assert_eq!(value_text(&json!(false)), "0");
        assert_eq!(value_text(&Value::Null), "");
    }

    /// The synonyms column becomes one `source: id,id` string per source.
    #[test]
    fn synonym_strings_name_each_source_once() {
        assert_eq!(
            synonym_strings("ClinVar::RCV000003481,VCV000003317--OMIM::607358.0011"),
            ["ClinVar: RCV000003481,VCV000003317", "OMIM: 607358.0011"]
        );
    }

    /// A warning text yields its message after the prefix, a `MSG:` line
    /// extends it, and every other line is dropped (`Runner.pm:321-357`).
    #[test]
    fn warning_messages_follow_the_runner() {
        assert_eq!(
            warning_messages("WARNING: line 2 skipped (x): Unable to parse\n"),
            ["line 2 skipped (x): Unable to parse"]
        );
        assert_eq!(
            warning_messages("WARNING: first\nMSG: more\nSTACK Foo::bar\nERROR: second\n"),
            ["first: more", "second"]
        );
        assert!(warning_messages("STACK only\n").is_empty());
    }

    fn record(line: Option<u64>, input: &str) -> Record {
        let mut value = Map::new();
        value.insert("input".to_string(), Value::String(input.to_string()));
        Record { line, value }
    }

    /// With the recoder's buffer of one, a skipped line's warning goes to the
    /// next object written and a trailing one is lost; with a wider buffer every
    /// warning of the one fill goes to its first object.
    #[test]
    fn warnings_attach_to_the_first_object_of_their_fill() {
        let records = [
            record(Some(1), "a"),
            record(Some(3), "c"),
            record(Some(4), "d"),
        ];
        let collected = vec![
            (Some(2), "WARNING: line 2 skipped (b): bad\n".to_string()),
            (Some(5), "WARNING: line 5 skipped (e): bad\n".to_string()),
            (None, "WARNING: unnumbered\n".to_string()),
        ];
        let lines: Vec<(u64, String)> = ["a", "b", "c", "d", "e"]
            .iter()
            .enumerate()
            .map(|(i, l)| (i as u64 + 1, l.to_string()))
            .collect();
        let attached = attach_warnings(&records, &lines, &collected, 1);
        assert_eq!(attached[&1], ["line 2 skipped (b): bad"]);
        assert_eq!(attached[&0], ["unnumbered"]);
        assert!(!attached.contains_key(&2));
        // An unnumbered warning that quotes an input line takes that line's number.
        let quoted = vec![(None, "WARNING: No variant found with ID 'b'\n".to_string())];
        assert_eq!(
            attach_warnings(&records, &lines, &quoted, 1)[&1],
            ["No variant found with ID 'b'"]
        );
        let attached = attach_warnings(&records, &lines, &collected, 5000);
        assert_eq!(
            attached[&0],
            [
                "line 2 skipped (b): bad",
                "line 5 skipped (e): bad",
                "unnumbered"
            ]
        );
        assert_eq!(attached.len(), 1);
    }

    /// `--pretty` lays the JSON out as `JSON->pretty` does: three-space
    /// indents, ` : ` after each key, a final newline; the compact form has
    /// neither spaces nor a newline. Keys are written in sorted order.
    #[test]
    fn pretty_output_follows_the_perl_json_layout() {
        let value = json!([{"T": {"input": "rs1", "id": ["rs1"], "n": []}}]);
        let mut pretty = Vec::new();
        write_json(&mut pretty, &value, true).unwrap();
        assert_eq!(
            String::from_utf8(pretty).unwrap(),
            "[\n   {\n      \"T\" : {\n         \"id\" : [\n            \"rs1\"\n         ],\n         \"input\" : \"rs1\",\n         \"n\" : []\n      }\n   }\n]\n"
        );
        let mut compact = Vec::new();
        write_json(&mut compact, &value, false).unwrap();
        assert_eq!(
            String::from_utf8(compact).unwrap(),
            "[{\"T\":{\"id\":[\"rs1\"],\"input\":\"rs1\",\"n\":[]}}]"
        );
    }

    /// The input of a run, given as a string.
    fn staged(data: &str) -> StagedInput {
        StagedInput {
            input_file: None,
            input_data: Some(data.to_string()),
            _stdin_copy: None,
        }
    }

    /// Each value reaches the annotation run joined to its option, so that
    /// run's parser takes a value that begins with a hyphen as the value it
    /// is: the input, the species, a transcript filter and the output path.
    #[test]
    fn values_reach_the_annotation_run_joined_to_their_options() {
        let args = args(&[
            "--id=--help",
            "--species=-s",
            "--transcript_filter",
            "-x is 1",
            "--cache_version",
            "116",
            "--buffer_size",
            "2",
        ]);
        let fields = Fields::settle(&args);
        let input = staged("--help");
        let output = Path::new("-out.json");
        let argv: Vec<String> = inner_argv(&args, &fields, &input, output, false, None)
            .iter()
            .map(|token| token.to_string_lossy().into_owned())
            .collect();
        for joined in [
            "--distance=0",
            "--output_file=-out.json",
            "--input_data=--help",
            "--species=-s",
            "--cache_version=116",
            "--buffer_size=2",
            "--transcript_filter=-x is 1",
        ] {
            assert!(
                argv.iter().any(|token| token == joined),
                "{joined}: {argv:?}"
            );
        }
        let config = inner_config(&args, &fields, &input, output, false, None).unwrap();
        assert_eq!(config.input_data.as_deref(), Some("--help"));
        assert_eq!(config.species, "-s");
        assert_eq!(config.transcript_filter, ["-x is 1"]);
        assert_eq!(config.output_file, "-out.json");
        assert_eq!(config.buffer_size, 2);
    }

    /// Arguments built in code can carry a value outside its option's list,
    /// which the tool's own parser refuses on a command line: the annotation
    /// run's refusal of it is one line that names the value.
    #[test]
    fn a_value_set_in_code_that_the_annotation_run_refuses_is_refused_in_one_line() {
        let mut args = args(&["--id", "rs699"]);
        args.format = Some("foo".to_string());
        let fields = Fields::settle(&args);
        let refusal = inner_config(
            &args,
            &fields,
            &staged("rs699"),
            Path::new("out.json"),
            false,
            None,
        )
        .unwrap_err();
        assert_eq!(
            refusal.to_string(),
            "ERROR: [recode] the annotation run refused its arguments: invalid value 'foo' for \
             '--format <FORMAT>'"
        );
    }
}
