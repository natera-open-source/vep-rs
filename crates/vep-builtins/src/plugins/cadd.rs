// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! CADD plugin: Combined Annotation Dependent Depletion scores.
//!
//! Tabix plugin with allele matching. Annotates sequence variants with the
//! PHRED-scaled and raw scores of the CADD SNV and indel files, and structural
//! variants with the scores of a CADD-SV file, whose rows name a class (`INS`,
//! `DEL`, `DUP`; CADD.pm:75-79) over a span. Results are stored at the variant level
//! (`plugin_data` on `InputVariant`) and print on every row of the record, as
//! Ensembl VEP's plugin declares `feature_types` of `Feature` and `Intergenic`.
//!
//! `CADD.pm` is the plugin's Ensembl VEP counterpart in Ensembl/VEP_plugins
//! release/116. Its parameters are either every one `key=value` (`snv=`,
//! `indels=`, `sv=`, `force_annotate=`, and any other key a file as well) or
//! every one a file by position, where a trailing `1` is `force_annotate`
//! (CADD.pm:118-132). A structural variant is not annotated when a file was
//! given under `snv=` or `indels=` unless `force_annotate` is set
//! (CADD.pm:218), and no variant is annotated from a region the files answer
//! with more than 100,000 rows unless it is (CADD.pm:230).

use std::borrow::Cow;

use indexmap::IndexMap;
use tracing::{debug, warn};

use crate::annotation_store::{check_table_or_store, open_best_store_with, AnnotationStore};
use crate::custom::grammar::truthy;
use crate::params::{Grammar, PluginParams};
use crate::pool_pass::map_in_order;
use crate::tabix::allele_match::{self, KeyPosition};
use crate::tabix::{
    find_index_path, header_lines, indexed_table_config, split_field_count, table_contig,
    BatchQueryResult, TabixAnnotator, TabixConfig, TabixRecord,
};
use crate::traits::{BuiltinPlugin, PluginError, PrefetchData};
use vep_core::allele_trim::{strip_shared_first_base, trim_directions, trim_sequences, CaseRule};
use vep_core::assembly::Assembly;
use vep_core::consequence::{RowKind, TranscriptConsequence};
use vep_core::variant::InputVariant;

const NAME: &str = "CADD";

/// The one flag token of the positional grammar (CADD.pm:119-120).
const POSITIONAL_FLAGS: [&str; 1] = ["1"];

/// Rows answered for one variant's region above which nothing is annotated
/// without `force_annotate` (CADD.pm:230).
const RECORD_CAP: usize = 100_000;

/// The header's text on the data's terms, appended to each field's description
/// (CADD.pm:100, 174).
const NON_COMMERCIAL_USE_CLAUSE: &str =
    "CADD is only available here for non-commercial use. See CADD website for more information.";

const PHRED_FIELD: &str = "CADD_PHRED";
const RAW_FIELD: &str = "CADD_RAW";
const PHRED_DESCRIPTION: &str = "PHRED-like scaled CADD score.";
const RAW_DESCRIPTION: &str = "Raw CADD score.";

/// The score columns a file may carry and the field each one writes
/// (CADD.pm:81-98): the SNV and indel files name `PHRED` and `RawScore`, a
/// CADD-SV file `CADD-SV_PHRED-score` and `CADD-SV_Raw-score`.
const PHRED_COLUMNS: [&str; 2] = ["PHRED", "CADD-SV_PHRED-score"];
const RAW_COLUMNS: [&str; 2] = ["RawScore", "CADD-SV_Raw-score"];

/// Ensembl VEP's warning when `force_annotate` is set (CADD.pm:137-139).
const FORCE_ANNOTATE_WARNING: &str =
    "Using snv and/or indels CADD annotation file with structural \
     variant can increase run time exponentially. Consider creating separate input files for \
     SNV/indels and SV and use appropriate CADD annotation file.";

/// A score file with the columns `parse_data` reads by name (CADD.pm:265-268),
/// resolved from its header once.
struct CaddSource {
    store: Box<dyn AnnotationStore>,
    /// The file as given, named in the error for a row the match dies on.
    path: String,
    /// `Pos`, or `Start` where the header has no `Pos`: the row's own position, matched as
    /// written whatever coordinate system the index declares.
    position_col: usize,
    /// `End`, when the file has the column: a Perl-true `End` cell spares the
    /// position the arithmetic that numifies it (`$data{"End"} || ($s +
    /// length($ref)) - 1`, CADD.pm:278), so the shared-base move then steps
    /// it as a string.
    end_col: Option<usize>,
    /// `Ref`, when the file has the column; without it the reference is `-`.
    ref_col: Option<usize>,
    /// `Alt`, when the file has the column: the alternate allele of a SNV or
    /// indel row.
    alt_col: Option<usize>,
    /// `Type`, when the file has the column: the class word of a CADD-SV row,
    /// and the alternate allele of any row whose `Alt` cell is Perl-false
    /// (`$data{"Alt"} || $data{"Type"}`, CADD.pm:268). A file has one of the
    /// two columns at least.
    type_col: Option<usize>,
    /// The PHRED-scaled score column, when the file has one.
    phred_col: Option<usize>,
    /// The raw score column, when the file has one.
    raw_col: Option<usize>,
}

/// The query configuration the header of a score file names, for a file that
/// has no tabix index: the `#Chrom` column, which opens the header line, the
/// position column and the `End` column where the header has one, with a
/// zero-based start where the position column is `Start`.
///
/// The plugin's own comparison gives each column its coordinate. A sequence
/// variant is matched at its start against the row's `Pos` cell (CADD.pm:209,
/// 243, 249, 266), so a `Pos` cell is a one-based position. A structural
/// variant is matched at its start less one against the row's `Start` cell
/// (CADD.pm:220, 243, 249, 266), so a `Start` cell is a zero-based start,
/// which `tabix` reads as the span from the base after it to the row's end.
/// The lookup runs from two bases before the matched position to the
/// variant's end (CADD.pm:227), and the plugin filters no returned row by
/// position (CADD.pm:110-111).
fn header_named_config(
    file_path: std::path::PathBuf,
    position_col: usize,
    position_is_start: bool,
    end_col: Option<usize>,
) -> TabixConfig {
    TabixConfig {
        file_path,
        chr_col: 0,
        start_col: position_col,
        end_col,
        ref_col: None,
        alt_col: None,
        zero_based: position_is_start,
    }
}

/// One side of an allele match: a reference and alternate allele at a position,
/// `-` for an empty allele, as `get_matched_variant_alleles` takes them.
#[derive(Clone, Copy)]
struct Alleles<'a> {
    reference: &'a [u8],
    alternate: &'a [u8],
    position: KeyPosition<'a>,
}

/// A side after `trim_sequences` (ensembl-variation Utils/Sequence.pm:965-1038)
/// in one trim order, compared as the key `"$ref_$alt_$pos"` is (:1207, 1232):
/// the alleles as text, the position as [`KeyPosition::equals`] reads it.
type Trimmed<'a> = (&'a [u8], &'a [u8], KeyPosition<'a>);

/// A side in every trim order it is compared in: the start-first order, and the
/// end-first order as well when either allele is longer than one base
/// (`_get_trim_directions`, Sequence.pm:1250-1258).
struct TrimmedForms<'a> {
    start_first: Trimmed<'a>,
    end_first: Option<Trimmed<'a>>,
}

impl<'a> Alleles<'a> {
    /// The common bases trimmed from both alleles, from the start and then the
    /// end, or with `end_first` the other way round; the position moves past
    /// each base trimmed from the start, and an allele trimmed to nothing is
    /// `-` (`trim_sequences` with `$empty_to_dash`, Sequence.pm:976-1035). The
    /// bases the alleles share are found under `case`.
    fn trimmed(self, end_first: bool, case: CaseRule) -> Trimmed<'a> {
        let trimmed = trim_sequences(self.reference, self.alternate, end_first, case);
        (
            trimmed.reference,
            trimmed.alt,
            self.position.moved_by(trimmed.left as u64),
        )
    }

    fn trimmed_forms(self, case: CaseRule) -> TrimmedForms<'a> {
        let mut orders = trim_directions(self.reference, self.alternate).iter();
        TrimmedForms {
            start_first: self.trimmed(*orders.next().expect("the start-first order"), case),
            end_first: orders
                .next()
                .map(|end_first| self.trimmed(*end_first, case)),
        }
    }
}

impl<'a> TrimmedForms<'a> {
    fn forms(&self) -> impl Iterator<Item = Trimmed<'a>> {
        std::iter::once(self.start_first).chain(self.end_first)
    }

    /// Whether the two sides name one allele change: a form of one equals a
    /// form of the other (`get_matched_variant_alleles`, Sequence.pm:1190-1241),
    /// the alleles compared under `case`.
    fn matches(&self, other: &TrimmedForms<'_>, case: CaseRule) -> bool {
        self.forms().any(|(reference, alternate, position)| {
            other
                .forms()
                .any(|(other_reference, other_alternate, other_position)| {
                    case.same_allele(other_reference, reference)
                        && case.same_allele(other_alternate, alternate)
                        && position.equals(other_position)
                })
        })
    }
}

/// The variant's side of the match, owning the forward-strand alleles of a
/// reverse-strand record.
struct VariantSide<'v> {
    reference: Cow<'v, [u8]>,
    alternate: Cow<'v, [u8]>,
    /// The position matched: the variant's start, or one base before a
    /// structural variant's (CADD.pm:209, 220).
    position: u64,
    /// How a table row's alleles are compared with this side's. A sequence
    /// variant's row is read in either letter case, where Ensembl VEP's key
    /// is compared as written (Sequence.pm:1207, 1232; the documented class
    /// `tabix_allele_case_insensitive_match`); a structural variant's class
    /// word is compared as written.
    case: CaseRule,
}

impl VariantSide<'_> {
    fn alleles(&self) -> Alleles<'_> {
        Alleles {
            reference: &self.reference,
            alternate: &self.alternate,
            position: KeyPosition::Number(self.position),
        }
    }

    /// The region the variant's rows are read from: two bases before the
    /// matched position to the variant's end (CADD.pm:227).
    fn region(&self, variant: &InputVariant) -> (u64, u64) {
        (self.position.saturating_sub(2), variant.end)
    }
}

/// A row's side as [`row_side`] builds it: the alleles borrowed from the
/// record, the position cell as the key reads it (stepped as text where
/// `parse_data` stepped it as a string) and the bases moved off it.
struct RowSide<'r> {
    reference: &'r [u8],
    alternate: &'r [u8],
    position: Cow<'r, str>,
    moved: u64,
}

impl RowSide<'_> {
    fn alleles(&self) -> Alleles<'_> {
        Alleles {
            reference: self.reference,
            alternate: self.alternate,
            position: KeyPosition::Cell {
                text: &self.position,
                moved: self.moved,
            },
        }
    }
}

/// CADD pathogenicity score plugin.
///
/// Annotates variants with CADD PHRED and raw scores by querying one or
/// more tabix-indexed score files. Allele matching is performed to ensure
/// the annotation corresponds to the correct ref/alt pair.
pub struct CaddPlugin {
    /// The score files in the order given, read after `init()`.
    sources: Vec<CaddSource>,
    /// `force_annotate`: annotate a structural variant from the SNV and indel
    /// files, and a variant whose region answers more than [`RECORD_CAP`] rows.
    force_annotate: bool,
    /// A file was given under `snv=` or `indels=` (CADD.pm:116, 128): without
    /// `force_annotate`, no structural variant is annotated (CADD.pm:218).
    sequence_file_named: bool,
    /// The run's assembly (the cache's own, else the value of `--assembly`),
    /// injected before `init()`; a data file whose header does not name it is
    /// rejected (CADD.pm:151-157).
    assembly: Option<Assembly>,
}

impl Default for CaddPlugin {
    fn default() -> Self {
        Self::new()
    }
}

impl CaddPlugin {
    /// Create an uninitialized CADD plugin.
    pub fn new() -> Self {
        Self {
            sources: Vec::new(),
            force_annotate: false,
            sequence_file_named: false,
            assembly: None,
        }
    }

    /// Open a score file: its header lines are checked for the assembly and the
    /// column names, its index read for the columns and coordinate system the
    /// region queries use, as `tabix` reads them. A file with no index is read
    /// from the `.vpd` store beside it alone, and only from a store built on
    /// the columns the header names (`header_named_config`).
    fn open_source(path: &str, requested: Option<Assembly>) -> Result<CaddSource, PluginError> {
        // `add_file` runs before any header is read (CADD.pm:135).
        check_table_or_store(std::path::Path::new(path))?;
        let header_lines = header_lines(std::path::Path::new(path))?;
        if let Some(assembly) = requested {
            // The assembly's name anywhere in a header line, as the plugin greps
            // for it (CADD.pm:151).
            if !header_lines
                .iter()
                .any(|line| line.contains(assembly.as_str()))
            {
                return Err(PluginError::Init(format!(
                    "{NAME}: the run's assembly is {assembly} but the data file '{path}' does \
                     not name {assembly} in its header. Supply the {assembly} CADD file."
                )));
            }
        }
        // The last header line opening `#Chrom` names the columns (CADD.pm:159-167).
        let columns: Vec<&str> = header_lines
            .iter()
            .rev()
            .find(|line| line.starts_with("#Chrom"))
            .ok_or_else(|| {
                PluginError::Init(format!(
                    "{NAME}: '#Chrom' was not found in the header of '{path}'"
                ))
            })?
            .split('\t')
            .collect();
        let column = |names: &[&str]| {
            names
                .iter()
                .find_map(|name| columns.iter().position(|column| column == name))
        };
        let phred_col = column(&PHRED_COLUMNS);
        let raw_col = column(&RAW_COLUMNS);
        if phred_col.is_none() && raw_col.is_none() {
            return Err(PluginError::Init(format!(
                "{NAME}: {path} does not have a known column to be included"
            )));
        }
        // `$data{"Pos"} || $data{"Start"}` and `$data{"Alt"} || $data{"Type"}`
        // (CADD.pm:266, 268): a file with neither can match no row.
        let position_col = column(&["Pos", "Start"]).ok_or_else(|| {
            PluginError::Init(format!(
                "{NAME}: {path} has neither a Pos nor a Start column"
            ))
        })?;
        let alt_col = column(&["Alt"]);
        let type_col = column(&["Type"]);
        if alt_col.is_none() && type_col.is_none() {
            return Err(PluginError::Init(format!(
                "{NAME}: {path} has neither an Alt nor a Type column"
            )));
        }
        let ref_col = column(&["Ref"]);
        let end_col = column(&["End"]);

        let indexed = find_index_path(std::path::Path::new(path)).is_ok();
        let config = indexed_table_config(std::path::Path::new(path), |index_path| {
            format!(
                "index {} carries no tabix header, so its records cannot be mapped to \
                 chromosomes; index the file with tabix",
                index_path.display()
            )
        })?;
        // With no index the columns above are the ones the store beside the
        // file was built on, whatever they are. The store is opened on the
        // columns the header names, so one built on other columns is refused
        // by the flags that rebuild it.
        let config = if indexed {
            config
        } else {
            header_named_config(
                config.file_path,
                position_col,
                column(&["Pos"]).is_none(),
                end_col,
            )
        };
        // The last score cell keeps its trailing whitespace, which Ensembl VEP's
        // `split /\t/` leaves in the value it prints (CADD.pm:261, 274).
        let mut store = open_best_store_with(
            std::path::Path::new(path),
            config,
            requested,
            TabixAnnotator::keeping_trailing_whitespace,
        )?;

        // The index names no allele column, so the frame check is told the
        // header's `Ref`.
        if let Some(frame) = store.mitochondrial_frame_mut() {
            frame.set_allele_column(ref_col);
        }

        debug!(
            file = path,
            position_col,
            ?ref_col,
            ?alt_col,
            ?type_col,
            ?phred_col,
            ?raw_col,
            "opened CADD source"
        );

        Ok(CaddSource {
            store,
            path: path.to_string(),
            position_col,
            end_col,
            ref_col,
            alt_col,
            type_col,
            phred_col,
            raw_col,
        })
    }

    /// The side of the match a variant presents (CADD.pm:207-225): a sequence
    /// variant's forward-strand alleles at its start, when its alternate allele
    /// is bases of `ACGT` or `-` (:214); a structural variant's class word
    /// (`INS`, `DEL` or `DUP`, :75-79) against a `-` reference one base before
    /// its start (:220-224). `None` for a variant the plugin does not annotate:
    /// a structural variant with a sequence file loaded and no `force_annotate`
    /// (:218), or of a class the files do not score.
    fn variant_side<'v>(&self, variant: &'v InputVariant) -> Option<VariantSide<'v>> {
        if variant.is_structural || variant.variant_class.is_structural() {
            if self.sequence_file_named && !self.force_annotate {
                return None;
            }
            let class: &'static [u8] = match variant.variant_class.so_term() {
                "insertion" => b"INS",
                "deletion" => b"DEL",
                "duplication" => b"DUP",
                _ => return None,
            };
            return Some(VariantSide {
                reference: Cow::Borrowed(b"-"),
                alternate: Cow::Borrowed(class),
                position: variant.start.checked_sub(1)?,
                case: CaseRule::AsWritten,
            });
        }
        let (reference, alternate) = variant.forward_strand_alleles();
        if alternate.is_empty()
            || !alternate
                .iter()
                .all(|base| matches!(base, b'A' | b'C' | b'G' | b'T' | b'-'))
        {
            return None;
        }
        Some(VariantSide {
            reference,
            alternate,
            position: variant.start,
            case: CaseRule::EitherCase,
        })
    }
}

/// The row's side of the match, as `parse_data` builds it (CADD.pm:257-298)
/// and `get_matched_variant_alleles` reads it: the position column as written
/// ([`KeyPosition::Cell`], compared as text until a shared base moves it, then
/// as the number `$s + length($ref)` (:278) made of it: `30881493.0` moved once
/// is 30881494, `30881493.5` is no position), `-` for a missing or empty
/// reference (:267), the `Alt` cell when it is Perl-true and the `Type` cell
/// otherwise (:268), and the first base dropped from both alleles with the
/// position moved past it when they differ in length and share it (:277-289).
/// A row whose `End` cell is Perl-true skips that arithmetic, so its move is
/// `$s++` on a cell nothing has numified ([`allele_match::perl_increment`]:
/// `0100` steps to the text `0101`).
///
/// `Ok(None)` for a row whose position is Perl-false after the move, which
/// steps `0` and the empty cell to 1 when it happens: the match returns nothing
/// for such a row before it reads an allele (Sequence.pm:1173). `Err` for a row
/// whose alternate allele is undefined, which `trim_sequences` dies on
/// (Sequence.pm:968, reached from :1230), so that Ensembl VEP's `run` dies with
/// the row: `split /\t/` (:261) leaves a cell at or past the row's last filled
/// one undefined ([`split_field_count`]), and `||` leaves the allele undefined
/// for an empty or `0` `Alt` cell when the file has no `Type` column or the
/// `Type` cell is undefined too. An empty cell before a filled one is the
/// defined empty string: a CADD-SV row with an empty `Type` cell is compared
/// and dies nowhere, as is a row of a file with both columns whose `Alt` cell
/// is empty and whose `Type` cell is defined.
fn row_side<'r>(
    record: &'r TabixRecord,
    source: &CaddSource,
    case: CaseRule,
) -> Result<Option<RowSide<'r>>, PluginError> {
    // `$data{"Pos"} || $data{"Start"}` (:266): a Perl-false cell, `0` or empty,
    // and a missing one alike are the undefined value, which `!$b->{pos}`
    // (Sequence.pm:1173) passes the row over on, unless the shared-base move
    // (:283) steps it, as `$s++` steps `0` and the undefined value to 1.
    let position_text = record
        .columns
        .get(source.position_col)
        .map_or("", String::as_str);
    let mut position = Cow::Borrowed(position_text);
    let mut moved = 0;
    let present = split_field_count(&record.columns);
    let cell = |column: Option<usize>| {
        let column = column?;
        record
            .columns
            .get(column)
            .filter(|_| column < present)
            .map(String::as_str)
    };
    // `$data{"Ref"} || "-"` (:267): a reference cell that is missing, empty or
    // `0` is `-`.
    let mut reference: &[u8] = cell(source.ref_col)
        .filter(|cell| truthy(cell))
        .map_or(b"-", str::as_bytes);
    let alternate = match (cell(source.alt_col), cell(source.type_col)) {
        (Some(alt), _) if truthy(alt) => alt,
        (_, Some(kind)) => kind,
        // An undefined allele never moves the position, so a Perl-false one
        // returns the match before `trim_sequences` dies on the allele.
        (_, None) if !truthy(position_text) => return Ok(None),
        (_, None) => {
            return Err(allele_match::undefined_allele_row(
                NAME,
                &source.path,
                record,
            ))
        }
    };
    let mut alternate: &[u8] = alternate.as_bytes();
    if let Some((stripped_reference, stripped_alternate)) =
        strip_shared_first_base(reference, alternate, case)
    {
        if cell(source.end_col).is_some_and(truthy) {
            position = Cow::Owned(allele_match::perl_increment(position_text));
        } else {
            moved = 1;
        }
        reference = stripped_reference;
        alternate = stripped_alternate;
    }
    if moved == 0 && !truthy(&position) {
        return Ok(None);
    }
    Ok(Some(RowSide {
        reference,
        alternate,
        position,
        moved,
    }))
}

/// The first row of `records`, the rows of each file in file order, whose
/// alleles match `forms`, with the file it comes from: `get_data` hands `run`
/// the rows of every file in order and `run` returns at the first match
/// (CADD.pm:227, 237-253). Every row before the match is put to
/// `get_matched_variant_alleles`, so the first of them whose alternate allele
/// is undefined ([`row_side`]) is the `Err`: the match dies on it
/// (Sequence.pm:968), `run` dies with it and the variant takes no field. A row
/// after the match is never read, and a row without a position is passed over
/// (Sequence.pm:1173).
fn first_match<'a>(
    sources: &'a [CaddSource],
    records: &[Vec<&'a TabixRecord>],
    forms: &TrimmedForms<'_>,
    case: CaseRule,
) -> Result<Option<(&'a CaddSource, &'a TabixRecord)>, PluginError> {
    for (source, rows) in sources.iter().zip(records) {
        for &record in rows {
            if let Some(row) = row_side(record, source, case)? {
                if forms.matches(&row.alleles().trimmed_forms(case), case) {
                    return Ok(Some((source, record)));
                }
            }
        }
    }
    Ok(None)
}

/// Ensembl VEP's warning for a region answering more than [`RECORD_CAP`] rows
/// without `force_annotate` (CADD.pm:231-233).
fn record_cap_warning(chr: &str, start: u64, end: u64, count: usize) -> String {
    format!(
        "too many match found ({count}) for CADD variant with location {chr}_{start}_{end}. \
         No CADD annotation will be made. Make sure you are not using SNVs/Indels CADD \
         annotation file with structural variant as input. If you still want to annotate \
         please use force_annotate=1."
    )
}

impl crate::traits::sealed::Sealed for CaddPlugin {}

impl BuiltinPlugin for CaddPlugin {
    fn name(&self) -> &str {
        NAME
    }

    /// The Ensembl VEP module declares `['Feature', 'Intergenic']`: every row,
    /// the intergenic row included.
    fn feature_types(&self) -> &[RowKind] {
        RowKind::ALL
    }

    /// The Ensembl VEP module declares `['VariationFeature',
    /// 'StructuralVariationFeature']` (CADD.pm:186-188): structural variants
    /// are scored too, from the structural-variant file (:218-224).
    fn annotates_structural_variants(&self) -> bool {
        true
    }

    /// The fields the loaded files carry, each described as Ensembl VEP's
    /// plugin describes it (CADD.pm:81-98, 174); both fields before any file is
    /// loaded.
    fn header_info(&self) -> Vec<(String, String)> {
        let (phred, raw) = if self.sources.is_empty() {
            (true, true)
        } else {
            (
                self.sources.iter().any(|s| s.phred_col.is_some()),
                self.sources.iter().any(|s| s.raw_col.is_some()),
            )
        };
        let mut fields = Vec::new();
        if phred {
            fields.push((
                PHRED_FIELD.to_string(),
                format!("{PHRED_DESCRIPTION} {NON_COMMERCIAL_USE_CLAUSE}"),
            ));
        }
        if raw {
            fields.push((
                RAW_FIELD.to_string(),
                format!("{RAW_DESCRIPTION} {NON_COMMERCIAL_USE_CLAUSE}"),
            ));
        }
        fields
    }

    fn set_assembly(&mut self, assembly: Option<Assembly>) {
        self.assembly = assembly;
    }

    fn set_reference_fasta(&mut self, fasta: Option<std::sync::Arc<vep_fasta::IndexedFasta>>) {
        for source in &mut self.sources {
            source.store.set_reference_fasta(fasta.clone());
        }
    }

    fn init(&mut self, params: &[String]) -> Result<(), PluginError> {
        let parsed = PluginParams::parse_any_key(
            NAME,
            Grammar::PositionalOrKeyValue {
                flags: &POSITIONAL_FLAGS,
            },
            params,
        )?;
        let mut files: Vec<&str> = Vec::new();
        if parsed.keyed().is_empty() {
            // Every token but a `1` is a file; a trailing `1` is force_annotate
            // (CADD.pm:119-120).
            files.extend(parsed.positional().iter().map(String::as_str));
            self.force_annotate = parsed.trailing_flag() == Some("1");
        } else {
            // Every key but `force_annotate` names a file, `snv` and `indels`
            // a sequence file (CADD.pm:124-129).
            for (key, value) in parsed.keyed() {
                match key.as_str() {
                    "force_annotate" => self.force_annotate = truthy(value),
                    "snv" | "indels" => {
                        self.sequence_file_named = true;
                        files.push(value);
                    }
                    _ => files.push(value),
                }
            }
        }
        if files.is_empty() {
            return Err(PluginError::Init(format!(
                "{NAME}: no CADD files specified; give snv=, indels= or sv= a path, or list \
                 the files by position"
            )));
        }
        if self.force_annotate {
            warn!("{FORCE_ANNOTATE_WARNING}");
        }
        for path in files {
            self.sources.push(Self::open_source(path, self.assembly)?);
        }

        debug!(
            sources = self.sources.len(),
            force_annotate = self.force_annotate,
            sequence_file_named = self.sequence_file_named,
            "CADD plugin initialized"
        );

        Ok(())
    }

    fn prefetch(&self, variants: &[InputVariant]) -> Result<Option<PrefetchData>, PluginError> {
        if self.sources.is_empty() || variants.is_empty() {
            return Ok(None);
        }

        let regions: Vec<(&str, u64, u64)> = variants
            .iter()
            .filter_map(|variant| {
                let (start, end) = self.variant_side(variant)?.region(variant);
                Some((table_contig(variant, self.assembly), start, end))
            })
            .collect();

        // The files are read side by side; the first failure in file order is
        // the batch's.
        let batch_results =
            map_in_order(&self.sources, |source| source.store.query_batch(&regions))
                .into_iter()
                .collect::<Result<Vec<BatchQueryResult>, PluginError>>()?;

        Ok(Some(PrefetchData::new(batch_results)))
    }

    /// Scores each variant from the first matching row of its region
    /// (`first_match`). A variant whose rows hold, before that row, one whose
    /// alternate allele is undefined is the variant Ensembl VEP's `run` dies on
    /// (CADD.pm:237-253, Sequence.pm:968).
    fn annotate(
        &self,
        variants: &mut [InputVariant],
        data: PrefetchData,
    ) -> Result<(), PluginError> {
        let batch_results: Vec<BatchQueryResult> = data.downcast()?;

        let mut failure = None;
        for variant in variants.iter_mut() {
            let Some(side) = self.variant_side(variant) else {
                continue;
            };
            let (start, end) = side.region(variant);
            let records: Vec<Vec<&TabixRecord>> = batch_results
                .iter()
                .map(|batch| batch.lookup(table_contig(variant, self.assembly), start, end))
                .collect();
            // The cap counts every row the files answered for the region, as
            // the plugin counts its fetched rows (CADD.pm:230).
            let answered: usize = records.iter().map(Vec::len).sum();
            if answered > RECORD_CAP && !self.force_annotate {
                warn!(
                    input_record = variant.input_record,
                    "{}",
                    record_cap_warning(&variant.original_chr, side.position, variant.end, answered)
                );
                continue;
            }
            // The first row of the first file that matches writes the fields its
            // file carries (CADD.pm:237-253, 270-275): every column of the
            // header is a key of the row's hash (:263), so a cell the row lacks
            // is the undefined value, the empty value here, and an empty cell
            // before a filled one the empty string, named in the record's
            // `plugin_empty_strings`.
            let forms = side.alleles().trimmed_forms(side.case);
            match first_match(&self.sources, &records, &forms, side.case) {
                Ok(Some((source, record))) => {
                    let present = split_field_count(&record.columns);
                    for (field, col) in
                        [(PHRED_FIELD, source.phred_col), (RAW_FIELD, source.raw_col)]
                    {
                        let Some(col) = col else {
                            continue;
                        };
                        let cell = record.columns.get(col).filter(|_| col < present);
                        variant.set_plugin_value(
                            field.to_string(),
                            cell.cloned().unwrap_or_default(),
                            cell.is_some(),
                            self.feature_types(),
                        );
                    }
                }
                Ok(None) => {}
                Err(error) => {
                    failure.get_or_insert(error);
                }
            }
        }

        match failure {
            Some(error) => Err(error),
            None => Ok(()),
        }
    }

    fn run_batch(&self, variants: &mut [InputVariant]) -> Result<(), PluginError> {
        // Combined prefetch and annotate for callers that skip `prefetch`.
        if let Some(data) = self.prefetch(variants)? {
            self.annotate(variants, data)?;
        }
        Ok(())
    }

    fn run(
        &self,
        _consequence: &TranscriptConsequence,
        _variant: &InputVariant,
    ) -> Result<IndexMap<String, String>, PluginError> {
        // Variant-level plugin: annotation happens in `run_batch`.
        Ok(IndexMap::new())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tabix::TabixConfig;
    use crate::test_fixtures::{write_bgzf_bytes, write_tabix_fixture, IndexSpec};
    use std::collections::{BTreeMap, HashMap};
    use std::path::{Path, PathBuf};
    use vep_core::coordinate::Strand;
    use vep_core::variant::VariantClass;

    const GRCH38_MARKER: &str = "##CADD GRCh38-v1.7 made-up fixture, every score invented";
    const SNV_HEADER: &str = "#Chrom\tPos\tRef\tAlt\tRawScore\tPHRED";
    const SV_MARKER: &str = "## CADD-SV GRCh38 made-up fixture, every score invented";
    const SV_HEADER: &str =
        "#Chrom\tStart\tEnd\tType\tName\tCADD-SV_PHRED-score\tCADD-SV_Raw-score";

    fn args(tokens: &[&str]) -> Vec<String> {
        tokens.iter().map(|t| t.to_string()).collect()
    }

    fn init_message(result: Result<(), PluginError>) -> String {
        match result {
            Err(PluginError::Init(message)) => message,
            other => panic!("expected PluginError::Init, got {other:?}"),
        }
    }

    /// Build a made-up CADD fixture with a two-line preamble: a comment line
    /// naming the build, which the plugin matches against the run's assembly
    /// (CADD.pm:151-157), then the `#Chrom` column-name row (CADD.pm:161-167).
    fn cadd_fixture(dir: &Path, name: &str, assembly_token: &str) -> PathBuf {
        let provenance =
            format!("##CADD {assembly_token}-v1.7 made-up fixture, every score invented");
        let lines = vec![
            provenance.as_str(),
            SNV_HEADER,
            "21\t5030000\tA\tC\t0.345678\t5.432",
            "21\t5030000\tA\tG\t0.456789\t6.543",
            "21\t5030001\tC\tA\t0.234567\t4.321",
        ];
        write_tabix_fixture(dir, name, &lines, IndexSpec::standard())
    }

    /// A made-up CADD-SV fixture with the `Start` and `Type` columns the plugin
    /// reads (CADD.pm:266, 268), indexed with `tabix -s 1 -b 2 -e 3 -0`: a
    /// deletion and a duplication over one span, a deletion one base
    /// further along, and an insertion inside that span.
    fn sv_fixture(dir: &Path) -> PathBuf {
        let lines = [
            SV_MARKER,
            SV_HEADER,
            "21\t5028999\t5047000\tDEL\t1000G_DEL_chr21:5028999-5047000\t12.3\t0.45",
            "21\t5028999\t5047000\tDUP\t1000G_DUP_chr21:5028999-5047000\t3.1\t-0.2",
            "21\t5029000\t5047000\tDEL\t1000G_DEL_chr21:5029000-5047000\t20.0\t1.5",
            "21\t5030100\t5030101\tINS\t1000G_INS_chr21:5030100-5030101\t5.5\t0.3",
        ];
        write_tabix_fixture(
            dir,
            "synthetic_cadd_sv.tsv.gz",
            &lines,
            IndexSpec::standard().with_span_cols(1, 2, 3).zero_based(),
        )
    }

    /// A structural variant as the VCF parser builds one from `POS=<start - 1>`,
    /// `END=<end>`: its start the base after the anchor.
    fn structural(class: VariantClass, start: u64, end: u64) -> InputVariant {
        let mut variant =
            InputVariant::new("21".into(), start, end, b"N".to_vec(), b"<SV>".to_vec());
        variant.variant_class = class;
        variant.is_structural = true;
        variant
    }

    fn scores(variant: &InputVariant) -> (Option<&str>, Option<&str>) {
        (
            variant.plugin_data.get(PHRED_FIELD).map(String::as_str),
            variant.plugin_data.get(RAW_FIELD).map(String::as_str),
        )
    }

    #[test]
    fn new_plugin_is_uninitialized() {
        let plugin = CaddPlugin::new();
        assert_eq!(plugin.name(), "CADD");
        assert!(plugin.sources.is_empty());
        assert!(!plugin.force_annotate);
        assert!(!plugin.sequence_file_named);
    }

    /// The header lines are read as the bytes `tabix -h` prints and the plugin
    /// greps (CADD.pm:151-167): a byte that is not UTF-8 (Latin-1 `é`, 0xE9)
    /// in the provenance line leaves the assembly marker and the column names
    /// readable, the byte written as U+FFFD.
    #[test]
    fn header_lines_that_are_not_utf8_are_read() {
        let dir = tempfile::tempdir().unwrap();
        let path = write_bgzf_bytes(
            dir.path(),
            "synthetic_snv_score.tsv.gz",
            b"##CADD GRCh38-v1.7 made-up fixtur\xE9\n#Chrom\tPos\tRef\tAlt\tRawScore\tPHRED\n21\t5030000\tA\tC\t0.345678\t5.432\n",
        );
        assert_eq!(
            header_lines(&path).unwrap(),
            [
                "##CADD GRCh38-v1.7 made-up fixtur\u{FFFD}",
                "#Chrom\tPos\tRef\tAlt\tRawScore\tPHRED"
            ]
        );
    }

    /// `die "No CADD files specified"` (CADD.pm:134): no parameter, and a
    /// positional list of the flag alone, name no file.
    #[test]
    fn init_requires_a_file() {
        for tokens in [&[][..], &["1"][..], &["force_annotate=1"][..]] {
            let mut plugin = CaddPlugin::new();
            let message = init_message(plugin.init(&args(tokens)));
            assert!(
                message.starts_with("CADD: no CADD files specified"),
                "{tokens:?}: {message}"
            );
        }
    }

    /// The header carries Ensembl VEP's descriptions with the non-commercial
    /// clause (CADD.pm:81-100, 174), and after loading only the fields the files
    /// carry.
    #[test]
    fn header_info_carries_the_non_commercial_clause() {
        let plugin = CaddPlugin::new();
        assert_eq!(
            plugin.header_info(),
            vec![
                (
                    "CADD_PHRED".to_string(),
                    "PHRED-like scaled CADD score. CADD is only available here for \
                     non-commercial use. See CADD website for more information."
                        .to_string()
                ),
                (
                    "CADD_RAW".to_string(),
                    "Raw CADD score. CADD is only available here for non-commercial use. See \
                     CADD website for more information."
                        .to_string()
                ),
            ]
        );

        let dir = tempfile::tempdir().unwrap();
        let lines = [
            GRCH38_MARKER,
            "#Chrom\tPos\tRef\tAlt\tPHRED",
            "21\t5030000\tA\tG\t6.543",
        ];
        let path = write_tabix_fixture(
            dir.path(),
            "synthetic_phred_only.tsv.gz",
            &lines,
            IndexSpec::standard(),
        );
        let mut plugin = CaddPlugin::new();
        plugin.init(&[path.display().to_string()]).unwrap();
        let fields: Vec<String> = plugin.header_info().into_iter().map(|(f, _)| f).collect();
        assert_eq!(fields, ["CADD_PHRED"]);
    }

    /// A file with none of the score columns is refused (CADD.pm:177), as is one
    /// whose header has no `#Chrom` line (CADD.pm:167).
    #[test]
    fn file_without_a_score_column_or_chrom_header_is_an_init_error() {
        let dir = tempfile::tempdir().unwrap();
        let lines = [
            GRCH38_MARKER,
            "#Chrom\tPos\tRef\tAlt\tOther",
            "21\t5030000\tA\tG\t1",
        ];
        let path = write_tabix_fixture(
            dir.path(),
            "synthetic_no_score.tsv.gz",
            &lines,
            IndexSpec::standard(),
        );
        let message = init_message(CaddPlugin::new().init(&[path.display().to_string()]));
        assert!(
            message.ends_with("does not have a known column to be included"),
            "{message}"
        );

        let lines = [
            GRCH38_MARKER,
            "#Chr\tPos\tRef\tAlt\tPHRED",
            "21\t5030000\tA\tG\t1",
        ];
        let path = write_tabix_fixture(
            dir.path(),
            "synthetic_no_chrom.tsv.gz",
            &lines,
            IndexSpec::standard(),
        );
        let message = init_message(CaddPlugin::new().init(&[path.display().to_string()]));
        assert!(
            message.contains("'#Chrom' was not found in the header"),
            "{message}"
        );
    }

    /// `add_file` checks every file for `<file>.tbi` before a header is read
    /// (CADD.pm:135): a table with no index is refused naming the missing
    /// `.tbi`, not for its header, which has no `#Chrom` line.
    #[test]
    fn file_without_its_tbi_is_refused_before_its_header_is_read() {
        let dir = tempfile::tempdir().unwrap();
        let lines = [
            GRCH38_MARKER,
            "#Chr\tPos\tRef\tAlt\tOther",
            "21\t5030000\tA\tG\t1",
        ];
        let path = write_tabix_fixture(
            dir.path(),
            "synthetic_unindexed.tsv.gz",
            &lines,
            IndexSpec::standard(),
        );
        std::fs::remove_file(dir.path().join("synthetic_unindexed.tsv.gz.tbi")).unwrap();
        let message = init_message(CaddPlugin::new().init(&[path.display().to_string()]));
        assert!(
            message.contains("synthetic_unindexed.tsv.gz.tbi not found"),
            "{message}"
        );
    }

    /// A `.vpd` store beside the table is read in place of tabix, so the table
    /// needs no `<file>.tbi`: with its index named on the stem alone (found by
    /// htslib for the column layout, refused by Ensembl VEP's check) the plugin
    /// loads and scores from the store; without the store the same layout is
    /// refused naming the missing `.tbi`. The store holds the made-up table's own row,
    /// bound to the table's content.
    #[test]
    fn a_vpd_store_beside_the_table_needs_no_tbi() {
        use crate::binary_store::{BinaryStoreWriter, ContentIdentity};
        use std::sync::Arc;

        let dir = tempfile::tempdir().unwrap();
        let path = cadd_fixture(dir.path(), "synthetic_cadd.tsv.gz", "GRCh38");
        std::fs::rename(
            dir.path().join("synthetic_cadd.tsv.gz.tbi"),
            dir.path().join("synthetic_cadd.tsv.tbi"),
        )
        .unwrap();
        let message = init_message(CaddPlugin::new().init(&[path.display().to_string()]));
        assert!(
            message.contains("synthetic_cadd.tsv.gz.tbi not found"),
            "{message}"
        );

        let header: Vec<String> = SNV_HEADER
            .trim_start_matches('#')
            .split('\t')
            .map(String::from)
            .collect();
        let record = TabixRecord {
            chr: "21".to_string(),
            start: 5030000,
            end: 5030000,
            columns: ["21", "5030000", "A", "G", "0.456789", "6.543"]
                .map(String::from)
                .to_vec(),
            header: Some(Arc::new(header.clone())),
        };
        BinaryStoreWriter::new(dir.path().join("synthetic_cadd.tsv.vpd"))
            .write(
                &header,
                &[("21".to_string(), vec![record])],
                &TabixConfig::default(),
                &ContentIdentity::of_file(&path).unwrap(),
            )
            .unwrap();
        let mut plugin = CaddPlugin::new();
        plugin
            .init(&[path.display().to_string()])
            .expect("the store replaces tabix");
        let mut variants = vec![InputVariant::new(
            "21".into(),
            5030000,
            5030000,
            b"A".to_vec(),
            b"G".to_vec(),
        )];
        plugin.run_batch(&mut variants).expect("run_batch");
        assert_eq!(scores(&variants[0]), (Some("6.543"), Some("0.456789")));
    }

    #[test]
    fn run_batch_skips_when_no_sources() {
        let plugin = CaddPlugin::new();
        let mut variants = vec![];
        assert!(plugin.run_batch(&mut variants).is_ok());
    }

    /// `snv=` and `indels=` mark a sequence file (CADD.pm:128); `force_annotate`
    /// takes Perl's truth of its value (CADD.pm:131); every other key names a
    /// data file, which marks no sequence file (`:124-129`).
    #[test]
    fn parse_named_params() {
        let mut plugin = CaddPlugin::new();
        // The files do not exist, so `init` fails after the parameters are parsed.
        let result = plugin.init(&args(&[
            "snv=/path/to/synthetic_snv.tsv.gz",
            "indels=/path/to/synthetic_indels.tsv.gz",
            "force_annotate=1",
        ]));
        assert!(result.is_err());
        assert!(plugin.force_annotate);
        assert!(plugin.sequence_file_named);

        let mut plugin = CaddPlugin::new();
        assert!(plugin
            .init(&args(&[
                "sv=/path/to/synthetic_sv.tsv.gz",
                "force_annotate=0"
            ]))
            .is_err());
        assert!(!plugin.force_annotate);
        assert!(!plugin.sequence_file_named);

        let dir = tempfile::tempdir().unwrap();
        let snv = cadd_fixture(dir.path(), "synthetic_snv.tsv.gz", "GRCh38");
        let snv = snv.display().to_string();
        for key in ["file", "foo"] {
            let mut plugin = CaddPlugin::new();
            plugin.init(&args(&[&format!("{key}={snv}")])).unwrap();
            assert_eq!(plugin.sources.len(), 1, "{key}");
            assert!(!plugin.sequence_file_named, "{key}");
            assert!(!plugin.force_annotate, "{key}");
        }
        let message = init_message(
            CaddPlugin::new().init(&args(&["file=/path/to/synthetic_missing.tsv.gz"])),
        );
        assert!(
            message.contains("/path/to/synthetic_missing.tsv.gz"),
            "the value is opened as a file: {message}"
        );
    }

    /// Under the positional grammar a trailing `1` is `force_annotate` and never
    /// a file; a `1` elsewhere is dropped from the files and sets nothing
    /// (CADD.pm:119-120); positional files mark no sequence file (CADD.pm:116).
    #[test]
    fn positional_trailing_one_is_force_annotate() {
        let dir = tempfile::tempdir().unwrap();
        let snv = cadd_fixture(dir.path(), "synthetic_snv.tsv.gz", "GRCh38");
        let snv = snv.display().to_string();

        let mut plugin = CaddPlugin::new();
        plugin.init(&args(&[&snv, "1"])).unwrap();
        assert!(plugin.force_annotate);
        assert!(!plugin.sequence_file_named);
        assert_eq!(plugin.sources.len(), 1);

        let mut plugin = CaddPlugin::new();
        plugin.init(&args(&["1", &snv])).unwrap();
        assert!(!plugin.force_annotate);
        assert_eq!(plugin.sources.len(), 1);

        let mut plugin = CaddPlugin::new();
        plugin.init(&args(&[&snv])).unwrap();
        assert!(!plugin.force_annotate);
    }

    /// A GRCh37 CADD file in a GRCh38 run must be refused: Ensembl VEP's CADD.pm
    /// greps the header for the assembly and dies (CADD.pm:151-157).
    #[test]
    fn rejects_a_data_file_from_the_other_assembly() {
        let dir = tempfile::tempdir().unwrap();
        let g37 = cadd_fixture(dir.path(), "synthetic_cadd_grch37.tsv.gz", "GRCh37");

        let mut plugin = CaddPlugin::new();
        plugin.set_assembly(Some(Assembly::Grch38));
        let message = init_message(plugin.init(&[format!("snv={}", g37.display())]));
        assert_eq!(
            message,
            format!(
                "CADD: the run's assembly is GRCh38 but the data file '{}' does not name \
                 GRCh38 in its header. Supply the GRCh38 CADD file.",
                g37.display()
            )
        );
    }

    /// The matching case must load and annotate, so the guard is not simply
    /// rejecting everything.
    #[test]
    fn accepts_and_annotates_a_matching_assembly_file() {
        let dir = tempfile::tempdir().unwrap();
        let g38 = cadd_fixture(dir.path(), "synthetic_cadd_grch38.tsv.gz", "GRCh38");

        let mut plugin = CaddPlugin::new();
        plugin.set_assembly(Some(Assembly::Grch38));
        plugin
            .init(&[format!("snv={}", g38.display())])
            .expect("a matching-assembly file must load");

        let mut variants = vec![InputVariant::new(
            "21".into(),
            5030000,
            5030000,
            b"A".to_vec(),
            b"G".to_vec(),
        )];
        plugin.run_batch(&mut variants).expect("run_batch");

        // The A>G allele at 5030000 must match, not the A>C on the preceding
        // line, which also exercises the multi-record read path.
        assert_eq!(scores(&variants[0]), (Some("6.543"), Some("0.456789")));
    }

    /// Ensembl VEP reads the rows from two bases before the record's start
    /// (CADD.pm:227), and its base class refuses a start of 0 and cannot parse
    /// a start of -1 (BaseVepTabixPlugin.pm:354, 596-602), so a record at
    /// position 1 or 2 is never scored there; the window is clamped at
    /// position 1 and each record's own row is read.
    #[test]
    fn records_at_the_first_two_positions_of_a_chromosome_are_scored() {
        let dir = tempfile::tempdir().unwrap();
        let lines = [
            GRCH38_MARKER,
            SNV_HEADER,
            "21\t1\tC\tG\t9.984677\t39.801",
            "21\t2\tC\tG\t7.465190\t22.363",
            "21\t3\tC\tG\t1.239896\t3.312",
        ];
        let table = write_tabix_fixture(
            dir.path(),
            "synthetic_cadd_start.tsv.gz",
            &lines,
            IndexSpec::standard(),
        );

        let mut plugin = CaddPlugin::new();
        plugin.set_assembly(Some(Assembly::Grch38));
        plugin.init(&[format!("snv={}", table.display())]).unwrap();

        let mut variants = vec![
            InputVariant::new("21".into(), 1, 1, b"C".to_vec(), b"G".to_vec()),
            InputVariant::new("21".into(), 2, 2, b"C".to_vec(), b"G".to_vec()),
            InputVariant::new("21".into(), 3, 3, b"C".to_vec(), b"G".to_vec()),
        ];
        plugin.run_batch(&mut variants).expect("run_batch");

        assert_eq!(scores(&variants[0]), (Some("39.801"), Some("9.984677")));
        assert_eq!(scores(&variants[1]), (Some("22.363"), Some("7.465190")));
        assert_eq!(scores(&variants[2]), (Some("3.312"), Some("1.239896")));
    }

    /// Of three rows at position 1, a substitution and an insertion and a
    /// deletion written with their anchor base, each record takes the row of
    /// its own alleles (CADD.pm:239, 277-289): the insertion after base 1 and
    /// the deletion of base 2 both start at 2 and read the rows at 1 through
    /// the window clamped at position 1, and the substitution at 2 scores its
    /// own row and not a row at 1.
    #[test]
    fn indel_rows_at_position_one_are_read_for_records_starting_at_two() {
        let dir = tempfile::tempdir().unwrap();
        let lines = [
            GRCH38_MARKER,
            SNV_HEADER,
            "21\t1\tC\tG\t9.984677\t39.801",
            "21\t1\tC\tCA\t0.111111\t11.111",
            "21\t1\tCC\tC\t0.222222\t22.222",
            "21\t2\tC\tG\t7.465190\t22.363",
            "21\t3\tC\tG\t1.239896\t3.312",
        ];
        let table = write_tabix_fixture(
            dir.path(),
            "synthetic_cadd_indels_at_start.tsv.gz",
            &lines,
            IndexSpec::standard(),
        );

        let mut plugin = CaddPlugin::new();
        plugin.set_assembly(Some(Assembly::Grch38));
        plugin.init(&[format!("snv={}", table.display())]).unwrap();

        // The records as the VCF parser builds them from `21 1 C G`,
        // `21 1 C CA`, `21 1 CC C`, `21 2 C G` and `21 3 C G`: the anchor base
        // dropped, an insertion's start one past its end.
        let mut variants = vec![
            InputVariant::new("21".into(), 1, 1, b"C".to_vec(), b"G".to_vec()),
            InputVariant::new("21".into(), 2, 1, b"-".to_vec(), b"A".to_vec()),
            InputVariant::new("21".into(), 2, 2, b"C".to_vec(), b"-".to_vec()),
            InputVariant::new("21".into(), 2, 2, b"C".to_vec(), b"G".to_vec()),
            InputVariant::new("21".into(), 3, 3, b"C".to_vec(), b"G".to_vec()),
        ];
        plugin.run_batch(&mut variants).expect("run_batch");

        assert_eq!(
            scores(&variants[0]),
            (Some("39.801"), Some("9.984677")),
            "substitution at 1"
        );
        assert_eq!(
            scores(&variants[1]),
            (Some("11.111"), Some("0.111111")),
            "insertion after base 1"
        );
        assert_eq!(
            scores(&variants[2]),
            (Some("22.222"), Some("0.222222")),
            "deletion of base 2"
        );
        assert_eq!(
            scores(&variants[3]),
            (Some("22.363"), Some("7.465190")),
            "substitution at 2"
        );
        assert_eq!(
            scores(&variants[4]),
            (Some("3.312"), Some("1.239896")),
            "substitution at 3"
        );
    }

    /// Ensembl VEP dies when no header line names the assembly (CADD.pm:151-157),
    /// so a file without its provenance line is refused where the plugin was
    /// handed an assembly and loads where it was handed none, with nothing to
    /// check against.
    #[test]
    fn header_without_assembly_token_is_an_init_error() {
        let dir = tempfile::tempdir().unwrap();
        let lines = [SNV_HEADER, "21\t5030000\tA\tG\t0.456789\t6.543"];
        let path = write_tabix_fixture(
            dir.path(),
            "synthetic_cadd_nomarker.tsv.gz",
            &lines,
            IndexSpec::standard(),
        );

        let mut plugin = CaddPlugin::new();
        plugin.set_assembly(Some(Assembly::Grch38));
        let message = init_message(plugin.init(&[format!("snv={}", path.display())]));
        assert!(
            message.starts_with("CADD: the run's assembly is GRCh38 but the data file")
                && message.contains("does not name GRCh38 in its header"),
            "{message}"
        );

        let mut plugin = CaddPlugin::new();
        plugin
            .init(&[format!("snv={}", path.display())])
            .expect("no assembly means nothing to check");
    }

    /// A record given on the reverse strand scores the row of its forward-strand
    /// alleles: `get_matched_variant_alleles` (CADD.pm:239) reverse-complements
    /// them, so `T/C -1` at 43774726 takes the `A G` row's 25.0 and 2.5, as
    /// Ensembl VEP 116.3 with CADD.pm prints for this record over these two rows,
    /// and not the `T C` row's 88.0.
    #[test]
    fn reverse_strand_record_scores_the_complemented_row() {
        let dir = tempfile::tempdir().unwrap();
        let lines = [
            GRCH38_MARKER,
            SNV_HEADER,
            "21\t43774726\tA\tG\t2.5\t25.0",
            "21\t43774726\tT\tC\t8.0\t88.0",
        ];
        let path = write_tabix_fixture(
            dir.path(),
            "synthetic_cadd_strand.tsv.gz",
            &lines,
            IndexSpec::standard(),
        );
        let mut plugin = CaddPlugin::new();
        plugin.init(&[format!("snv={}", path.display())]).unwrap();

        let scores_on = |strand: Strand| {
            let mut variant = InputVariant::new(
                "21".into(),
                43774726,
                43774726,
                b"T".to_vec(),
                b"C".to_vec(),
            );
            variant.strand = strand;
            let mut variants = vec![variant];
            plugin.run_batch(&mut variants).unwrap();
            let data = &variants[0].plugin_data;
            (data.get(PHRED_FIELD).cloned(), data.get(RAW_FIELD).cloned())
        };
        assert_eq!(
            scores_on(Strand::Reverse),
            (Some("25.0".to_string()), Some("2.5".to_string()))
        );
        assert_eq!(
            scores_on(Strand::Forward),
            (Some("88.0".to_string()), Some("8.0".to_string()))
        );
    }

    /// Without an assembly there is nothing to contradict, so any file loads.
    #[test]
    fn no_requested_assembly_skips_the_guard() {
        let dir = tempfile::tempdir().unwrap();
        let g37 = cadd_fixture(dir.path(), "synthetic_cadd_grch37b.tsv.gz", "GRCh37");

        let mut plugin = CaddPlugin::new();
        plugin
            .init(&[format!("snv={}", g37.display())])
            .expect("no assembly means no mismatch to detect");
    }

    /// An indel row is written VCF-style with its anchor base, which `parse_data`
    /// drops with the position moved past it (CADD.pm:277-289); a record whose
    /// alleles share bases at both ends is trimmed to the change it names
    /// (Sequence.pm:1199-1208), so `ACG/ATG` takes the `C/T` row inside it.
    #[test]
    fn indel_rows_and_shared_bases_match_after_trimming() {
        let dir = tempfile::tempdir().unwrap();
        let lines = [
            GRCH38_MARKER,
            SNV_HEADER,
            "21\t1000\tA\tAT\t0.1\t1.5",
            "21\t1000\tAT\tA\t0.2\t2.5",
            "21\t1001\tC\tT\t0.3\t3.5",
        ];
        let path = write_tabix_fixture(
            dir.path(),
            "synthetic_cadd_indels.tsv.gz",
            &lines,
            IndexSpec::standard(),
        );
        let mut plugin = CaddPlugin::new();
        plugin
            .init(&[format!("indels={}", path.display())])
            .unwrap();

        let mut variants = vec![
            InputVariant::new("21".into(), 1001, 1000, b"-".to_vec(), b"T".to_vec()),
            InputVariant::new("21".into(), 1001, 1001, b"T".to_vec(), b"-".to_vec()),
            InputVariant::new("21".into(), 1000, 1002, b"ACG".to_vec(), b"ATG".to_vec()),
            InputVariant::new("21".into(), 1001, 1001, b"C".to_vec(), b"A".to_vec()),
        ];
        plugin.run_batch(&mut variants).unwrap();
        assert_eq!(
            scores(&variants[0]),
            (Some("1.5"), Some("0.1")),
            "insertion"
        );
        assert_eq!(scores(&variants[1]), (Some("2.5"), Some("0.2")), "deletion");
        assert_eq!(
            scores(&variants[2]),
            (Some("3.5"), Some("0.3")),
            "shared bases"
        );
        assert_eq!(scores(&variants[3]), (None, None), "another allele");
    }

    /// Ensembl VEP annotates an allele of `ACGT-` alone (CADD.pm:214), so an
    /// ambiguous allele of the record takes no score, and so would a
    /// lower-case one, which the parsers of both engines upper-case before
    /// any plugin runs (Parser.pm:554).
    #[test]
    fn lower_case_or_ambiguous_alternate_allele_is_not_scored() {
        let dir = tempfile::tempdir().unwrap();
        let path = cadd_fixture(dir.path(), "synthetic_cadd_case.tsv.gz", "GRCh38");
        let mut plugin = CaddPlugin::new();
        plugin.init(&[format!("snv={}", path.display())]).unwrap();
        let mut variants = vec![
            InputVariant::new("21".into(), 5030000, 5030000, b"a".to_vec(), b"g".to_vec()),
            InputVariant::new("21".into(), 5030000, 5030000, b"A".to_vec(), b"N".to_vec()),
        ];
        plugin.run_batch(&mut variants).unwrap();
        assert_eq!(scores(&variants[0]), (None, None));
        assert_eq!(scores(&variants[1]), (None, None));
    }

    /// A row whose allele cells are written in lower case scores the
    /// upper-case record it describes, where Ensembl VEP's key is compared as
    /// written (`get_matched_variant_alleles`, Sequence.pm:1207, 1232, on the
    /// cells CADD.pm:239-251 hands it): a substitution, and a deletion written
    /// with its anchor base, whose shared first base is found in either case
    /// (CADD.pm:279-289). Another allele at the row's position takes nothing.
    #[test]
    fn a_row_with_lower_case_alleles_scores_its_record() {
        let dir = tempfile::tempdir().unwrap();
        let lines = [
            "##CADD GRCh38-v1.7 made-up fixture, every score invented",
            SNV_HEADER,
            "21\t5030000\ta\tc\t0.345678\t5.432",
            "21\t5030010\tat\tA\t0.5\t7.0",
        ];
        let path = write_tabix_fixture(
            dir.path(),
            "synthetic_cadd_lower_case.tsv.gz",
            &lines,
            IndexSpec::standard(),
        );
        let mut plugin = CaddPlugin::new();
        plugin.init(&[format!("snv={}", path.display())]).unwrap();
        let mut variants = vec![
            InputVariant::new("21".into(), 5030000, 5030000, b"A".to_vec(), b"C".to_vec()),
            InputVariant::new("21".into(), 5030011, 5030011, b"T".to_vec(), b"-".to_vec()),
            InputVariant::new("21".into(), 5030000, 5030000, b"A".to_vec(), b"G".to_vec()),
        ];
        plugin.run_batch(&mut variants).unwrap();
        assert_eq!(scores(&variants[0]), (Some("5.432"), Some("0.345678")));
        assert_eq!(scores(&variants[1]), (Some("7.0"), Some("0.5")));
        assert_eq!(scores(&variants[2]), (None, None));
    }

    /// A CADD-SV row matches a structural variant by its class word and by its
    /// `Start`, one base before the variant's start (CADD.pm:220-224): the
    /// deletion over 5029000-5047000 takes the `DEL` row at 5028999 and not the
    /// `DUP` row at the same span or the `DEL` row one base along; a class the
    /// table lacks (CADD.pm:75-79) takes nothing. The plugin is handed the
    /// structural variants at all because CADD.pm:186-188 lists
    /// `StructuralVariationFeature` among its `variant_feature_types`.
    #[test]
    fn sv_file_scores_deletion_by_type_and_start() {
        let dir = tempfile::tempdir().unwrap();
        let sv = sv_fixture(dir.path());
        let mut plugin = CaddPlugin::new();
        plugin.set_assembly(Some(Assembly::Grch38));
        plugin.init(&[format!("sv={}", sv.display())]).unwrap();
        assert_eq!(plugin.sources[0].ref_col, None);
        assert!(plugin.annotates_structural_variants());

        let mut variants = vec![
            structural(VariantClass::StructuralDeletion, 5029000, 5047000),
            structural(VariantClass::Duplication, 5029000, 5047000),
            structural(VariantClass::StructuralDeletion, 5029002, 5047000),
            structural(VariantClass::Inversion, 5029000, 5047000),
            structural(VariantClass::StructuralInsertion, 5030101, 5030101),
        ];
        plugin.run_batch(&mut variants).unwrap();
        assert_eq!(
            scores(&variants[0]),
            (Some("12.3"), Some("0.45")),
            "deletion"
        );
        assert_eq!(
            scores(&variants[1]),
            (Some("3.1"), Some("-0.2")),
            "duplication"
        );
        assert_eq!(
            scores(&variants[2]),
            (None, None),
            "deletion two bases along"
        );
        assert_eq!(scores(&variants[3]), (None, None), "inversion");
        assert_eq!(
            scores(&variants[4]),
            (Some("5.5"), Some("0.3")),
            "insertion"
        );
    }

    /// A CADD-SV row whose class word is written in lower case scores no
    /// structural variant: the word is compared as written with the variant's
    /// own (`$SV_TERMS{$so_term}`, CADD.pm:223, in the key of Sequence.pm:1207,
    /// 1232), where a sequence variant's row is read in either letter case.
    /// The row in capitals over another span scores its deletion.
    #[test]
    fn a_class_word_in_lower_case_scores_no_structural_variant() {
        let dir = tempfile::tempdir().unwrap();
        let lines = [
            SV_MARKER,
            SV_HEADER,
            "21\t5028999\t5047000\tdel\tsynthetic_lower_case_word\t12.3\t0.45",
            "21\t5050099\t5060000\tDEL\tsynthetic_upper_case_word\t20.0\t1.5",
        ];
        let path = write_tabix_fixture(
            dir.path(),
            "synthetic_cadd_sv_lower_case.tsv.gz",
            &lines,
            IndexSpec::standard().with_span_cols(1, 2, 3).zero_based(),
        );
        let mut plugin = CaddPlugin::new();
        plugin.init(&[format!("sv={}", path.display())]).unwrap();
        let mut variants = vec![
            structural(VariantClass::StructuralDeletion, 5029000, 5047000),
            structural(VariantClass::StructuralDeletion, 5050100, 5060000),
        ];
        plugin.run_batch(&mut variants).unwrap();
        assert_eq!(scores(&variants[0]), (None, None), "lower-case word");
        assert_eq!(
            scores(&variants[1]),
            (Some("20.0"), Some("1.5")),
            "upper-case word"
        );
    }

    /// With a file under `snv=` or `indels=`, a structural variant is annotated
    /// only under `force_annotate` (CADD.pm:218); files given by position mark
    /// no sequence file (CADD.pm:116), so the variant is annotated without it.
    /// A sequence variant of the same batch is scored in every case.
    #[test]
    fn sv_input_is_skipped_with_snv_file_unless_forced() {
        let dir = tempfile::tempdir().unwrap();
        let snv = cadd_fixture(dir.path(), "synthetic_snv.tsv.gz", "GRCh38");
        let sv = sv_fixture(dir.path());
        let (snv, sv) = (snv.display().to_string(), sv.display().to_string());

        let run = |tokens: &[&str]| {
            let mut plugin = CaddPlugin::new();
            plugin.set_assembly(Some(Assembly::Grch38));
            plugin.init(&args(tokens)).unwrap();
            let mut variants = vec![
                structural(VariantClass::StructuralDeletion, 5029000, 5047000),
                InputVariant::new("21".into(), 5030000, 5030000, b"A".to_vec(), b"G".to_vec()),
            ];
            plugin.run_batch(&mut variants).unwrap();
            (
                variants[0].plugin_data.get(PHRED_FIELD).cloned(),
                variants[1].plugin_data.get(PHRED_FIELD).cloned(),
            )
        };
        let (snv_score, sv_score) = (Some("6.543".to_string()), Some("12.3".to_string()));
        assert_eq!(
            run(&[&format!("snv={snv}"), &format!("sv={sv}")]),
            (None, snv_score.clone())
        );
        assert_eq!(
            run(&[
                &format!("snv={snv}"),
                &format!("sv={sv}"),
                "force_annotate=1"
            ]),
            (sv_score.clone(), snv_score.clone())
        );
        assert_eq!(run(&[&snv, &sv]), (sv_score.clone(), snv_score));
        assert_eq!(run(&[&format!("sv={sv}")]), (sv_score, None));
    }

    /// A region the files answer with more than 100,000 rows is not annotated
    /// without `force_annotate`, with Ensembl VEP's warning (CADD.pm:230-235);
    /// with it, the matching row among them scores the variant.
    #[test]
    fn over_100000_records_without_force_is_a_warning_and_no_score() {
        let dir = tempfile::tempdir().unwrap();
        let sv = sv_fixture(dir.path());
        let variant = structural(VariantClass::StructuralDeletion, 5029000, 5200000);

        // The made-up rows a dense file would answer over the span, built as the batch
        // query would index them, the matching row first among them.
        let columns = |start: u64, class: &str| -> Vec<String> {
            vec![
                String::new(),
                start.to_string(),
                "5200000".to_string(),
                class.to_string(),
                "name".to_string(),
                "12.3".to_string(),
                "0.45".to_string(),
            ]
        };
        let mut by_start: BTreeMap<u64, Vec<TabixRecord>> = BTreeMap::new();
        for start in 5029000..5029000 + RECORD_CAP as u64 + 1 {
            by_start.insert(
                start,
                vec![TabixRecord {
                    chr: "21".into(),
                    start,
                    end: 5200000,
                    columns: columns(start - 1, if start == 5029000 { "DEL" } else { "DUP" }),
                    header: None,
                }],
            );
        }
        let records: HashMap<String, BTreeMap<u64, Vec<TabixRecord>>> =
            HashMap::from([("21".to_string(), by_start)]);

        let annotate = |tokens: &[&str]| {
            let mut plugin = CaddPlugin::new();
            plugin.init(&args(tokens)).unwrap();
            let mut variants = vec![variant.clone()];
            let data = PrefetchData::new(vec![BatchQueryResult::from_records(records.clone())]);
            plugin.annotate(&mut variants, data).unwrap();
            variants[0].plugin_data.get(PHRED_FIELD).cloned()
        };
        let sv = sv.display().to_string();
        assert_eq!(annotate(&[&format!("sv={sv}")]), None);
        assert_eq!(
            annotate(&[&format!("sv={sv}"), "force_annotate=1"]),
            Some("12.3".to_string())
        );
        assert_eq!(
            record_cap_warning("21", 5028999, 5200000, RECORD_CAP + 1),
            "too many match found (100001) for CADD variant with location 21_5028999_5200000. \
             No CADD annotation will be made. Make sure you are not using SNVs/Indels CADD \
             annotation file with structural variant as input. If you still want to annotate \
             please use force_annotate=1."
        );

        // The warning names the input record it is about, for a run that
        // holds a batch's warnings to place it with that record.
        let warnings = std::sync::Arc::new(crate::test_fixtures::Warnings::default());
        let mut plugin = CaddPlugin::new();
        plugin.init(&args(&[&format!("sv={sv}")])).unwrap();
        let mut variants = vec![variant.clone()];
        variants[0].input_record = 17;
        let data = PrefetchData::new(vec![BatchQueryResult::from_records(records.clone())]);
        tracing::subscriber::with_default(std::sync::Arc::clone(&warnings), || {
            plugin.annotate(&mut variants, data).unwrap();
        });
        assert_eq!(
            *warnings.0.lock().unwrap(),
            [record_cap_warning("21", 5028999, 5200000, RECORD_CAP + 1)]
        );
        assert_eq!(*warnings.1.lock().unwrap(), [17]);
    }

    /// CADD.pm:261 splits the line htslib hands over, the line end alone
    /// removed, and :274 copies each score cell as split: a trailing space stays
    /// in the last column of the test's header (`PHRED` in its SNV file, `CADD-SV_Raw-score`
    /// in its CADD-SV file) and a whitespace-only cell is a value of one space. A
    /// line ending in a tab leaves the value empty, which Ensembl VEP prints as
    /// `CADD_PHRED=` (its `split` drops the trailing empty field and the
    /// undefined value prints as nothing after the `=`).
    #[test]
    fn a_trailing_space_on_the_last_score_is_printed_as_ensembl_vep_prints_it() {
        let dir = tempfile::tempdir().unwrap();
        let snv = write_tabix_fixture(
            dir.path(),
            "synthetic_cadd_trailing.tsv.gz",
            &[
                GRCH38_MARKER,
                SNV_HEADER,
                "21\t5030000\tA\tG\t0.456789\t6.543 ",
                "21\t5030001\tC\tA\t0.234567\t ",
                "21\t5030002\tG\tT\t0.5\t",
            ],
            IndexSpec::standard(),
        );
        let mut plugin = CaddPlugin::new();
        plugin.init(&[format!("snv={}", snv.display())]).unwrap();
        let mut variants = vec![
            InputVariant::new("21".into(), 5030000, 5030000, b"A".to_vec(), b"G".to_vec()),
            InputVariant::new("21".into(), 5030001, 5030001, b"C".to_vec(), b"A".to_vec()),
            InputVariant::new("21".into(), 5030002, 5030002, b"G".to_vec(), b"T".to_vec()),
        ];
        plugin.run_batch(&mut variants).unwrap();
        assert_eq!(scores(&variants[0]), (Some("6.543 "), Some("0.456789")));
        assert_eq!(scores(&variants[1]), (Some(" "), Some("0.234567")));
        assert_eq!(scores(&variants[2]), (Some(""), Some("0.5")));

        let sv = write_tabix_fixture(
            dir.path(),
            "synthetic_cadd_sv_trailing.tsv.gz",
            &[
                SV_MARKER,
                SV_HEADER,
                "21\t5028999\t5047000\tDEL\t1000G_DEL_chr21:5028999-5047000\t12.3\t0.45 ",
                "21\t5030100\t5035000\tDUP\t1000G_DUP_chr21:5030100-5035000\t10.9\t ",
            ],
            IndexSpec::standard().with_span_cols(1, 2, 3).zero_based(),
        );
        let mut plugin = CaddPlugin::new();
        plugin.init(&[format!("sv={}", sv.display())]).unwrap();
        let mut variants = vec![
            structural(VariantClass::StructuralDeletion, 5029000, 5047000),
            structural(VariantClass::Duplication, 5030101, 5035000),
        ];
        plugin.run_batch(&mut variants).unwrap();
        assert_eq!(scores(&variants[0]), (Some("12.3"), Some("0.45 ")));
        assert_eq!(scores(&variants[1]), (Some("10.9"), Some(" ")));
    }

    /// `trim_sequences` in both orders (Sequence.pm:976-1030): each order trims
    /// both ends, the position moves past the bases trimmed from the start
    /// alone, and an allele trimmed away is `-`; a deletion inside a repeat has
    /// one form per order, and either form matches a row written that way.
    #[test]
    fn trimming_follows_trim_sequences() {
        let alleles = Alleles {
            reference: b"ACGT",
            alternate: b"AGT",
            position: KeyPosition::Number(100),
        };
        assert_eq!(
            alleles.trimmed(false, CaseRule::AsWritten),
            (&b"C"[..], &b"-"[..], KeyPosition::Number(101))
        );
        assert_eq!(
            alleles.trimmed(true, CaseRule::AsWritten),
            (&b"C"[..], &b"-"[..], KeyPosition::Number(101))
        );

        let repeat = Alleles {
            reference: b"ATA",
            alternate: b"A",
            position: KeyPosition::Number(100),
        };
        assert_eq!(
            repeat.trimmed(false, CaseRule::AsWritten),
            (&b"TA"[..], &b"-"[..], KeyPosition::Number(101))
        );
        assert_eq!(
            repeat.trimmed(true, CaseRule::AsWritten),
            (&b"AT"[..], &b"-"[..], KeyPosition::Number(100))
        );
        let forms = repeat.trimmed_forms(CaseRule::AsWritten);
        assert_eq!(
            forms.start_first,
            (&b"TA"[..], &b"-"[..], KeyPosition::Number(101))
        );
        assert_eq!(
            forms.end_first,
            Some((&b"AT"[..], &b"-"[..], KeyPosition::Number(100)))
        );
        let left_row = Alleles {
            reference: b"AT",
            alternate: b"-",
            position: KeyPosition::Number(100),
        };
        let right_row = Alleles {
            reference: b"TA",
            alternate: b"-",
            position: KeyPosition::Number(101),
        };
        assert!(forms.matches(
            &left_row.trimmed_forms(CaseRule::AsWritten),
            CaseRule::AsWritten
        ));
        assert!(forms.matches(
            &right_row.trimmed_forms(CaseRule::AsWritten),
            CaseRule::AsWritten
        ));

        let snv = Alleles {
            reference: b"A",
            alternate: b"G",
            position: KeyPosition::Number(100),
        };
        assert_eq!(snv.trimmed_forms(CaseRule::AsWritten).end_first, None);
        assert!(!snv.trimmed_forms(CaseRule::AsWritten).matches(
            &left_row.trimmed_forms(CaseRule::AsWritten),
            CaseRule::AsWritten
        ));
    }

    /// The error `run_batch` returns for a row Ensembl VEP's match dies on, with
    /// the plugin's fields absent from the failed variant.
    fn undefined_row_failure(result: Result<(), PluginError>, failed: &InputVariant) -> String {
        let message = match result {
            Err(PluginError::Run(message)) => message,
            other => panic!("expected PluginError::Run, got {other:?}"),
        };
        assert!(
            message.contains("has no reference or alternate allele cell"),
            "{message}"
        );
        assert!(failed.plugin_data.is_empty(), "{:?}", failed.plugin_data);
        message
    }

    /// A made-up CADD-SV table whose `DEL` row over 5028999-5047000 has a row of its
    /// three coordinate cells alone (`Type` undefined) `short_first` or after
    /// it, and an `INS` row outside that span.
    fn sv_fixture_with_a_short_row(dir: &Path, name: &str, short_first: bool) -> PathBuf {
        let deletion = "21\t5028999\t5047000\tDEL\t1000G_DEL_chr21:5028999-5047000\t12.3\t0.45";
        let insertion = "21\t5050100\t5050101\tINS\t1000G_INS_chr21:5050100-5050101\t5.5\t0.3";
        let short = "21\t5028999\t5047000";
        let lines = if short_first {
            [SV_MARKER, SV_HEADER, short, deletion, insertion]
        } else {
            [SV_MARKER, SV_HEADER, deletion, short, insertion]
        };
        write_tabix_fixture(
            dir,
            name,
            &lines,
            IndexSpec::standard().with_span_cols(1, 2, 3).zero_based(),
        )
    }

    /// The deletion the `DEL` row of [`sv_fixture_with_a_short_row`] scores and
    /// the insertion its `INS` row scores.
    fn sv_variants() -> Vec<InputVariant> {
        vec![
            structural(VariantClass::StructuralDeletion, 5029000, 5047000),
            structural(VariantClass::StructuralInsertion, 5050101, 5050101),
        ]
    }

    /// `run` puts every row of the region to `get_matched_variant_alleles` in
    /// file order and returns at the first match (CADD.pm:237-253), and the
    /// match dies on a row whose alternate allele is undefined
    /// (Sequence.pm:968, from :1230): the output line gets none of the plugin's
    /// fields. In a file with no `Type` column `$data{"Alt"} || $data{"Type"}` (:268)
    /// is undefined for an `Alt` cell past the row's last filled one (`split /\t/`,
    /// :261), and for an empty or `0` `Alt` cell as well. The variant whose
    /// rows hold such a row before its
    /// match fails, the allele no row matches included, whose loop reaches an
    /// undefined row after the other alleles' rows; a variant whose region holds
    /// no such row is scored.
    #[test]
    fn an_undefined_alt_before_the_matching_row_fails_the_variant() {
        let dir = tempfile::tempdir().unwrap();
        // Each made-up row below is written into a made-up table of its own.
        for (name, undefined_row) in [
            ("short", "21\t5030000\tA"),
            ("trailing_tab", "21\t5030000\tA\t"),
            ("empty_alt", "21\t5030000\tA\t\t0.5\t9.0"),
            ("zero_alt", "21\t5030000\tA\t0\t0.5\t9.0"),
        ] {
            let lines = [
                GRCH38_MARKER,
                SNV_HEADER,
                undefined_row,
                "21\t5030000\tA\tG\t0.456789\t6.543",
                "21\t5030005\tC\tA\t0.234567\t4.321",
            ];
            let path = write_tabix_fixture(
                dir.path(),
                &format!("synthetic_cadd_{name}.tsv.gz"),
                &lines,
                IndexSpec::standard(),
            );
            let mut plugin = CaddPlugin::new();
            plugin.init(&[format!("snv={}", path.display())]).unwrap();
            let mut variants = vec![
                InputVariant::new("21".into(), 5030000, 5030000, b"A".to_vec(), b"G".to_vec()),
                InputVariant::new("21".into(), 5030005, 5030005, b"C".to_vec(), b"A".to_vec()),
            ];
            let message = undefined_row_failure(plugin.run_batch(&mut variants), &variants[0]);
            assert!(
                message.contains(&format!("{} at 21:5030000", path.display())),
                "{name}: {message}"
            );
            assert_eq!(
                scores(&variants[1]),
                (Some("4.321"), Some("0.234567")),
                "{name}: the variant five bases along reads no undefined row"
            );
        }

        let path = write_tabix_fixture(
            dir.path(),
            "synthetic_cadd_short_last.tsv.gz",
            &[
                GRCH38_MARKER,
                SNV_HEADER,
                "21\t5030000\tA\tC\t0.345678\t5.432",
                "21\t5030000\tA\tG\t0.456789\t6.543",
                "21\t5030000\tA",
            ],
            IndexSpec::standard(),
        );
        let mut plugin = CaddPlugin::new();
        plugin.init(&[format!("snv={}", path.display())]).unwrap();
        let mut variants = vec![InputVariant::new(
            "21".into(),
            5030000,
            5030000,
            b"A".to_vec(),
            b"T".to_vec(),
        )];
        undefined_row_failure(plugin.run_batch(&mut variants), &variants[0]);
    }

    /// In a CADD-SV file the alternate allele is the `Type` cell (CADD.pm:268),
    /// undefined only past the row's last filled cell: a row of its three
    /// coordinate cells alone before the matching row fails the structural
    /// variant, and a structural variant outside the row's span is scored.
    #[test]
    fn a_type_cell_past_the_last_filled_one_fails_the_structural_variant() {
        let dir = tempfile::tempdir().unwrap();
        let path =
            sv_fixture_with_a_short_row(dir.path(), "synthetic_cadd_sv_short_before.tsv.gz", true);
        let mut plugin = CaddPlugin::new();
        plugin.init(&[format!("sv={}", path.display())]).unwrap();
        let mut variants = sv_variants();
        let message = undefined_row_failure(plugin.run_batch(&mut variants), &variants[0]);
        assert!(
            message.contains(&format!("{} at 21:5029000", path.display())),
            "{message}"
        );
        assert_eq!(scores(&variants[1]), (Some("5.5"), Some("0.3")));
    }

    /// `run` returns at the first match (CADD.pm:252), so a row after it is
    /// never put to the match: the alleles whose rows precede the undefined row
    /// are scored from a SNV file, and the deletion whose `DEL` row precedes the
    /// short row of a CADD-SV file is scored, as is the insertion outside its
    /// span.
    #[test]
    fn an_undefined_alt_after_the_matching_row_is_never_read() {
        let dir = tempfile::tempdir().unwrap();
        let snv = write_tabix_fixture(
            dir.path(),
            "synthetic_cadd_short_after.tsv.gz",
            &[
                GRCH38_MARKER,
                SNV_HEADER,
                "21\t5030000\tA\tC\t0.345678\t5.432",
                "21\t5030000\tA\tG\t0.456789\t6.543",
                "21\t5030000\tA",
            ],
            IndexSpec::standard(),
        );
        let mut plugin = CaddPlugin::new();
        plugin.init(&[format!("snv={}", snv.display())]).unwrap();
        let mut variants = vec![
            InputVariant::new("21".into(), 5030000, 5030000, b"A".to_vec(), b"G".to_vec()),
            InputVariant::new("21".into(), 5030000, 5030000, b"A".to_vec(), b"C".to_vec()),
        ];
        plugin.run_batch(&mut variants).unwrap();
        assert_eq!(scores(&variants[0]), (Some("6.543"), Some("0.456789")));
        assert_eq!(scores(&variants[1]), (Some("5.432"), Some("0.345678")));

        let sv =
            sv_fixture_with_a_short_row(dir.path(), "synthetic_cadd_sv_short_after.tsv.gz", false);
        let mut plugin = CaddPlugin::new();
        plugin.init(&[format!("sv={}", sv.display())]).unwrap();
        let mut variants = sv_variants();
        plugin.run_batch(&mut variants).unwrap();
        assert_eq!(scores(&variants[0]), (Some("12.3"), Some("0.45")));
        assert_eq!(scores(&variants[1]), (Some("5.5"), Some("0.3")));
    }

    /// An empty cell before a filled one is the defined empty string to
    /// `split /\t/` (CADD.pm:261): a CADD-SV row with an empty `Type` cell
    /// dies nowhere and matches no class, and in a file with both an `Alt` and a
    /// `Type` column an empty or `0` `Alt` cell takes the `Type` cell (:268),
    /// defined and matching no allele; the matching row after either scores
    /// the variant.
    #[test]
    fn an_empty_alt_cell_before_a_filled_one_is_a_defined_allele() {
        let dir = tempfile::tempdir().unwrap();
        let sv = write_tabix_fixture(
            dir.path(),
            "synthetic_cadd_sv_empty_type.tsv.gz",
            &[
                SV_MARKER,
                SV_HEADER,
                "21\t5028999\t5047000\t\t1000G_DEL_chr21:5028999-5047000\t1.0\t0.1",
                "21\t5028999\t5047000\tDEL\t1000G_DEL_chr21:5028999-5047000\t12.3\t0.45",
            ],
            IndexSpec::standard().with_span_cols(1, 2, 3).zero_based(),
        );
        let mut plugin = CaddPlugin::new();
        plugin.init(&[format!("sv={}", sv.display())]).unwrap();
        let mut variants = vec![structural(
            VariantClass::StructuralDeletion,
            5029000,
            5047000,
        )];
        plugin.run_batch(&mut variants).unwrap();
        assert_eq!(scores(&variants[0]), (Some("12.3"), Some("0.45")));

        let annotated = write_tabix_fixture(
            dir.path(),
            "synthetic_cadd_annotated.tsv.gz",
            &[
                GRCH38_MARKER,
                "#Chrom\tPos\tRef\tAlt\tType\tRawScore\tPHRED",
                "21\t5030000\tA\t\tSNV\t0.1\t1.0",
                "21\t5030000\tA\t0\tSNV\t0.2\t2.0",
                "21\t5030000\tA\tG\tSNV\t0.456789\t6.543",
            ],
            IndexSpec::standard(),
        );
        let mut plugin = CaddPlugin::new();
        plugin
            .init(&[format!("snv={}", annotated.display())])
            .unwrap();
        let mut variants = vec![InputVariant::new(
            "21".into(),
            5030000,
            5030000,
            b"A".to_vec(),
            b"G".to_vec(),
        )];
        plugin.run_batch(&mut variants).unwrap();
        assert_eq!(scores(&variants[0]), (Some("6.543"), Some("0.456789")));
    }

    /// Every column of the header is a key of a matched row's hash
    /// (CADD.pm:263, 270-275): a row short of its `PHRED` cell, or of both
    /// score cells, writes the keys with the undefined value, the empty value
    /// here, as does a trailing empty cell, which the `split /\t/` of `:261`
    /// drops; an empty `RawScore` cell before a filled `PHRED` writes the empty
    /// string, named in the record's `plugin_empty_strings`.
    #[test]
    fn a_matched_row_short_of_a_score_cell_writes_the_key_with_the_empty_value() {
        let dir = tempfile::tempdir().unwrap();
        let path = write_tabix_fixture(
            dir.path(),
            "synthetic_cadd_short_rows.tsv.gz",
            &[
                GRCH38_MARKER,
                SNV_HEADER,
                "21\t5030000\tA\tG\t0.456789",
                "21\t5030001\tA\tG",
                "21\t5030002\tA\tG\t\t6.543",
                "21\t5030003\tA\tG\t0.456789\t",
            ],
            IndexSpec::standard(),
        );
        let mut plugin = CaddPlugin::new();
        plugin.init(&[format!("snv={}", path.display())]).unwrap();
        let mut variants: Vec<InputVariant> = (5030000..=5030003)
            .map(|pos| InputVariant::new("21".into(), pos, pos, b"A".to_vec(), b"G".to_vec()))
            .collect();
        plugin.run_batch(&mut variants).unwrap();
        assert_eq!(scores(&variants[0]), (Some(""), Some("0.456789")));
        assert_eq!(scores(&variants[1]), (Some(""), Some("")));
        assert_eq!(scores(&variants[2]), (Some("6.543"), Some("")));
        assert_eq!(scores(&variants[3]), (Some(""), Some("0.456789")));
        let marks: Vec<&[String]> = variants
            .iter()
            .map(|variant| variant.plugin_empty_strings.as_slice())
            .collect();
        assert_eq!(marks, [&[][..], &[], &[RAW_FIELD.to_string()], &[]]);
    }

    /// The cells `parse_data` reads (CADD.pm:265-268) as the match takes them:
    /// a whitespace-only `Alt` cell is a defined allele; an empty, `0` or
    /// missing `Alt` cell with no `Type` cell is the undefined allele the match
    /// dies on; a `Type` cell is undefined past the last filled cell alone; a
    /// missing or empty `Ref` cell is `-` (:267); and a row without a position
    /// is passed over before any allele is read (Sequence.pm:1173), whatever
    /// its allele cells hold, unless the shared-base move steps it.
    #[test]
    fn row_side_reads_the_cells_as_parse_data_does() {
        let dir = tempfile::tempdir().unwrap();
        // Every row of this test is made up, built here from its cells.
        let record = |columns: &[&str]| TabixRecord {
            chr: "21".to_string(),
            start: 5030000,
            end: 5030000,
            columns: columns.iter().map(|c| c.to_string()).collect(),
            header: None,
        };
        fn alleles(
            row: Result<Option<RowSide<'_>>, PluginError>,
        ) -> Option<(Vec<u8>, Vec<u8>, String, u64)> {
            row.unwrap().map(|row| {
                (
                    row.reference.to_vec(),
                    row.alternate.to_vec(),
                    row.position.to_string(),
                    row.moved,
                )
            })
        }

        let snv = cadd_fixture(dir.path(), "synthetic_cadd_cells.tsv.gz", "GRCh38");
        let mut plugin = CaddPlugin::new();
        plugin.init(&[format!("snv={}", snv.display())]).unwrap();
        let source = &plugin.sources[0];
        assert_eq!(
            alleles(row_side(
                &record(&["21", "5030000", "A", " ", "0.5", "9.0"]),
                source,
                CaseRule::AsWritten
            )),
            Some((b"A".to_vec(), b" ".to_vec(), "5030000".to_string(), 0))
        );
        assert_eq!(
            alleles(row_side(
                &record(&["21", "5030000", "", "G", "0.5", "9.0"]),
                source,
                CaseRule::AsWritten
            )),
            Some((b"-".to_vec(), b"G".to_vec(), "5030000".to_string(), 0))
        );
        for columns in [
            &["21", "5030000", "A", "", "0.5", "9.0"][..],
            &["21", "5030000", "A", "0", "0.5", "9.0"][..],
            &["21", "5030000", "A"][..],
            &["21", "5030000", "A", ""][..],
        ] {
            let error = row_side(&record(columns), source, CaseRule::AsWritten).err();
            assert!(
                matches!(&error, Some(PluginError::Run(message)) if message.contains("at 21:5030000 has no reference or alternate allele cell")),
                "{columns:?}: {error:?}"
            );
        }
        assert_eq!(
            alleles(row_side(
                &record(&["21", "", "A"]),
                source,
                CaseRule::AsWritten
            )),
            None
        );
        assert_eq!(
            alleles(row_side(
                &record(&["21", "", "A", "G"]),
                source,
                CaseRule::AsWritten
            )),
            None
        );
        // A `0` reference is `-` (`$data{"Ref"} || "-"`, :267); a `0` position is
        // passed over before any allele is read, as `||` (:266) and `!$b->{pos}`
        // (Sequence.pm:1173) pass it, whatever the alternate cell holds, unless
        // the shared-base move steps it.
        assert_eq!(
            alleles(row_side(
                &record(&["21", "5030000", "0", "G", "0.5", "9.0"]),
                source,
                CaseRule::AsWritten
            )),
            Some((b"-".to_vec(), b"G".to_vec(), "5030000".to_string(), 0))
        );
        assert_eq!(
            alleles(row_side(
                &record(&["21", "0", "A"]),
                source,
                CaseRule::AsWritten
            )),
            None
        );

        let sv = sv_fixture(dir.path());
        let mut plugin = CaddPlugin::new();
        plugin.init(&[format!("sv={}", sv.display())]).unwrap();
        let source = &plugin.sources[0];
        assert_eq!(
            alleles(row_side(
                &record(&["21", "5028999", "5047000", "", "name", "1.0", "0.1"]),
                source,
                CaseRule::AsWritten
            )),
            Some((b"-".to_vec(), b"".to_vec(), "5028999".to_string(), 0))
        );
        assert!(matches!(
            row_side(
                &record(&["21", "5028999", "5047000"]),
                source,
                CaseRule::AsWritten
            ),
            Err(PluginError::Run(_))
        ));
        assert!(matches!(
            row_side(
                &record(&["21", "5028999", "5047000", "", ""]),
                source,
                CaseRule::AsWritten
            ),
            Err(PluginError::Run(_))
        ));
    }

    /// A row's Pos cell is keyed by `get_matched_variant_alleles` as text
    /// (ensembl-variation Utils/Sequence.pm:1232) until a shared leading base
    /// moves it, when `$s++` (CADD.pm:283) steps the number `$s + length($ref)`
    /// (`:278`) made of it: for a substitution `2000` scores and `2000.0` and
    /// `2000.5` do not, nor does `2001x` at 2001; for the deletion `GA>G` at
    /// `3000.0` and `5000x` the moved start is the variant's, while `4000.5`
    /// moves to 4001.5 and matches nothing. Through made-up records built here, an
    /// exponent cell is its number once moved and nothing before, and a cell
    /// nothing numeric leads matches nothing either way.
    #[test]
    fn a_position_cell_is_matched_as_the_key_holds_it() {
        let dir = tempfile::tempdir().unwrap();
        let lines = [
            GRCH38_MARKER,
            SNV_HEADER,
            "21\t2000\tA\tG\t0.1\t1.0",
            "21\t2000.0\tA\tC\t0.2\t2.0",
            "21\t2000.5\tA\tT\t0.3\t3.0",
            "21\t2001x\tA\tG\t0.4\t4.0",
            "21\t3000.0\tGA\tG\t0.5\t5.0",
            "21\t4000.5\tGA\tG\t0.6\t6.0",
            "21\t5000x\tGA\tG\t0.7\t7.0",
        ];
        let path = write_tabix_fixture(
            dir.path(),
            "synthetic_cadd_cells.tsv.gz",
            &lines,
            IndexSpec::standard(),
        );
        let mut plugin = CaddPlugin::new();
        plugin.init(&[format!("snv={}", path.display())]).unwrap();
        let mut variants = vec![
            InputVariant::new("21".into(), 2000, 2000, b"A".to_vec(), b"G".to_vec()),
            InputVariant::new("21".into(), 2000, 2000, b"A".to_vec(), b"C".to_vec()),
            InputVariant::new("21".into(), 2000, 2000, b"A".to_vec(), b"T".to_vec()),
            InputVariant::new("21".into(), 2001, 2001, b"A".to_vec(), b"G".to_vec()),
            InputVariant::new("21".into(), 3001, 3001, b"A".to_vec(), b"-".to_vec()),
            InputVariant::new("21".into(), 4001, 4001, b"A".to_vec(), b"-".to_vec()),
            InputVariant::new("21".into(), 5001, 5001, b"A".to_vec(), b"-".to_vec()),
        ];
        plugin.run_batch(&mut variants).unwrap();
        assert_eq!(scores(&variants[0]), (Some("1.0"), Some("0.1")));
        assert_eq!(
            scores(&variants[1]),
            (None, None),
            "2000.0 is not the text 2000"
        );
        assert_eq!(
            scores(&variants[2]),
            (None, None),
            "2000.5 is not the text 2000"
        );
        assert_eq!(
            scores(&variants[3]),
            (None, None),
            "2001x is not the text 2001"
        );
        assert_eq!(
            scores(&variants[4]),
            (Some("5.0"), Some("0.5")),
            "3000.0 moved once is 3001"
        );
        assert_eq!(
            scores(&variants[5]),
            (None, None),
            "4000.5 moved once is 4001.5"
        );
        assert_eq!(
            scores(&variants[6]),
            (Some("7.0"), Some("0.7")),
            "5000x moved once is 5001"
        );

        let source = &plugin.sources[0];
        let record = |pos: &str, reference: &str, alt: &str| TabixRecord {
            chr: "21".to_string(),
            start: 6000,
            end: 6000,
            columns: ["21", pos, reference, alt, "0.9", "9.0"]
                .iter()
                .map(|cell| cell.to_string())
                .collect(),
            header: None,
        };
        let matches =
            |pos: &str, reference: &str, alt: &str, start: u64, var_ref: &[u8], var_alt: &[u8]| {
                let variant = VariantSide {
                    reference: Cow::Borrowed(var_ref),
                    alternate: Cow::Borrowed(var_alt),
                    position: start,
                    case: CaseRule::AsWritten,
                };
                let forms = variant.alleles().trimmed_forms(CaseRule::AsWritten);
                row_side(&record(pos, reference, alt), source, CaseRule::AsWritten)
                    .unwrap()
                    .is_some_and(|row| {
                        forms.matches(
                            &row.alleles().trimmed_forms(CaseRule::AsWritten),
                            CaseRule::AsWritten,
                        )
                    })
            };
        assert!(!matches("6e3", "A", "G", 6000, b"A", b"G"));
        assert!(matches("6e3", "GA", "G", 6001, b"A", b"-"));
        assert!(!matches("abc", "A", "G", 6000, b"A", b"G"));
        assert!(!matches("abc", "GA", "G", 6001, b"A", b"-"));
        assert!(!matches("0.0", "A", "G", 0, b"A", b"G"));
        // A Perl-false cell is the undefined value until the move steps it:
        // `undef++` and `"0"++` are 1 (Sequence.pm:1173 reads the moved value).
        assert!(matches("0", "GA", "G", 1, b"A", b"-"));
        assert!(matches("", "GA", "G", 1, b"A", b"-"));
        assert!(!matches("0", "A", "G", 0, b"A", b"G"));
        // A Perl-true cell that is no integer beside an undefined alternate fails
        // the variant, as the match dies on the allele after passing the cell.
        assert!(matches!(
            row_side(&record("2000.5", "A", ""), source, CaseRule::AsWritten),
            Err(PluginError::Run(_))
        ));

        // With a Perl-true `End` cell the move is `$s++` on the cell as a string
        // (:278 skips the arithmetic): `0100` steps to the text `0101`, which no
        // start equals, where the same cell without `End` steps to 101.
        let with_end = write_tabix_fixture(
            dir.path(),
            "synthetic_cadd_end.tsv.gz",
            &[
                GRCH38_MARKER,
                "#Chrom\tPos\tEnd\tRef\tAlt\tRawScore\tPHRED",
                "21\t0100\t101\tGA\tG\t0.1\t1.0",
                "21\t0200\t\tGA\tG\t0.2\t2.0",
            ],
            IndexSpec::standard(),
        );
        let mut plugin = CaddPlugin::new();
        plugin
            .init(&[format!("snv={}", with_end.display())])
            .unwrap();
        let source = &plugin.sources[0];
        // A made-up row in that table's shape, its `Pos` and `End` cells given.
        let end_row = |pos: &str, end: &str| TabixRecord {
            chr: "21".to_string(),
            start: 100,
            end: 100,
            columns: ["21", pos, end, "GA", "G", "0.1", "1.0"]
                .iter()
                .map(|cell| cell.to_string())
                .collect(),
            header: None,
        };
        let stepped = |pos: &str, end: &str| {
            row_side(&end_row(pos, end), source, CaseRule::AsWritten)
                .unwrap()
                .map(|row| (row.position.to_string(), row.moved))
        };
        assert_eq!(stepped("0100", "101"), Some(("0101".to_string(), 0)));
        assert_eq!(stepped("0100", ""), Some(("0100".to_string(), 1)));
        assert_eq!(stepped("0100", "0"), Some(("0100".to_string(), 1)));
        let mut variants = vec![
            InputVariant::new("21".into(), 101, 101, b"A".to_vec(), b"-".to_vec()),
            InputVariant::new("21".into(), 201, 201, b"A".to_vec(), b"-".to_vec()),
        ];
        plugin.run_batch(&mut variants).unwrap();
        assert_eq!(scores(&variants[0]), (None, None), "0101 is no start");
        assert_eq!(
            scores(&variants[1]),
            (Some("2.0"), Some("0.2")),
            "0200 moved once is 201"
        );
    }

    /// Writes beside `table`, under the name the plugin looks for, a store of
    /// one made-up row bound to the table's content and built on the columns
    /// of `built_on`. A store is accepted or refused for its columns and its
    /// content before a row of it is read.
    fn store_built_on(table: &Path, built_on: TabixConfig) -> PathBuf {
        use crate::binary_store::{BinaryStoreWriter, ContentIdentity};
        use std::sync::Arc;

        let store = PathBuf::from(table.to_str().unwrap().replace(".tsv.gz", ".tsv.vpd"));
        let header: Vec<String> = ["Chrom", "Start", "Name"].map(String::from).to_vec();
        let record = TabixRecord {
            chr: "21".to_string(),
            start: 100,
            end: 100,
            columns: ["21", "100", "B"].map(String::from).to_vec(),
            header: Some(Arc::new(header.clone())),
        };
        BinaryStoreWriter::new(store.clone())
            .write(
                &header,
                &[("21".to_string(), vec![record])],
                &TabixConfig {
                    file_path: table.to_path_buf(),
                    ..built_on
                },
                &ContentIdentity::of_file(table).unwrap(),
            )
            .unwrap();
        store
    }

    /// The columns the header of a file names, read as the plugin's own
    /// comparison reads them: a `Pos` cell is the base's own position
    /// (CADD.pm:209), a `Start` cell the base before a structural variant's
    /// first (CADD.pm:220), and an `End` column ends the row's span.
    #[test]
    fn header_named_config_places_a_pos_cell_one_based_and_a_start_cell_zero_based() {
        let columns = |position_col, position_is_start, end_col| {
            let config = header_named_config(
                PathBuf::from("scores.tsv.gz"),
                position_col,
                position_is_start,
                end_col,
            );
            (
                config.chr_col,
                config.start_col,
                config.end_col,
                config.ref_col,
                config.alt_col,
                config.zero_based,
            )
        };
        // A header with a `Pos` column.
        assert_eq!(columns(1, false, None), (0, 1, None, None, None, false));
        // A header with a `Start` and an `End` column.
        assert_eq!(columns(1, true, Some(2)), (0, 1, Some(2), None, None, true));
        // A header with a `Pos` and an `End` column.
        assert_eq!(
            columns(1, false, Some(2)),
            (0, 1, Some(2), None, None, false)
        );
    }

    /// A file with no index is read from the store beside it only where the
    /// store was built on the columns the header names. A store built on
    /// other columns would be read by them: its rows placed one base along,
    /// or as points where the file holds spans. It is refused instead, and
    /// the line that ends the load names the missing index, the store, its
    /// flags, the flags that rebuild it and the flag the conversion needs on
    /// a file with no index.
    #[test]
    fn a_file_with_no_index_is_read_from_a_store_built_on_the_columns_its_header_names() {
        let positions = TabixConfig::default();
        let zero_based_positions = TabixConfig {
            zero_based: true,
            ..TabixConfig::default()
        };
        let zero_based_spans = TabixConfig {
            end_col: Some(2),
            zero_based: true,
            ..TabixConfig::default()
        };
        let one_based_spans = TabixConfig {
            end_col: Some(2),
            ..TabixConfig::default()
        };
        let load = |table: &Path, built_on: &TabixConfig| {
            let store = store_built_on(table, built_on.clone());
            let outcome = CaddPlugin::open_source(table.to_str().unwrap(), None)
                .err()
                .map(|error| error.to_string());
            (store, outcome)
        };
        let refusal = |table: &Path, store: &Path, built_on: &str, named: &str| {
            format!(
                "plugin initialization failed: tabix index {table}.tbi not found; a plugin data \
                 file needs its .tbi index beside it, so index it with tabix without -C; the \
                 store {store} beside the table is refused: store columns `{built_on}` differ \
                 from the table's `{named}`; convert the table again with vep-cache-converter \
                 convert-plugin-data and the table's flags; a data file with no index is \
                 converted with --no-verify",
                table = table.display(),
                store = store.display()
            )
        };

        let dir = tempfile::tempdir().unwrap();
        let snv = cadd_fixture(dir.path(), "synthetic_snv.tsv.gz", "GRCh38");
        std::fs::remove_file(dir.path().join("synthetic_snv.tsv.gz.tbi")).unwrap();
        assert_eq!(load(&snv, &positions).1, None, "a store of positions");
        let (store, outcome) = load(&snv, &zero_based_positions);
        assert_eq!(
            outcome,
            Some(refusal(
                &snv,
                &store,
                "--chr-col 0 --start-col 1 --zero-based",
                "--chr-col 0 --start-col 1"
            ))
        );

        let sv = sv_fixture(dir.path());
        std::fs::remove_file(dir.path().join("synthetic_cadd_sv.tsv.gz.tbi")).unwrap();
        assert_eq!(
            load(&sv, &zero_based_spans).1,
            None,
            "a store of zero-based spans"
        );
        let named = "--chr-col 0 --start-col 1 --end-col 2 --zero-based";
        for (built_on, flags) in [
            (&positions, "--chr-col 0 --start-col 1"),
            (&one_based_spans, "--chr-col 0 --start-col 1 --end-col 2"),
        ] {
            let (store, outcome) = load(&sv, built_on);
            assert_eq!(outcome, Some(refusal(&sv, &store, flags, named)));
        }
    }
}
