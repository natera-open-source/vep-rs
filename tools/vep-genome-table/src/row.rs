// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! One row of the genome table: its nine columns, the grammar of each and the genome key.

/// The table's columns, in file order.
pub const COLUMNS: [&str; 9] = [
    "genome",
    "assembly",
    "transcripts",
    "set",
    "cache_version",
    "seed",
    "records",
    "outcome",
    "sha256",
];

/// The values of the `transcripts` column.
pub const TRANSCRIPTS: [&str; 3] = ["ensembl", "refseq", "merged"];

/// The sets named whole; every other set is a collection, in one of the forms
/// [`COLLECTION_FORMS`] names.
pub const NAMED_SETS: [&str; 6] = [
    "vertebrates",
    "metazoa",
    "plants",
    "fungi",
    "protists",
    "pangenome",
];

/// The forms of a collection set: `<N>` is an unsigned integer without a leading zero,
/// `<name>` lower-case letters and digits.
pub const COLLECTION_FORMS: [&str; 3] = [
    "bacteria_<N>_collection",
    "fungi_<name>_collection",
    "protists_<name>_collection",
];

const BACTERIA_PREFIX: &str = "bacteria_";
const FUNGI_PREFIX: &str = "fungi_";
const PROTISTS_PREFIX: &str = "protists_";
const COLLECTION_SUFFIX: &str = "_collection";
const KEY_SEPARATOR: &str = "__";

/// The header line of both table files and of a row file, without its newline.
pub fn header() -> String {
    COLUMNS.join("\t")
}

/// The genome key `<set>__<genome>__<assembly>__<transcripts>`. It is built from a row's
/// fields and never parsed back, since a genome or an assembly may itself hold `__`.
pub fn genome_key(set: &str, genome: &str, assembly: &str, transcripts: &str) -> String {
    [set, genome, assembly, transcripts].join(KEY_SEPARATOR)
}

/// The three kinds of outcome a row states.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OutcomeKind {
    Exact,
    Differs,
    Excluded,
}

/// The kind of `outcome`, or `None` when it is not `exact`, `differs:<name>` or
/// `excluded:<name>` with a name of lower-case letters, digits and hyphens.
pub fn outcome_kind(outcome: &str) -> Option<OutcomeKind> {
    if outcome == "exact" {
        return Some(OutcomeKind::Exact);
    }
    let (kind, name) = outcome.split_once(':')?;
    let kind = match kind {
        "differs" => OutcomeKind::Differs,
        "excluded" => OutcomeKind::Excluded,
        _ => return None,
    };
    let named = !name.is_empty()
        && name
            .bytes()
            .all(|b| b.is_ascii_lowercase() || b.is_ascii_digit() || b == b'-');
    named.then_some(kind)
}

/// The text between `prefix` and the collection suffix, or `None` when `set` lacks either.
fn collection_name<'a>(set: &'a str, prefix: &str) -> Option<&'a str> {
    set.strip_prefix(prefix)?.strip_suffix(COLLECTION_SUFFIX)
}

/// Lower-case letters and digits, one or more: the name of a fungal or a protist collection.
fn is_collection_name(text: &str) -> bool {
    !text.is_empty()
        && text
            .bytes()
            .all(|b| b.is_ascii_lowercase() || b.is_ascii_digit())
}

/// Whether `set` is one of the named sets or a collection in one of [`COLLECTION_FORMS`].
pub fn is_set(set: &str) -> bool {
    NAMED_SETS.contains(&set)
        || collection_name(set, BACTERIA_PREFIX).is_some_and(is_count)
        || collection_name(set, FUNGI_PREFIX).is_some_and(is_collection_name)
        || collection_name(set, PROTISTS_PREFIX).is_some_and(is_collection_name)
}

/// The sets a refusal names: the named sets, then the collection forms.
pub fn set_forms() -> String {
    let [bacteria, fungi, protists] = COLLECTION_FORMS;
    format!(
        "{}, {bacteria}, {fungi} or {protists}",
        NAMED_SETS.join(", ")
    )
}

/// Whether a row of `set` belongs to the compressed bacterial file.
pub fn is_bacterial_set(set: &str) -> bool {
    set.starts_with(BACTERIA_PREFIX)
}

/// Letters, digits, `.`, `_` and `-`, one or more: the grammar of a genome, an assembly and
/// a set's release.
pub fn is_name(text: &str) -> bool {
    !text.is_empty()
        && text
            .bytes()
            .all(|b| b.is_ascii_alphanumeric() || matches!(b, b'.' | b'_' | b'-'))
}

/// An unsigned decimal integer that fits 64 bits, written without a sign or a leading zero.
pub fn is_count(text: &str) -> bool {
    !text.is_empty()
        && (text == "0" || !text.starts_with('0'))
        && text.bytes().all(|b| b.is_ascii_digit())
        && text.parse::<u64>().is_ok()
}

/// An unsigned decimal integer above zero, in the form [`is_count`] takes.
pub fn is_positive(text: &str) -> bool {
    is_count(text) && text != "0"
}

fn is_sha256(text: &str) -> bool {
    text.len() == 64
        && text
            .bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}

/// One row of the table. Every field keeps the text the file holds, so a row is written
/// back byte for byte.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Row {
    pub genome: String,
    pub assembly: String,
    pub transcripts: String,
    pub set: String,
    pub cache_version: String,
    pub seed: String,
    pub records: String,
    pub outcome: String,
    pub sha256: String,
}

impl Row {
    /// Reads one tab-separated line, without its newline, and holds it to the field rules:
    /// nine fields, each in its column's grammar, and on an `exact` row a digest and
    /// records above zero. The error says what the line breaks.
    pub fn parse(line: &str) -> Result<Row, String> {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() != COLUMNS.len() {
            return Err(format!(
                "{} fields where the table has {} columns",
                fields.len(),
                COLUMNS.len()
            ));
        }
        let row = Row {
            genome: fields[0].to_owned(),
            assembly: fields[1].to_owned(),
            transcripts: fields[2].to_owned(),
            set: fields[3].to_owned(),
            cache_version: fields[4].to_owned(),
            seed: fields[5].to_owned(),
            records: fields[6].to_owned(),
            outcome: fields[7].to_owned(),
            sha256: fields[8].to_owned(),
        };
        row.validate()?;
        Ok(row)
    }

    fn validate(&self) -> Result<(), String> {
        if !is_name(&self.genome) {
            return Err(format!(
                "genome \"{}\" is not letters, digits, '.', '_' and '-'",
                self.genome
            ));
        }
        for suffix in ["_refseq", "_merged"] {
            if self.genome.ends_with(suffix) {
                return Err(format!(
                    "genome \"{}\" carries the suffix {suffix}, which the transcripts column states",
                    self.genome
                ));
            }
        }
        if !is_name(&self.assembly) {
            return Err(format!(
                "assembly \"{}\" is not letters, digits, '.', '_' and '-'",
                self.assembly
            ));
        }
        if !TRANSCRIPTS.contains(&self.transcripts.as_str()) {
            return Err(format!(
                "transcripts \"{}\" is not ensembl, refseq or merged",
                self.transcripts
            ));
        }
        if !is_set(&self.set) {
            return Err(format!("set \"{}\" is not {}", self.set, set_forms()));
        }
        if !is_positive(&self.cache_version) {
            return Err(format!(
                "cache_version \"{}\" is not an integer above zero",
                self.cache_version
            ));
        }
        if !is_count(&self.seed) {
            return Err(format!("seed \"{}\" is not an unsigned integer", self.seed));
        }
        if !is_count(&self.records) {
            return Err(format!(
                "records \"{}\" is not an unsigned integer",
                self.records
            ));
        }
        let Some(kind) = outcome_kind(&self.outcome) else {
            return Err(format!(
                "outcome \"{}\" is not exact, differs:<name> or excluded:<name>",
                self.outcome
            ));
        };
        if self.sha256 != "-" && !is_sha256(&self.sha256) {
            return Err(format!(
                "sha256 \"{}\" is not 64 lower-case hex digits or -",
                self.sha256
            ));
        }
        if kind == OutcomeKind::Exact {
            if self.sha256 == "-" {
                return Err("outcome exact needs a sha256, not -".to_owned());
            }
            if self.records == "0" {
                return Err("outcome exact needs records above zero".to_owned());
            }
        }
        Ok(())
    }

    /// The row's genome key.
    pub fn key(&self) -> String {
        genome_key(&self.set, &self.genome, &self.assembly, &self.transcripts)
    }

    /// Whether the row belongs to the compressed bacterial file.
    pub fn is_bacterial(&self) -> bool {
        is_bacterial_set(&self.set)
    }

    /// The fields in column order.
    pub fn fields(&self) -> [&str; 9] {
        [
            &self.genome,
            &self.assembly,
            &self.transcripts,
            &self.set,
            &self.cache_version,
            &self.seed,
            &self.records,
            &self.outcome,
            &self.sha256,
        ]
    }

    /// The field of `column`, or `None` when the table has no such column.
    pub fn field(&self, column: &str) -> Option<&str> {
        let index = COLUMNS.iter().position(|name| *name == column)?;
        Some(self.fields()[index])
    }

    /// The row as one tab-separated line, without its newline.
    pub fn line(&self) -> String {
        self.fields().join("\t")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn digest() -> String {
        "0123456789abcdef".repeat(4)
    }

    fn line(fields: [&str; 9]) -> String {
        fields.join("\t")
    }

    #[test]
    fn a_row_is_written_back_as_it_was_read() {
        let text = line([
            "homo_sapiens",
            "GRCh38",
            "merged",
            "vertebrates",
            "116",
            "7",
            "30",
            "exact",
            &digest(),
        ]);
        let row = Row::parse(&text).unwrap();
        assert_eq!(row.line(), text);
        assert_eq!(row.key(), "vertebrates__homo_sapiens__GRCh38__merged");
        assert_eq!(row.field("records"), Some("30"));
        assert_eq!(row.field("key"), None);
    }

    #[test]
    fn the_key_keeps_a_leading_underscore_of_the_genome() {
        assert_eq!(
            genome_key(
                "bacteria_12_collection",
                "_genus_species",
                "ASM1v1",
                "ensembl"
            ),
            "bacteria_12_collection___genus_species__ASM1v1__ensembl"
        );
    }

    #[test]
    fn the_set_grammar_takes_the_named_sets_and_numbered_collections() {
        for set in NAMED_SETS {
            assert!(is_set(set), "{set}");
        }
        assert!(is_set("bacteria_0_collection"));
        assert!(is_set("bacteria_128_collection"));
        for set in [
            "bacteria",
            "bacteria__collection",
            "bacteria_01_collection",
            "bacteria_x_collection",
            "Vertebrates",
            "",
        ] {
            assert!(!is_set(set), "{set}");
        }
        assert!(is_bacterial_set("bacteria_3_collection"));
        assert!(!is_bacterial_set("pangenome"));
    }

    #[test]
    fn the_set_grammar_takes_the_fungal_and_the_protist_collections() {
        for set in [
            "fungi_ascomycota1_collection",
            "fungi_a_collection",
            "protists_alveolata1_collection",
            "protists_7_collection",
        ] {
            assert!(is_set(set), "{set}");
            assert!(!is_bacterial_set(set), "{set}");
        }
        for set in [
            "fungi__collection",
            "fungi_collection",
            "fungi_Ascomycota1_collection",
            "fungi_asco_mycota_collection",
            "fungi_ascomycota1",
            "protists__collection",
            "protists_alveolata-1_collection",
            "protists_alveolata1_collections",
            "plants_one_collection",
            "archaea_1_collection",
        ] {
            assert!(!is_set(set), "{set}");
        }
    }

    #[test]
    fn a_refusal_names_every_form_of_set() {
        assert_eq!(
            set_forms(),
            "vertebrates, metazoa, plants, fungi, protists, pangenome, \
             bacteria_<N>_collection, fungi_<name>_collection or protists_<name>_collection"
        );
    }

    #[test]
    fn the_outcome_grammar_takes_three_forms() {
        assert_eq!(outcome_kind("exact"), Some(OutcomeKind::Exact));
        assert_eq!(
            outcome_kind("differs:unclassified"),
            Some(OutcomeKind::Differs)
        );
        assert_eq!(
            outcome_kind("excluded:no-transcripts-2"),
            Some(OutcomeKind::Excluded)
        );
        for outcome in [
            "",
            "exact:",
            "differs",
            "differs:",
            "differs:Upper",
            "differs:a_b",
            "other:x",
        ] {
            assert_eq!(outcome_kind(outcome), None, "{outcome}");
        }
    }

    #[test]
    fn each_field_rule_names_its_column() {
        let good = [
            "homo_sapiens",
            "GRCh38",
            "ensembl",
            "vertebrates",
            "116",
            "7",
            "30",
            "exact",
            "",
        ];
        let digest = digest();
        let cases: [(usize, &str, &str); 12] = [
            (0, "homo sapiens", "genome \"homo sapiens\""),
            (0, "homo_sapiens_refseq", "carries the suffix _refseq"),
            (0, "mus_musculus_merged", "carries the suffix _merged"),
            (1, "", "assembly \"\""),
            (2, "both", "transcripts \"both\""),
            (3, "bacteria", "set \"bacteria\""),
            (4, "0", "cache_version \"0\""),
            (5, "-1", "seed \"-1\""),
            (6, "030", "records \"030\""),
            (6, "0", "needs records above zero"),
            (7, "same", "outcome \"same\""),
            (8, "ABCDEF", "sha256 \"ABCDEF\""),
        ];
        for (index, value, expected) in cases {
            let mut fields = good;
            fields[8] = &digest;
            fields[index] = value;
            let error = Row::parse(&line(fields)).unwrap_err();
            assert!(error.contains(expected), "{value}: {error}");
        }
    }

    #[test]
    fn only_an_exact_row_must_carry_a_digest_and_records() {
        let row = |records: &str, outcome: &str| {
            Row::parse(&line([
                "homo_sapiens",
                "GRCh38",
                "ensembl",
                "pangenome",
                "107",
                "7",
                records,
                outcome,
                "-",
            ]))
        };
        assert!(row("0", "excluded:no-transcripts").is_ok());
        assert!(row("30", "differs:unclassified").is_ok());
        assert_eq!(
            row("30", "exact").unwrap_err(),
            "outcome exact needs a sha256, not -"
        );
    }

    #[test]
    fn a_line_of_another_width_is_refused_with_its_field_count() {
        let error = Row::parse("a\tb\tc").unwrap_err();
        assert_eq!(error, "3 fields where the table has 9 columns");
    }
}
