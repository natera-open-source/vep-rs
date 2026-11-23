// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Chromosome synonyms: the names a chromosome goes by in an input, a cache or
//! a FASTA (`21`, `chr21`, `NC_000021.9`, `CM000683.2`), read from a synonyms
//! file as Ensembl VEP reads it (`BaseVEP.pm` `chromosome_synonyms`,
//! release/116.3 lines 637-665) and applied to a source's own chromosome names
//! as `get_source_chr_name` applies them (lines 689-738).

use std::collections::HashMap;
use std::path::Path;

/// A chromosome synonym table.
///
/// Each later token of a line is a synonym of the line's first token, in both
/// directions, and of nothing else: the table is symmetric but not transitive,
/// so two synonyms of one reference name do not reach each other, as in
/// Ensembl VEP's table.
#[derive(Debug, Clone, Default)]
pub struct ChromosomeSynonyms {
    /// A name's synonyms, in the order the file first named them.
    synonyms: HashMap<String, Vec<String>>,
}

impl ChromosomeSynonyms {
    /// An empty table: every lookup falls back to the `chr` rule alone.
    pub fn new() -> Self {
        Self::default()
    }

    /// Reads a synonyms file as `chromosome_synonyms` reads it
    /// (`BaseVEP.pm:646-656`): each line is split on runs of whitespace, its
    /// first field is the reference name and every later field a synonym of it,
    /// recorded both ways. Perl's `split(/\s+/)` keeps a leading empty field, so
    /// a line opening with whitespace names the empty string as its reference
    /// and links the line's names to it, not to each other. A blank line carries
    /// no name and a line with one name names nothing. The error is the file's
    /// own, for the caller to word as Ensembl VEP does (`Could not read synonyms
    /// file <file>: <reason>`).
    pub fn load_from_file(path: &Path) -> std::io::Result<Self> {
        Ok(Self::parse(&std::fs::read_to_string(path)?))
    }

    /// The table a synonyms file's text describes ([`Self::load_from_file`]).
    pub fn parse(text: &str) -> Self {
        let mut table = Self::new();
        for line in text.lines() {
            let (reference, rest) = line.split_once(is_perl_space).unwrap_or((line, ""));
            for synonym in rest.split(is_perl_space).filter(|name| !name.is_empty()) {
                table.add_synonym(reference, synonym);
            }
        }
        table
    }

    /// Records `a` and `b` as synonyms of each other.
    pub fn add_synonym(&mut self, a: &str, b: &str) {
        self.link(a, b);
        self.link(b, a);
    }

    fn link(&mut self, from: &str, to: &str) {
        let list = self.synonyms.entry(from.to_string()).or_default();
        if !list.iter().any(|s| s == to) {
            list.push(to.to_string());
        }
    }

    /// The synonyms recorded for `chr`, in file order; empty for a name the
    /// table does not know.
    pub fn synonyms_of(&self, chr: &str) -> &[String] {
        self.synonyms.get(chr).map(Vec::as_slice).unwrap_or(&[])
    }

    /// The synonyms of `chr` a source holds, in the order the file first names
    /// them; none for a name the source holds itself. Two or more show a name
    /// that picks out no one sequence of the source: Ensembl VEP then takes the
    /// first in the order `keys` returns the name's synonyms (`BaseVEP.pm:711-716`),
    /// which changes between runs of one command.
    pub fn held_synonyms(&self, chr: &str, holds: impl Fn(&str) -> bool) -> Vec<&str> {
        if holds(chr) {
            return Vec::new();
        }
        let listed = self.synonyms_of(chr).iter();
        listed.map(String::as_str).filter(|s| holds(s)).collect()
    }

    /// [`Self::source_chr_name`] for a source whose names are distinct sequences,
    /// a cache's chromosomes: a name the table gives to two or more of them
    /// ([`Self::held_synonyms`]) stays as written, since the table does not say
    /// which of them it stands for, and so matches none.
    pub fn sequence_chr_name(&self, chr: &str, holds: impl Fn(&str) -> bool) -> String {
        let mut held = self.synonyms_of(chr).iter().filter(|s| holds(s));
        if !holds(chr) && held.nth(1).is_some() {
            return chr.to_string();
        }
        self.source_chr_name(chr, holds)
    }

    /// Whether the table records any synonym.
    pub fn is_empty(&self) -> bool {
        self.synonyms.is_empty()
    }

    /// Every name the table records a synonym for, each once, in no order.
    pub fn names(&self) -> impl Iterator<Item = &str> {
        self.synonyms.keys().map(String::as_str)
    }

    /// The name a source knows `chr` by, given which names the source holds
    /// (`BaseVEP.pm:697-734`): `chr` itself when the source holds it; else the
    /// first of its synonyms the source holds; else `chr` with a leading `chr`
    /// (in any case) removed, or with `chr` prepended, when the source holds
    /// that; else `chr` unchanged. A bare `M` whose source holds `MT` reads as
    /// `MT`, the one rename `Parser.pm:794-796` makes on the record, tried
    /// there as here after the synonyms and the `chr` rule.
    pub fn source_chr_name(&self, chr: &str, holds: impl Fn(&str) -> bool) -> String {
        if holds(chr) {
            return chr.to_string();
        }
        if let Some(synonym) = self.synonyms_of(chr).iter().find(|s| holds(s)) {
            return synonym.clone();
        }
        if let Some(stripped) = strip_chr_prefix(chr) {
            if holds(stripped) {
                return stripped.to_string();
            }
        } else {
            let prefixed = format!("chr{chr}");
            if holds(&prefixed) {
                return prefixed;
            }
            if chr == "M" && holds("MT") {
                return "MT".to_string();
            }
        }
        chr.to_string()
    }

    /// Every query name under which the reference sequence reader reaches a
    /// FASTA record called `name`. The reader tries, for a query it does not
    /// find, the query's `chr`-stripped or `chr`-prefixed spelling, that
    /// spelling's synonyms, then the query's own synonyms
    /// (`Variation/Utils/FastaSequence.pm:380-389`); so `name` is reached from
    /// its own `chr` spellings, from each of its synonyms, and from each
    /// synonym's `chr` spelling.
    pub fn aliases_of(&self, name: &str) -> Vec<String> {
        let mut out: Vec<String> = Vec::new();
        let mut push = |alias: String| {
            if alias != name && !out.contains(&alias) {
                out.push(alias);
            }
        };
        match strip_chr_prefix(name) {
            Some(stripped) => push(stripped.to_string()),
            None => {
                push(format!("chr{name}"));
                push(format!("CHR{name}"));
            }
        }
        for synonym in self.synonyms_of(name) {
            push(synonym.clone());
            // The query `chr<synonym>` reaches `name` through its stripped
            // spelling's synonyms; the query `<stripped synonym>` through
            // `chr<query>`'s synonyms.
            match strip_chr_prefix(synonym) {
                Some(stripped) => push(stripped.to_string()),
                None => push(format!("chr{synonym}")),
            }
        }
        out
    }
}

/// The characters Perl's `\s` matches in a file read as bytes: ASCII
/// whitespace plus the vertical tab, which `char::is_ascii_whitespace`
/// leaves out.
fn is_perl_space(c: char) -> bool {
    c.is_ascii_whitespace() || c == '\u{0B}'
}

/// `chr` without a leading `chr` in any case (`BaseVEP.pm:722-724`), or `None`
/// when it has none.
fn strip_chr_prefix(chr: &str) -> Option<&str> {
    let prefix = chr.get(..3)?;
    prefix.eq_ignore_ascii_case("chr").then(|| &chr[3..])
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    /// Lines of a release-116 GRCh38 cache's `chr_synonyms.txt`, which Ensembl
    /// VEP's dumper writes as one sorted pair per line
    /// (`Pipeline/DumpVEP/CreateDumpJobs.pm:200-209`).
    const CACHE_LINES: &str =
        "21\tNC_000021.9\n21\tchr21\n21\tCM000683.2\nMT\tchrM\nMT\tNC_012920.1\nX\tchrX\n";

    fn ensembl_cache(name: &str) -> bool {
        matches!(name, "21" | "22" | "X" | "Y" | "MT" | "KI270728.1")
    }

    /// `BaseVEP.pm:646-656`: the first token of a whitespace-split line is the
    /// reference, every later token its synonym, both directions, so a
    /// three-token line links each name to the reference and the reference to
    /// both; the two synonyms are not linked to each other, and a comment line
    /// is a line whose first token is `#`. `names` lists each linked name once.
    #[test]
    fn whitespace_split_first_token_reference() {
        let table = ChromosomeSynonyms::parse("21  NC_000021.9\tCM000683.2\n\n# chr21\nX\n");
        assert_eq!(table.synonyms_of("21"), ["NC_000021.9", "CM000683.2"]);
        assert_eq!(table.synonyms_of("NC_000021.9"), ["21"]);
        assert_eq!(table.synonyms_of("CM000683.2"), ["21"]);
        assert_eq!(table.synonyms_of("#"), ["chr21"]);
        assert!(table.synonyms_of("X").is_empty());
        assert!(table.synonyms_of("Y").is_empty());
        assert!(!table.is_empty());
        assert!(ChromosomeSynonyms::parse("").is_empty());
        let mut names: Vec<&str> = table.names().collect();
        names.sort_unstable();
        assert_eq!(names, ["#", "21", "CM000683.2", "NC_000021.9", "chr21"]);
    }

    /// `BaseVEP.pm:648-650`: `split(/\s+/, "  21\tNC_000021.9")` is `('', '21',
    /// 'NC_000021.9')`, so the reference is the empty string and both names are
    /// its synonyms, never each other's: `NC_000021.9` resolves to no cache name
    /// and `Parser.pm:542-546` skips its record. Trailing whitespace and a
    /// whitespace-only line yield no field, as Perl's `split` drops trailing
    /// empty fields; a tab or a vertical tab separates as a space does.
    #[test]
    fn leading_whitespace_makes_the_reference_empty() {
        let table = ChromosomeSynonyms::parse("  21\tNC_000021.9\n   \n21 chr21 \t\n\x0BX\tchrX\n");
        assert_eq!(table.synonyms_of("21"), ["", "chr21"]);
        assert_eq!(table.synonyms_of("NC_000021.9"), [""]);
        assert_eq!(table.synonyms_of(""), ["21", "NC_000021.9", "X", "chrX"]);
        assert_eq!(table.synonyms_of("chr21"), ["21"]);
        assert_eq!(table.synonyms_of("X"), [""]);
        assert_eq!(
            table.source_chr_name("NC_000021.9", ensembl_cache),
            "NC_000021.9"
        );
        assert_eq!(table.source_chr_name("chr21", ensembl_cache), "21");
        assert!(ChromosomeSynonyms::parse("   \n\t\n").is_empty());
    }

    /// `BaseVEP.pm:706-731`: `NC_000021.9` is not a cache chromosome, its synonym
    /// `21` is, so the source name is `21`; `CM000683.2` the same; a name the
    /// cache holds is returned as is; a name with no synonym the cache holds
    /// falls to the `chr` rule, and one that resolves nowhere is returned
    /// unchanged (`Parser.pm:542-546` then skips the variant).
    #[test]
    fn nc_accession_resolves_to_cache_name() {
        let table = ChromosomeSynonyms::parse(CACHE_LINES);
        let name = |chr: &str| table.source_chr_name(chr, ensembl_cache);
        assert_eq!(name("NC_000021.9"), "21");
        assert_eq!(name("CM000683.2"), "21");
        assert_eq!(name("21"), "21");
        assert_eq!(name("KI270728.1"), "KI270728.1");
        assert_eq!(name("chrM"), "MT");
        assert_eq!(name("NC_000022.11"), "NC_000022.11");
    }

    /// The `chr` rule without a synonyms file (`BaseVEP.pm:719-731`): `chr22` and
    /// `CHR22` strip to `22`, `21` gains `chr` only when the source holds
    /// `chr21`, `M` becomes `MT` when the source holds it (`Parser.pm:794-796`),
    /// and `chrM` resolves to nothing, as Ensembl VEP leaves it.
    #[test]
    fn chr_prefix_rule_runs_when_no_synonym_matches() {
        let table = ChromosomeSynonyms::new();
        assert_eq!(table.source_chr_name("chr22", ensembl_cache), "22");
        assert_eq!(table.source_chr_name("CHR22", ensembl_cache), "22");
        assert_eq!(table.source_chr_name("Chr22", ensembl_cache), "22");
        assert_eq!(table.source_chr_name("M", ensembl_cache), "MT");
        assert_eq!(table.source_chr_name("chrM", ensembl_cache), "chrM");
        assert_eq!(table.source_chr_name("21", |name| name == "chr21"), "chr21");
        assert_eq!(table.source_chr_name("chr", ensembl_cache), "chr");
    }

    /// A synonym the source holds wins over the `chr` rule, in file order
    /// (`BaseVEP.pm:711-716` runs before `:721-730`).
    #[test]
    fn a_held_synonym_wins_over_the_chr_rule() {
        let table = ChromosomeSynonyms::parse("chr21\tKI270728.1\n");
        assert_eq!(table.source_chr_name("chr21", ensembl_cache), "KI270728.1");
    }

    /// A file that gives one name to two sequences a source holds: the name's held
    /// synonyms are both, in the order the file first names them, and a name the
    /// source holds, or one the file gives to one held sequence, has no such pair.
    /// Ensembl VEP's dumper writes such a file where a synonym several sequences
    /// share links their other names (`BaseVEP.pm:599-618`).
    #[test]
    fn a_name_given_to_several_held_sequences_lists_each_in_file_order() {
        let table = ChromosomeSynonyms::parse("X\tshared\n21\tshared\n22\tshared\n21\tonly_21\n");
        assert_eq!(
            table.held_synonyms("shared", ensembl_cache),
            ["X", "21", "22"]
        );
        assert_eq!(table.held_synonyms("only_21", ensembl_cache), ["21"]);
        assert_eq!(table.held_synonyms("shared", |name| name == "22"), ["22"]);
        assert!(table.held_synonyms("unlisted", ensembl_cache).is_empty());
        assert!(table.held_synonyms("21", |name| name == "21").is_empty());
    }

    /// Among a cache's chromosomes a name given to two or more of them is read as
    /// none and stays as written; a name given to one is read as it, a held name
    /// as itself, and the `chr` rule still follows the synonyms.
    #[test]
    fn a_name_given_to_several_held_sequences_is_read_as_none_of_them() {
        let table = ChromosomeSynonyms::parse("X\tshared\n21\tshared\n21\tonly_21\n");
        assert_eq!(table.sequence_chr_name("shared", ensembl_cache), "shared");
        assert_eq!(table.source_chr_name("shared", ensembl_cache), "X");
        assert_eq!(table.sequence_chr_name("only_21", ensembl_cache), "21");
        assert_eq!(table.sequence_chr_name("shared", |name| name == "21"), "21");
        assert_eq!(table.sequence_chr_name("X", ensembl_cache), "X");
        assert_eq!(table.sequence_chr_name("chr22", ensembl_cache), "22");
        assert_eq!(table.sequence_chr_name("M", ensembl_cache), "MT");
    }

    /// The queries that reach a FASTA record (`FastaSequence.pm:380-389`): a
    /// record `NC_000021.9` is found from `21` (its synonym) and from `chr21`
    /// (stripped to `21`, whose synonyms hold it); a record `21` from `chr21`,
    /// from `NC_000021.9` and from `chrNC_000021.9` (stripped, then `21`'s
    /// synonym); a record `chr21` from `21` (`chr` prepended) and from nothing
    /// else, since `21`'s other synonyms name `21`, not `chr21`.
    #[test]
    fn aliases_follow_the_fasta_sequence_reader() {
        let table = ChromosomeSynonyms::parse(CACHE_LINES);
        assert_eq!(
            table.aliases_of("NC_000021.9"),
            ["chrNC_000021.9", "CHRNC_000021.9", "21", "chr21"]
        );
        assert_eq!(table.aliases_of("chr21"), ["21"]);
        assert_eq!(
            table.aliases_of("21"),
            [
                "chr21",
                "CHR21",
                "NC_000021.9",
                "chrNC_000021.9",
                "CM000683.2",
                "chrCM000683.2"
            ]
        );
        assert_eq!(ChromosomeSynonyms::new().aliases_of("X"), ["chrX", "CHRX"]);
    }

    /// `load_from_file` reads the file as `parse` reads its text, and a missing
    /// file is the reader's own error for the caller to word.
    #[test]
    fn cache_synonyms_loaded_from_json_cache() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("chr_synonyms.txt");
        std::fs::File::create(&path)
            .unwrap()
            .write_all(CACHE_LINES.as_bytes())
            .unwrap();
        let table = ChromosomeSynonyms::load_from_file(&path).unwrap();
        assert_eq!(table.source_chr_name("NC_000021.9", ensembl_cache), "21");
        assert_eq!(table.synonyms_of("MT"), ["chrM", "NC_012920.1"]);
        let missing = ChromosomeSynonyms::load_from_file(&dir.path().join("absent.txt"));
        assert_eq!(missing.unwrap_err().kind(), std::io::ErrorKind::NotFound);
    }

    /// Organelle sequence names in Ensembl release 116 and Ensembl Genomes
    /// release 63 caches.
    const ORGANELLES: [&str; 15] = [
        "MT",
        "Mt",
        "Mito",
        "MtDNA",
        "mitochondrion_genome",
        "Pt",
        "Chloro",
        "chloroplast",
        "OVOC_MITOCHONDRIAL",
        "SRAE_MITOCHONDRIAL",
        "SM_V7_MITO",
        "TMUE_MITO",
        "mito_A_fumigatus_Af293",
        "mito_C_glabrata_CBS138",
        "mito_C_parapsilosis_CDC317",
    ];

    /// A cache name of each shape beside one synonym of it. The accessions are
    /// the ones `docs/cache-species.md` lists for the name; `chrUn_GL456210v1`,
    /// the bacterial names (a chromosome named by a word, plasmids named with a
    /// dot or a hyphen) and every `SYN` synonym are made up in the shapes such
    /// caches carry.
    const NAME_SHAPES: [(&str, &str, &str); 15] = [
        ("a number", "4", "NC_003075.7"),
        ("a roman numeral", "I", "BK006935.2"),
        ("a chromosome arm", "2L", "AE014134.6"),
        ("a long digit run", "211000022278760", "SYN_SCAFFOLD_760"),
        (
            "an accession with a version",
            "GL456210.1",
            "chrUn_GL456210v1",
        ),
        (
            "a name with underscores",
            "Unmapped_Scaffold_8_D1580_D1567",
            "SYN_SCAFFOLD_8.1",
        ),
        (
            "a chromosome named by a word",
            "Chromosome",
            "NZ_SYN000000.1",
        ),
        ("a plasmid", "pSYN1", "NZ_SYN000001.1"),
        ("a plasmid with a dot", "pSYN1.2", "NZ_SYN000002.1"),
        ("a plasmid with a hyphen", "pSYN-b", "NZ_SYN000003.1"),
        ("a mitochondrion", "MT", "NC_005089.1"),
        ("a plant mitochondrion", "Mt", "NC_001284.2"),
        ("a chloroplast", "Pt", "NC_000932.1"),
        ("a yeast mitochondrion", "Mito", "AJ011856.1"),
        ("a fly mitochondrion", "mitochondrion_genome", "NC_024511.2"),
    ];

    /// A synonym line links its two names both ways (`BaseVEP.pm:652-655`), so
    /// for every name shape a record written with the synonym resolves to the
    /// cache's name, and a source named by the synonym, a reference or a custom
    /// file, is reached from the cache's name (`:706-716` over that source's own
    /// names). Without the line neither name reaches the other, and two
    /// synonyms of one name do not reach each other.
    #[test]
    fn every_name_shape_resolves_through_a_synonym_in_both_directions() {
        let organelles: Vec<(String, String)> = ORGANELLES
            .iter()
            .enumerate()
            .map(|(i, name)| (name.to_string(), format!("SYN_ORGANELLE_{i}.1")))
            .collect();
        let pairs = NAME_SHAPES
            .iter()
            .map(|&(shape, cache, synonym)| (shape, cache, synonym))
            .chain(
                organelles
                    .iter()
                    .map(|(cache, synonym)| ("an organelle", cache.as_str(), synonym.as_str())),
            );
        let none = ChromosomeSynonyms::new();
        for (shape, cache, synonym) in pairs {
            let table = ChromosomeSynonyms::parse(&format!("{cache}\t{synonym}\n"));
            let cache_holds = |name: &str| name == cache;
            let source_holds = |name: &str| name == synonym;
            assert_eq!(
                table.source_chr_name(synonym, cache_holds),
                cache,
                "{shape}"
            );
            assert_eq!(
                table.source_chr_name(cache, source_holds),
                synonym,
                "{shape}"
            );
            assert_eq!(
                none.source_chr_name(synonym, cache_holds),
                synonym,
                "{shape}"
            );
            assert_eq!(none.source_chr_name(cache, source_holds), cache, "{shape}");
        }

        let table = ChromosomeSynonyms::parse("4\tCP002687.1\tNC_003075.7\n");
        let cache_holds = |name: &str| name == "4";
        assert_eq!(table.source_chr_name("CP002687.1", cache_holds), "4");
        assert_eq!(table.source_chr_name("NC_003075.7", cache_holds), "4");
        assert_eq!(
            table.source_chr_name("CP002687.1", |name| name == "NC_003075.7"),
            "CP002687.1"
        );
    }

    /// A sequence whose own name begins with `chr`, in any case, is returned
    /// as it is when the source holds it (`BaseVEP.pm:706`), also where the
    /// source holds the stripped spelling too. On a miss the name is stripped
    /// (`:722-726`) and never prefixed, since `:728` is the `elsif` of the
    /// strip; a record written with `chr` before the whole name reaches it by
    /// the strip; and the stripped spelling reaches it only where the name
    /// opens with a lower-case `chr`, the one prefix `:728` adds. A synonym the
    /// source holds wins over the stripped and the prefixed spelling
    /// (`:711-716` runs first).
    #[test]
    fn a_name_that_itself_begins_with_chr_is_held_before_the_chr_rule() {
        let none = ChromosomeSynonyms::new();
        for name in [
            "chromosome",
            "Chromosome",
            "CHROMOSOME",
            "chr_1",
            "Chr1",
            "CHR_I",
            "chrUn_scaffold_7",
        ] {
            let stripped = &name[3..];
            let prefixed = format!("chr{name}");
            assert_eq!(none.source_chr_name(name, |held| held == name), name);
            assert_eq!(
                none.source_chr_name(name, |held| held == name || held == stripped),
                name
            );
            assert_eq!(
                none.source_chr_name(name, |held| held == stripped),
                stripped
            );
            assert_eq!(none.source_chr_name(name, |held| held == prefixed), name);
            assert_eq!(none.source_chr_name(&prefixed, |held| held == name), name);
            let reached = if name.starts_with("chr") {
                name
            } else {
                stripped
            };
            assert_eq!(
                none.source_chr_name(stripped, |held| held == name),
                reached,
                "{name}"
            );
            assert_eq!(
                none.source_chr_name(stripped, |held| held == name || held == stripped),
                stripped,
                "{name}"
            );
        }

        let table = ChromosomeSynonyms::parse("Chromosome\tNZ_SYN000000.1\n");
        assert_eq!(
            table.source_chr_name("Chromosome", |held| {
                held == "NZ_SYN000000.1" || held == "omosome"
            }),
            "NZ_SYN000000.1"
        );
        let table = ChromosomeSynonyms::parse("pSYN1\tNZ_SYN000001.1\n");
        assert_eq!(
            table.source_chr_name("pSYN1", |held| {
                held == "NZ_SYN000001.1" || held == "chrpSYN1"
            }),
            "NZ_SYN000001.1"
        );
    }

    /// An organelle is reached under its own name (`BaseVEP.pm:706`), under that
    /// name behind `chr` (`:722-726`) and under a name a synonym line links to
    /// it (`:711-716`), and under no other: without a line, `MT`, `M`, `chrM`
    /// and `chrMT` each stay as written, and a line linking one of them reaches
    /// the organelle from that spelling alone, since the `chr` rule tests the
    /// source's own names and never a synonym of the stripped spelling. The
    /// lookup reads a bare `M` as `MT` only where the source holds `MT` itself;
    /// a record's `M` is renamed before the lookup also where `MT` is a synonym
    /// of a held name (`Parser.pm:794-796` over the set `:782` widens). Where
    /// the source holds `MT`, `chrMT` strips to it, `chrM` strips to an `M` the
    /// source lacks, and a name differing in case is another name.
    #[test]
    fn an_organelle_is_reached_by_its_own_name_or_a_synonym_and_by_nothing_else() {
        const SPELLINGS: [&str; 4] = ["MT", "M", "chrM", "chrMT"];
        let none = ChromosomeSynonyms::new();
        for organelle in ORGANELLES {
            let holds = |name: &str| name == organelle || name == "1";
            assert_eq!(none.source_chr_name(organelle, holds), organelle);
            for prefix in ["chr", "Chr", "CHR"] {
                assert_eq!(
                    none.source_chr_name(&format!("{prefix}{organelle}"), holds),
                    organelle
                );
            }
            if organelle == "MT" {
                continue;
            }
            for spelling in SPELLINGS {
                assert_eq!(
                    none.source_chr_name(spelling, holds),
                    spelling,
                    "{spelling} on {organelle}"
                );
                let table = ChromosomeSynonyms::parse(&format!("{organelle}\t{spelling}\n"));
                assert_eq!(
                    table.source_chr_name(spelling, holds),
                    organelle,
                    "{spelling} linked to {organelle}"
                );
                for other in SPELLINGS.into_iter().filter(|other| *other != spelling) {
                    assert_eq!(
                        table.source_chr_name(other, holds),
                        other,
                        "{other} where {spelling} is linked to {organelle}"
                    );
                }
            }
        }

        let holds = |name: &str| name == "MT" || name == "1";
        assert_eq!(none.source_chr_name("M", holds), "MT");
        assert_eq!(none.source_chr_name("chrMT", holds), "MT");
        assert_eq!(none.source_chr_name("chrM", holds), "chrM");
        for other_case in ["Mt", "mt", "m", "chrm"] {
            assert_eq!(none.source_chr_name(other_case, holds), other_case);
        }
    }

    /// The names listed for a reference record (`FastaSequence.pm:380-389`),
    /// for a plasmid and for each organelle name: without a synonym, the
    /// record's name behind `chr` and behind `CHR`, two spellings of the prefix
    /// `:381-382` strips from a query; with one, a record named as the cache
    /// names it also lists the synonym and the synonym behind `chr`, or without
    /// its own `chr` where it opens with one (`:386`), and a record named by
    /// the synonym lists the cache's name and that name behind `chr`.
    #[test]
    fn the_reference_reader_aliases_cover_each_shape() {
        let none = ChromosomeSynonyms::new();
        assert_eq!(none.aliases_of("pSYN1.2"), ["chrpSYN1.2", "CHRpSYN1.2"]);
        for organelle in ORGANELLES {
            assert_eq!(
                none.aliases_of(organelle),
                [format!("chr{organelle}"), format!("CHR{organelle}")]
            );
        }

        let table = ChromosomeSynonyms::parse(
            "pSYN1.2\tNZ_SYN000002.1\nMito\tAJ011856.1\nMT\tchrM\nMT\tNC_005089.1\n",
        );
        assert_eq!(
            table.aliases_of("pSYN1.2"),
            [
                "chrpSYN1.2",
                "CHRpSYN1.2",
                "NZ_SYN000002.1",
                "chrNZ_SYN000002.1"
            ]
        );
        assert_eq!(
            table.aliases_of("NZ_SYN000002.1"),
            [
                "chrNZ_SYN000002.1",
                "CHRNZ_SYN000002.1",
                "pSYN1.2",
                "chrpSYN1.2"
            ]
        );
        assert_eq!(
            table.aliases_of("Mito"),
            ["chrMito", "CHRMito", "AJ011856.1", "chrAJ011856.1"]
        );
        assert_eq!(
            table.aliases_of("AJ011856.1"),
            ["chrAJ011856.1", "CHRAJ011856.1", "Mito", "chrMito"]
        );
        assert_eq!(
            table.aliases_of("MT"),
            [
                "chrMT",
                "CHRMT",
                "chrM",
                "M",
                "NC_005089.1",
                "chrNC_005089.1"
            ]
        );
        assert_eq!(table.aliases_of("Pt"), ["chrPt", "CHRPt"]);
    }
}
