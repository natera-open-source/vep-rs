// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! What a build leaves out of the cache it writes: the names `info.json`
//! lists under `left_out`, which a run reads to say what its command asks for
//! and the cache does not carry, and the lines the build ends on, one for
//! each, in the form of an installation's closing lines.

use vep_cache_builder::count::{counted, one_or};
use vep_cache_builder::gff3::Reading;
use vep_core::cache_info::left_out;

/// The output columns a cache built without the core table dumps holds no
/// value for, in the order `docs/cache-setup.md` lists them.
pub const CORE_TABLE_COLUMNS: [&str; 13] = [
    "SYMBOL_SOURCE",
    "HGNC_ID",
    "MANE_SELECT",
    "MANE_PLUS_CLINICAL",
    "APPRIS",
    "miRNA",
    "RefSeq",
    "SWISSPROT",
    "TREMBL",
    "UNIPARC",
    "UNIPROT_ISOFORM",
    "DOMAINS",
    "GENE_PHENO",
];

/// The facts of a build that settle what its cache leaves out.
pub struct Build<'a> {
    pub species: &'a str,
    pub assembly: &'a str,
    /// Whether the build read the core table dumps.
    pub core_dumps: bool,
    /// How the build read its transcript model.
    pub reading: Reading,
    /// The transcripts written with a translation that names no protein
    /// (`BuildStats::unnamed_translations`).
    pub unnamed_translations: usize,
    /// Whether the build wrote variation files.
    pub variation: bool,
    /// Whether the variation files hold the 1000 Genomes `AF` column.
    pub panel_frequency: bool,
    /// Whether the regulatory tree holds regulatory features.
    pub regulatory_features: bool,
    /// Whether the regulatory tree holds motif features.
    pub motif_features: bool,
    /// Whether the build read the tables the PubMed identifiers of the known
    /// variants come from (`VariantCellFacts::pubmed_tables_read`); a table
    /// that holds no row for the genome counts as read.
    pub pubmed: bool,
    /// Whether the build read the tables the variant synonyms come from
    /// (`VariantCellFacts::synonym_tables_read`).
    pub var_synonyms: bool,
    /// Whether the build read the failed-variants GVF or joined the
    /// `failed_variation` table to its known variants.
    pub failed_variants: bool,
    /// Whether every regulatory feature took its activity from an activity
    /// file (`RegulatoryStats::activity_on_every_feature`).
    pub cell_type_activity: bool,
    /// Whether every matrix a motif names carries its counts and every motif
    /// its sequence (`RegulatoryStats::counts_on_every_matrix` and
    /// `sequence_on_every_motif`).
    pub motif_scores: bool,
}

/// What a build leaves out of its cache.
#[derive(Debug, Default, PartialEq, Eq)]
pub struct LeftOut {
    /// The names `info.json` lists under `left_out`.
    pub names: Vec<&'static str>,
    /// The lines every build ends on, each a warning of its log.
    pub lines: Vec<String>,
    /// The lines on the cells of the variation files and of the regulatory
    /// tree that hold no value, which a build named one input at a time
    /// ends on after `lines`; an installation names those cells in its own
    /// closing lines, each with the reason its site gave.
    pub cell_lines: Vec<String>,
}

/// The route to a cell the build reads from no file.
const CONVERSION_CARRIES: &str =
    "a conversion of Ensembl VEP's cache carries them (docs/cache-conversion.md)";

/// `names` as a closing line joins the parts it has none of: `A`, `A or B`,
/// `A, B or C`.
pub(crate) fn any_of(names: &[&str]) -> String {
    match names.split_last() {
        None => String::new(),
        Some((only, [])) => only.to_string(),
        Some((last, rest)) => format!("{} or {last}", rest.join(", ")),
    }
}

/// What `build` leaves out. Without the core table dumps: the columns of
/// [`CORE_TABLE_COLUMNS`] and the mature miRNA arms `mature_miRNA_variant` is
/// read from. With the model read as a file source
/// ([`Reading::GffSource`]): a source name for `SOURCE`. With the 1000
/// Genomes columns: the panel-wide frequency as the release states it, in
/// whose place `AF` is the count recovered from the continental frequencies.
/// A translation that names no protein draws a line and no name, since no
/// flag of a run asks for what it lacks alone.
///
/// With variation files: each of the PubMed identifiers, the variant
/// synonyms and the failed-variant flags whose source the build did not read
/// (`Build::pubmed`, `Build::var_synonyms`, `Build::failed_variants`), and
/// for `homo_sapiens`, the one genome whose Ensembl VEP cache holds them,
/// the gnomAD frequencies and, without the 1000 Genomes columns, the 1000
/// Genomes frequencies; the build's own warning names a 1000 Genomes file
/// its inputs lack. With a regulatory tree: the cell-type activity of its
/// regulatory features, where a feature was written without the activity an
/// activity file gives it, and the matrix counts of its motif features,
/// where a matrix was written without its counts or a motif without its
/// sequence. A cell whose source the build read draws neither a name nor a
/// place in a line.
pub fn of(build: &Build<'_>) -> LeftOut {
    let cache = format!("the cache of {} {}", build.species, build.assembly);
    let human = build.species == crate::download::HUMAN_SPECIES;
    let mut left = LeftOut::default();
    if !build.core_dumps {
        let (last, first) = CORE_TABLE_COLUMNS
            .split_last()
            .expect("the list holds a column");
        left.names.push(left_out::CORE_TABLES);
        left.lines.push(format!(
            "No {} or {last} value in {cache}, and no mature_miRNA_variant from it: the build \
             read no core table dumps; --core-dumps reads the release's (docs/cache-setup.md)",
            first.join(", ")
        ));
    }
    if build.reading == Reading::GffSource {
        left.names.push(left_out::SOURCE);
        left.lines.push(format!(
            "No SOURCE in {cache}: its transcripts, read from a GTF or GFF3 with no protein \
             FASTA and no core table dumps, name no source, so a run writes no SOURCE field for \
             them (docs/cache-from-gtf-or-gff3.md)"
        ));
    }
    if build.unnamed_translations > 0 {
        left.lines.push(format!(
            "No ENSP or HGVSp for {} of {cache}: the CDS line that opens {} coding sequence \
             names no protein_id (docs/cache-from-gtf-or-gff3.md)",
            counted(build.unnamed_translations, "transcript"),
            one_or(build.unnamed_translations, "its", "each one's")
        ));
    }
    if build.panel_frequency {
        left.names.push(left_out::PANEL_FREQUENCY);
        left.lines.push(format!(
            "The 1000 Genomes AF in {cache} is each allele's count over the panel, recovered \
             from the continental frequencies, and can differ from the AF of Ensembl VEP's own \
             cache, which a conversion carries \
             (docs/cache-variation-data.md#the-1000-genomes-frequencies)"
        ));
    }
    if build.variation {
        if human && !build.panel_frequency {
            left.names.push(left_out::THOUSAND_GENOMES_FREQUENCIES);
        }
        let mut cells = Vec::new();
        if human {
            left.names.push(left_out::GNOMAD_FREQUENCIES);
            cells.push("gnomAD frequencies");
        }
        for (name, read, cell) in [
            (left_out::PUBMED, build.pubmed, "PubMed identifiers"),
            (
                left_out::VAR_SYNONYMS,
                build.var_synonyms,
                "variant synonyms",
            ),
            (
                left_out::FAILED_VARIANTS,
                build.failed_variants,
                "failed-variant flags",
            ),
        ] {
            if !read {
                left.names.push(name);
                cells.push(cell);
            }
        }
        if !cells.is_empty() {
            left.cell_lines.push(format!(
                "No {} in {cache}: the build reads them from no file; {CONVERSION_CARRIES}",
                any_of(&cells)
            ));
        }
    }
    let mut cells = Vec::new();
    if build.regulatory_features && !build.cell_type_activity {
        left.names.push(left_out::CELL_TYPE_ACTIVITY);
        cells.push("cell-type activity of regulatory features");
    }
    if build.motif_features && !build.motif_scores {
        left.names.push(left_out::MOTIF_SCORES);
        cells.push("HIGH_INF_POS or MOTIF_SCORE_CHANGE of motif features");
    }
    if !cells.is_empty() {
        left.cell_lines.push(format!(
            "No {} in {cache}: the build reads them from no file; {CONVERSION_CARRIES}",
            cells.join(", ")
        ));
    }
    left
}

#[cfg(test)]
mod tests {
    use super::*;

    fn build() -> Build<'static> {
        Build {
            species: "homo_sapiens",
            assembly: "GRCh38",
            core_dumps: true,
            reading: Reading::EnsemblCache,
            unnamed_translations: 0,
            variation: false,
            panel_frequency: false,
            regulatory_features: false,
            motif_features: false,
            pubmed: false,
            var_synonyms: false,
            failed_variants: false,
            cell_type_activity: false,
            motif_scores: false,
        }
    }

    /// A build of the cache's content with the core table dumps and neither
    /// variation files nor a regulatory tree leaves nothing out: no name and
    /// no line.
    #[test]
    fn a_build_with_the_core_tables_leaves_nothing_out() {
        assert_eq!(of(&build()), LeftOut::default());
    }

    /// Each fact draws its own line, in the order core tables, source,
    /// unnamed translations, panel frequency, and its own name; the unnamed
    /// translations draw no name. The lines are held whole.
    #[test]
    fn each_part_left_out_draws_its_name_and_its_line() {
        let left = of(&Build {
            core_dumps: false,
            reading: Reading::GffSource,
            unnamed_translations: 2,
            variation: true,
            panel_frequency: true,
            ..build()
        });
        assert_eq!(
            left.names,
            [
                "core_tables",
                "source",
                "panel_frequency",
                "gnomad_frequencies",
                "pubmed",
                "var_synonyms",
                "failed_variants"
            ]
        );
        assert_eq!(
            left.lines,
            [
                "No SYMBOL_SOURCE, HGNC_ID, MANE_SELECT, MANE_PLUS_CLINICAL, APPRIS, miRNA, \
                 RefSeq, SWISSPROT, TREMBL, UNIPARC, UNIPROT_ISOFORM, DOMAINS or GENE_PHENO value \
                 in the cache of homo_sapiens GRCh38, and no mature_miRNA_variant from it: the \
                 build read no core table dumps; --core-dumps reads the release's \
                 (docs/cache-setup.md)",
                "No SOURCE in the cache of homo_sapiens GRCh38: its transcripts, read from a GTF \
                 or GFF3 with no protein FASTA and no core table dumps, name no source, so a run \
                 writes no SOURCE field for them (docs/cache-from-gtf-or-gff3.md)",
                "No ENSP or HGVSp for 2 transcripts of the cache of homo_sapiens GRCh38: the CDS \
                 line that opens each one's coding sequence names no protein_id \
                 (docs/cache-from-gtf-or-gff3.md)",
                "The 1000 Genomes AF in the cache of homo_sapiens GRCh38 is each allele's count \
                 over the panel, recovered from the continental frequencies, and can differ from \
                 the AF of Ensembl VEP's own cache, which a conversion carries \
                 (docs/cache-variation-data.md#the-1000-genomes-frequencies)",
            ]
        );
    }

    /// One unnamed translation is counted in the singular, and a build that
    /// leaves out that alone writes no name.
    #[test]
    fn one_unnamed_translation_is_counted_in_the_singular() {
        let left = of(&Build {
            unnamed_translations: 1,
            ..build()
        });
        assert!(left.names.is_empty());
        assert_eq!(
            left.lines,
            [
                "No ENSP or HGVSp for 1 transcript of the cache of homo_sapiens GRCh38: the CDS \
              line that opens its coding sequence names no protein_id \
              (docs/cache-from-gtf-or-gff3.md)"
            ]
        );
    }

    /// The cells of the variation files: a build with them that read none of
    /// their sources names the PubMed identifiers, the variant synonyms and
    /// the failed-variant flags; `homo_sapiens` alone names the gnomAD
    /// frequencies and, without the 1000 Genomes columns, the 1000 Genomes
    /// frequencies. The cell line of a build named one input at a time names
    /// them in one line.
    #[test]
    fn variation_files_name_the_cells_no_file_gives() {
        let human = of(&Build {
            variation: true,
            ..build()
        });
        assert_eq!(
            human.names,
            [
                "thousand_genomes_frequencies",
                "gnomad_frequencies",
                "pubmed",
                "var_synonyms",
                "failed_variants"
            ]
        );
        assert!(human.lines.is_empty());
        assert_eq!(
            human.cell_lines,
            [
                "No gnomAD frequencies, PubMed identifiers, variant synonyms or failed-variant \
              flags in the cache of homo_sapiens GRCh38: the build reads them from no file; a \
              conversion of Ensembl VEP's cache carries them (docs/cache-conversion.md)"
            ]
        );
        let other = of(&Build {
            species: "danio_rerio",
            assembly: "GRCz11",
            variation: true,
            ..build()
        });
        assert_eq!(other.names, ["pubmed", "var_synonyms", "failed_variants"]);
        assert_eq!(
            other.cell_lines,
            [
                "No PubMed identifiers, variant synonyms or failed-variant flags in the cache of \
              danio_rerio GRCz11: the build reads them from no file; a conversion of Ensembl \
              VEP's cache carries them (docs/cache-conversion.md)"
            ]
        );
    }

    /// A cell whose source the build read draws no name and no place in the
    /// line: each of the PubMed identifiers, the variant synonyms and the
    /// failed-variant flags alone, and the three together, which leave
    /// `homo_sapiens` the gnomAD frequencies and another species no line.
    /// Without the fact the name stays.
    #[test]
    fn a_cell_of_the_variation_files_the_build_read_draws_no_name() {
        let line = |cells: &str| {
            format!(
                "No {cells} in the cache of homo_sapiens GRCh38: the build reads them from no \
                 file; a conversion of Ensembl VEP's cache carries them \
                 (docs/cache-conversion.md)"
            )
        };
        let variation = Build {
            variation: true,
            ..build()
        };
        let pubmed = of(&Build {
            pubmed: true,
            ..variation
        });
        assert_eq!(
            pubmed.names,
            [
                "thousand_genomes_frequencies",
                "gnomad_frequencies",
                "var_synonyms",
                "failed_variants"
            ]
        );
        assert_eq!(
            pubmed.cell_lines,
            [line(
                "gnomAD frequencies, variant synonyms or failed-variant flags"
            )]
        );
        let synonyms = of(&Build {
            var_synonyms: true,
            ..variation
        });
        assert_eq!(
            synonyms.names,
            [
                "thousand_genomes_frequencies",
                "gnomad_frequencies",
                "pubmed",
                "failed_variants"
            ]
        );
        assert_eq!(
            synonyms.cell_lines,
            [line(
                "gnomAD frequencies, PubMed identifiers or failed-variant flags"
            )]
        );
        let failed = of(&Build {
            failed_variants: true,
            ..variation
        });
        assert_eq!(
            failed.names,
            [
                "thousand_genomes_frequencies",
                "gnomad_frequencies",
                "pubmed",
                "var_synonyms"
            ]
        );
        assert_eq!(
            failed.cell_lines,
            [line(
                "gnomAD frequencies, PubMed identifiers or variant synonyms"
            )]
        );
        let every_cell = Build {
            pubmed: true,
            var_synonyms: true,
            failed_variants: true,
            ..variation
        };
        let human = of(&every_cell);
        assert_eq!(
            human.names,
            ["thousand_genomes_frequencies", "gnomad_frequencies"]
        );
        assert_eq!(human.cell_lines, [line("gnomAD frequencies")]);
        let other = of(&Build {
            species: "danio_rerio",
            assembly: "GRCz11",
            ..every_cell
        });
        assert_eq!(other, LeftOut::default());
    }

    /// A cell of the regulatory tree whose source the build read draws no
    /// name and no place in the line: the activity alone leaves the motif
    /// cells, the matrix counts alone the activity, and both leave nothing.
    /// A fact of a part the tree does not hold changes nothing.
    #[test]
    fn a_cell_of_the_regulatory_tree_the_build_read_draws_no_name() {
        let line = |cells: &str| {
            format!(
                "No {cells} in the cache of homo_sapiens GRCh38: the build reads them from no \
                 file; a conversion of Ensembl VEP's cache carries them \
                 (docs/cache-conversion.md)"
            )
        };
        let tree = Build {
            regulatory_features: true,
            motif_features: true,
            ..build()
        };
        let activity = of(&Build {
            cell_type_activity: true,
            ..tree
        });
        assert_eq!(activity.names, ["motif_scores"]);
        assert_eq!(
            activity.cell_lines,
            [line("HIGH_INF_POS or MOTIF_SCORE_CHANGE of motif features")]
        );
        let scores = of(&Build {
            motif_scores: true,
            ..tree
        });
        assert_eq!(scores.names, ["cell_type_activity"]);
        assert_eq!(
            scores.cell_lines,
            [line("cell-type activity of regulatory features")]
        );
        let both = of(&Build {
            cell_type_activity: true,
            motif_scores: true,
            ..tree
        });
        assert_eq!(both, LeftOut::default());
        let no_tree = of(&Build {
            cell_type_activity: true,
            motif_scores: true,
            pubmed: true,
            ..build()
        });
        assert_eq!(no_tree, LeftOut::default());
    }

    /// `names` joined as a closing line names what a cache has none of.
    #[test]
    fn names_are_joined_with_commas_and_a_last_or() {
        assert_eq!(any_of(&[]), "");
        assert_eq!(any_of(&["A"]), "A");
        assert_eq!(any_of(&["A", "B"]), "A or B");
        assert_eq!(any_of(&["A", "B", "C", "D"]), "A, B, C or D");
    }

    /// The cells of a regulatory tree whose build read no source of them:
    /// the activity where it holds regulatory features, the two motif cells
    /// where it holds motif features, each with its name, in one line.
    #[test]
    fn a_regulatory_tree_names_the_cells_no_file_gives() {
        let line = |cells: &str| {
            format!(
                "No {cells} in the cache of homo_sapiens GRCh38: the build reads them from no \
                 file; a conversion of Ensembl VEP's cache carries them \
                 (docs/cache-conversion.md)"
            )
        };
        let both = of(&Build {
            regulatory_features: true,
            motif_features: true,
            ..build()
        });
        assert_eq!(both.names, ["cell_type_activity", "motif_scores"]);
        assert_eq!(
            both.cell_lines,
            [line(
                "cell-type activity of regulatory features, HIGH_INF_POS or MOTIF_SCORE_CHANGE \
                 of motif features"
            )]
        );
        let features = of(&Build {
            regulatory_features: true,
            ..build()
        });
        assert_eq!(features.names, ["cell_type_activity"]);
        assert_eq!(
            features.cell_lines,
            [line("cell-type activity of regulatory features")]
        );
        let motifs = of(&Build {
            motif_features: true,
            ..build()
        });
        assert_eq!(motifs.names, ["motif_scores"]);
        assert_eq!(
            motifs.cell_lines,
            [line("HIGH_INF_POS or MOTIF_SCORE_CHANGE of motif features")]
        );
    }
}
