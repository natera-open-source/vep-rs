// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The variation files the binary writes with `--include-variations` are the
//! ones Ensembl VEP's cache carries: `variations/<chr>/all_vars.gz` with its
//! `.csi`, each record in the variation database's form, the 1000 Genomes
//! frequencies joined from the indexed frequency file, the somatic file's
//! records flagged, and, with `--variation-dumps`, the ClinVar allele and
//! somatic columns of the variation database dumps; `vep`'s own variation
//! source reads them back.

use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::process::Command;

use noodles::bgzf;
use noodles::core::Position;
use noodles::csi::binning_index::index::header::Format;
use noodles::csi::binning_index::index::reference_sequence::bin::Chunk;
use noodles::csi::binning_index::index::reference_sequence::index::BinnedIndex;
use noodles::csi::binning_index::index::Header;
use noodles::csi::binning_index::Indexer;
use serde_json::Value;
use vep_cli::variation_source::VariationSource;

fn fixture_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/grch38_chr21")
}

const VCF_HEADER: &str = "##fileformat=VCFv4.2\n\
    ##INFO=<ID=dbSNP_156,Number=0,Type=Flag,Description=\"dbSNP\">\n\
    ##INFO=<ID=COSMIC_102,Number=0,Type=Flag,Description=\"COSMIC\">\n\
    ##INFO=<ID=TSA,Number=1,Type=String,Description=\"Type of sequence alteration\">\n\
    ##INFO=<ID=E_Phenotype_or_Disease,Number=0,Type=Flag,Description=\"Phenotype_or_Disease\">\n\
    ##INFO=<ID=CLIN_pathogenic,Number=0,Type=Flag,Description=\"pathogenic\">\n\
    ##INFO=<ID=AFR,Number=A,Type=Float,Description=\"AFR\">\n\
    ##INFO=<ID=AMR,Number=A,Type=Float,Description=\"AMR\">\n\
    ##INFO=<ID=EAS,Number=A,Type=Float,Description=\"EAS\">\n\
    ##INFO=<ID=EUR,Number=A,Type=Float,Description=\"EUR\">\n\
    ##INFO=<ID=SAS,Number=A,Type=Float,Description=\"SAS\">\n\
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

/// A BGZF-compressed VCF named as the Ensembl FTP names it, with a `.csi` on
/// its CHROM and POS columns as Ensembl ships one, so the build reads the
/// genome-wide files through their index; `records` are sorted by chromosome
/// then position and name every chromosome in `chromosomes` order.
fn write_indexed_vcf(dir: &Path, name: &str, chromosomes: &[&str], records: &str) {
    write_indexed_vcf_with_header(dir, name, VCF_HEADER, chromosomes, records);
}

/// [`write_indexed_vcf`] with the header text given.
fn write_indexed_vcf_with_header(
    dir: &Path,
    name: &str,
    header: &str,
    chromosomes: &[&str],
    records: &str,
) {
    let path = dir.join(name);
    let mut writer = bgzf::io::Writer::new(std::fs::File::create(&path).unwrap());
    writer.write_all(header.as_bytes()).unwrap();
    let mut indexer = Indexer::<BinnedIndex>::new(14, 5);
    for line in records.lines() {
        let fields: Vec<&str> = line.split('\t').collect();
        let chromosome = chromosomes
            .iter()
            .position(|c| *c == fields[0])
            .expect("record chromosome is listed");
        let pos: usize = fields[1].parse().unwrap();
        let begin = writer.virtual_position();
        writeln!(writer, "{line}").unwrap();
        let end = writer.virtual_position();
        let start = Position::try_from(pos).unwrap();
        let last = Position::try_from(pos + fields[3].len() - 1).unwrap();
        indexer
            .add_record(
                Some((chromosome, start, last, true)),
                Chunk::new(begin, end),
            )
            .unwrap();
    }
    writer.try_finish().unwrap();
    let header = Header::builder()
        .set_format(Format::Vcf)
        .set_reference_sequence_name_index(0)
        .set_start_position_index(1)
        .set_end_position_index(None)
        .set_reference_sequence_names(chromosomes.iter().map(|c| c.as_bytes().into()).collect())
        .build();
    let index = indexer.set_header(header).build(chromosomes.len());
    noodles::csi::fs::write(format!("{}.csi", path.display()), &index).unwrap();
}

/// The variation database tables of one ClinVar span, 21:101-102: a germline
/// feature of rs1 (risk allele `-`, significance `pathogenic`, reference allele
/// `TG`) and a somatic feature with a classification.
fn write_dumps(dir: &Path) {
    let tables: &[(&str, &str)] = &[
        (
            "source",
            "1\tdbSNP\t156\t\\N\t\\N\t\\N\tgermline\tvariation\n32\tClinVar\t202509\t\\N\t\\N\t\\N\tmixed\tphenotype_feature\n",
        ),
        ("seq_region", "131537\t21\t4\n"),
        (
            "attrib_type",
            "480\tclinvar_clin_sig\tClinVar clinical significance\t\\N\n483\trisk_allele\tRisk allele\t\\N\n577\tsomatic_clin_sig\tClinVar somatic classification\t\\N\n583\tclinvar_ref_allele\tClinVar reference allele\t\\N\n",
        ),
        ("phenotype", "7\t\\N\t\\N\tBreast cancer\t\\N\n"),
        (
            "phenotype_feature",
            "1001\t7\t32\t\\N\tVariation\trs1\t1\t131537\t101\t102\t1\tGermline\n1002\t7\t32\t\\N\tVariation\trs1\t1\t131537\t101\t102\t1\tSomatic\n",
        ),
        (
            "phenotype_feature_attrib",
            "1001\t480\tpathogenic\n1001\t483\t-\n1001\t583\tTG\n1002\t577\tTier I - Strong\n",
        ),
    ];
    for (table, rows) in tables {
        std::fs::write(dir.join(format!("{table}.txt")), rows).unwrap();
    }
}

/// A `TG` deletion written as `ATG/A` at 100 is cached as `TG/-` at 101-102 with
/// the evidence flag read and the ClinVar columns joined on that span; an SNV on
/// another span carries no ClinVar column but its 1000 Genomes frequencies, read
/// through the frequency file's index (`AFR=0.0008`, one of 1,322 chromosomes,
/// is one of the panel's 5,008: `0.0001997`); the somatic file's record on
/// chromosome 21 is flagged and the one on chromosome 1 is not written.
/// `info.json` names the columns and the `tabix` layout, and `vep`'s variation
/// source reads every record back through the `.csi`, an insertion by its start.
#[test]
fn variation_files_carry_ensembl_form_and_read_back_through_the_runtime_source() {
    let fixture = fixture_dir();
    let inputs = tempfile::tempdir().unwrap();
    write_indexed_vcf(
        inputs.path(),
        "homo_sapiens-chr21.vcf.gz",
        &["21"],
        "21\t100\trs1\tATG\tA\t.\t.\tdbSNP_156;TSA=deletion;E_Phenotype_or_Disease;CLIN_pathogenic\n\
         21\t300\trs2\tA\tG\t.\t.\tdbSNP_156;TSA=SNV\n\
         21\t400\trs3\tG\tGG\t.\t.\tdbSNP_156;TSA=insertion\n",
    );
    write_indexed_vcf(
        inputs.path(),
        "1000GENOMES-phase_3.vcf.gz",
        &["1", "21"],
        "1\t300\trs9\tA\tG\t.\t.\tdbSNP_156;TSA=SNV;AFR=0.5;AMR=0.5;EAS=0.5;EUR=0.5;SAS=0.5\n\
         21\t300\trs2\tA\tG\t.\t.\tdbSNP_156;TSA=SNV;AFR=0.0008;AMR=0;EAS=0;EUR=0;SAS=0\n",
    );
    write_indexed_vcf(
        inputs.path(),
        "homo_sapiens_somatic.vcf.gz",
        &["1", "21"],
        "1\t100\tCOSV0000000\tT\t.\t.\t.\tCOSMIC_102;TSA=SNV;E_Phenotype_or_Disease\n\
         21\t200\tCOSV1\tC\t.\t.\t.\tCOSMIC_102;TSA=SNV;E_Phenotype_or_Disease\n",
    );
    write_dumps(inputs.path());
    let out = tempfile::tempdir().unwrap();

    let output = Command::new(env!("CARGO_BIN_EXE_vep-cache-builder"))
        .args([
            "--release",
            "116",
            "--assembly",
            "GRCh38",
            "--include-variations",
        ])
        .arg("--gff3")
        .arg(fixture.join("slice.gff3"))
        .arg("--fasta")
        .arg(fixture.join("peptides.fa"))
        .arg("--output-dir")
        .arg(out.path())
        .arg("--variation-vcf-dir")
        .arg(inputs.path())
        .arg("--variation-dumps")
        .arg(inputs.path())
        .output()
        .expect("vep-cache-builder runs");
    let log = String::from_utf8_lossy(&output.stderr);
    assert!(output.status.success(), "{log}");
    // The build ends on what the cache leaves out: the columns of the core
    // table dumps, which this build reads none of, the `AF` it recovers from
    // the continental frequencies, and the cells of the variation files it
    // reads from no file; `info.json` lists each by name.
    let closing: Vec<&str> = log
        .lines()
        .filter_map(|line| line.split_once(" WARN vep_cache_builder: "))
        .map(|(_, warning)| warning)
        .collect();
    assert_eq!(
        closing,
        [
            "No SYMBOL_SOURCE, HGNC_ID, MANE_SELECT, MANE_PLUS_CLINICAL, APPRIS, miRNA, RefSeq, \
             SWISSPROT, TREMBL, UNIPARC, UNIPROT_ISOFORM, DOMAINS or GENE_PHENO value in the \
             cache of homo_sapiens GRCh38, and no mature_miRNA_variant from it: the build read \
             no core table dumps; --core-dumps reads the release's (docs/cache-setup.md)",
            "The 1000 Genomes AF in the cache of homo_sapiens GRCh38 is each allele's count over \
             the panel, recovered from the continental frequencies, and can differ from the AF \
             of Ensembl VEP's own cache, which a conversion carries \
             (docs/cache-variation-data.md#the-1000-genomes-frequencies)",
            "No gnomAD frequencies, PubMed identifiers, variant synonyms or failed-variant flags \
             in the cache of homo_sapiens GRCh38: the build reads them from no file; a \
             conversion of Ensembl VEP's cache carries them (docs/cache-conversion.md)",
        ],
        "{log}"
    );
    assert_eq!(
        log.lines()
            .rev()
            .take(3)
            .filter(|line| line.contains(" WARN vep_cache_builder: "))
            .count(),
        3,
        "the three warnings are the build's last lines: {log}"
    );
    // One ClinVar span of each kind and one frequency record on the
    // chromosome, each counted in the singular.
    for line in [
        format!(
            "Read ClinVar features of 1 germline variant and 1 somatic span from {}",
            inputs.path().display()
        ),
        format!(
            "chr 21: 1 frequency record read from {} (0 skipped)",
            inputs.path().join("1000GENOMES-phase_3.vcf.gz").display()
        ),
    ] {
        assert_eq!(
            log.lines().filter(|logged| logged.ends_with(&line)).count(),
            1,
            "{line}: {log}"
        );
    }

    assert!(!out.path().join("variations/1").exists());
    let all_vars = out.path().join("variations/21/all_vars.gz");
    let mut text = String::new();
    bgzf::io::Reader::new(std::fs::File::open(&all_vars).unwrap())
        .read_to_string(&mut text)
        .unwrap();
    assert_eq!(
        text,
        "21\trs1\t.\t.\t101\t102\tTG/-\t.\tpathogenic\t1\t-:pathogenic\tTG\tBreast cancer(Tier I - Strong)\t.\t.\t.\t.\t.\t.\t.\t.\n\
         21\tCOSV1\t.\t1\t200\t.\tCOSMIC_MUTATION\t.\t.\t1\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\n\
         21\trs2\t.\t.\t300\t.\tA/G\t.\t.\t.\t.\t.\t.\t.\t.\tG:0.0001997\tG:0.0008\tG:0\tG:0\tG:0\tG:0\n\
         21\trs3\t.\t.\t402\t401\t-/G\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\n"
    );

    let info: Value =
        serde_json::from_str(&std::fs::read_to_string(out.path().join("info.json")).unwrap())
            .unwrap();
    assert_eq!(info["var_type"], "tabix");
    assert_eq!(
        info["left_out"],
        serde_json::json!([
            "core_tables",
            "panel_frequency",
            "gnomad_frequencies",
            "pubmed",
            "var_synonyms",
            "failed_variants"
        ])
    );
    let columns: Vec<&str> = info["variation_cols"]
        .as_array()
        .expect("variation_cols is written with the records")
        .iter()
        .map(|c| c.as_str().unwrap())
        .collect();
    assert_eq!(
        columns,
        [
            "chr",
            "variation_name",
            "failed",
            "somatic",
            "start",
            "end",
            "allele_string",
            "strand",
            "clin_sig",
            "phenotype_or_disease",
            "clin_sig_allele",
            "clin_sig_ref_allele",
            "clinical_impact",
            "pubmed",
            "var_synonyms",
            "AF",
            "AFR",
            "AMR",
            "EAS",
            "EUR",
            "SAS",
        ]
    );
    assert_eq!(
        columns.len(),
        text.lines().next().unwrap().split('\t').count()
    );

    let source = VariationSource::open_json_cache(out.path().to_str().unwrap(), false)
        .expect("the cache opens")
        .expect("the cache has a variation source");
    let names = |start: u64, end: u64| -> Vec<String> {
        source
            .records_overlapping("21", start, end)
            .unwrap()
            .iter()
            .map(|r| r.variation_name.clone())
            .collect()
    };
    assert_eq!(names(101, 102), ["rs1"]);
    assert_eq!(names(200, 200), ["COSV1"]);
    assert_eq!(names(402, 401), ["rs3"]);
    assert_eq!(names(1000, 1000), Vec::<String>::new());

    let snv = &source.records_overlapping("21", 300, 300).unwrap()[0];
    assert_eq!(snv.variation_name, "rs2");
    assert_eq!(snv.frequencies["AF"], "G:0.0001997");
    assert_eq!(snv.frequencies["AFR"], "G:0.0008");
    assert_eq!(snv.frequencies["SAS"], "G:0");
    assert_eq!(snv.frequencies.len(), 6);

    let deletion = &source.records_overlapping("21", 101, 102).unwrap()[0];
    assert_eq!((deletion.start, deletion.end), (101, 102));
    assert_eq!(deletion.allele_string, "TG/-");
    assert_eq!(deletion.clin_sig.as_deref(), Some("pathogenic"));
    assert_eq!(deletion.phenotype_or_disease, 1);
    let clinvar = deletion.clinvar.as_deref().expect("ClinVar columns read");
    assert_eq!(clinvar.clin_sig_allele.as_deref(), Some("-:pathogenic"));
    assert_eq!(clinvar.clin_sig_ref_allele.as_deref(), Some("TG"));
    assert_eq!(
        clinvar.clinical_impact.as_deref(),
        Some("Breast cancer(Tier I - Strong)")
    );

    let somatic = &source.records_overlapping("21", 200, 200).unwrap()[0];
    assert_eq!(somatic.somatic, 1);
    assert_eq!(somatic.allele_string, "COSMIC_MUTATION");
    assert!(somatic.frequencies.is_empty());
}

/// The germline ClinVar records of eleven made-up variants of chromosome 21,
/// one line a record: its `phenotype_feature_id`, the variant's name and
/// position, the record's allele and its significance as the tables spell it.
/// The lines stand as the test's `phenotype_feature` table holds its rows: by
/// falling identifier, the records of one variant among those of the others.
const CLINVAR_RECORDS: &[(u32, &str, u32, &str, &str)] = &[
    (2032, "rs9930000010", 1900, "T", "likely benign"),
    (2031, "rs9930000001", 1000, "C", "likely pathogenic"),
    (2030, "rs9930000008", 1700, "A", "likely pathogenic"),
    (2029, "rs9930000007", 1600, "A", "benign"),
    (2028, "rs9930000001", 1000, "C", "uncertain significance"),
    (2027, "rs9930000005", 1400, "A", "drug response"),
    (2026, "rs9930000004", 1300, "G", "likely pathogenic"),
    (2025, "rs9930000001", 1000, "A", "pathogenic"),
    (2024, "rs9930000002", 1100, "G", "uncertain significance"),
    (2023, "rs9930000001", 1000, "A", "likely pathogenic"),
    (2022, "rs9930000001", 1000, "A", "pathogenic"),
    (2021, "rs9930000002", 1100, "G", "benign"),
    (2020, "rs9930000003", 1200, "T", "likely pathogenic"),
    (2019, "rs9930000004", 1300, "C", "benign"),
    (2018, "rs9930000005", 1400, "A", "pathogenic"),
    (2017, "rs9930000006", 1500, "C", "likely benign"),
    (2016, "rs9930000007", 1600, "G", "uncertain significance"),
    (2015, "rs9930000008", 1700, "A", "uncertain significance"),
    (2014, "rs9930000009", 1800, "T", "pathogenic"),
    (2013, "rs9930000010", 1900, "T", "benign"),
    (2012, "rs9930000011", 2000, "G", "pathogenic"),
    (2011, "rs9930000011", 2000, "G", "pathogenic"),
    (2010, "rs9930000010", 1900, "T", "not provided"),
    (2009, "rs9930000009", 1800, "T", "likely pathogenic"),
    (2008, "rs9930000008", 1700, "A", "pathogenic"),
    (2007, "rs9930000007", 1600, "A", "likely benign"),
    (2006, "rs9930000006", 1500, "C", "benign"),
    (2005, "rs9930000005", 1400, "A", "risk factor"),
    (2004, "rs9930000004", 1300, "G", "pathogenic"),
    (2003, "rs9930000003", 1200, "T", "uncertain significance"),
    (2002, "rs9930000002", 1100, "G", "likely benign"),
    (2001, "rs9930000001", 1000, "A", "pathogenic"),
];

/// `clin_sig_allele` lists a variant's `allele:term` pairs in the order of its
/// ClinVar records' `phenotype_feature_id`, each different pair once, at the
/// place of the first record that carries it, whatever order the rows of the
/// `phenotype_feature` table stand in: the table of [`CLINVAR_RECORDS`] holds
/// them by falling identifier, the records of the variants interleaved.
/// Ensembl VEP joins the pairs from the keys of a hash
/// (`Database/Variation.pm:200-213`), which holds them in no order of its own.
/// `rs9930000001` has six records, by identifier four of `A` reading
/// `pathogenic`, `pathogenic`, `likely pathogenic` and `pathogenic` and two of
/// `C` reading `uncertain significance` and `likely pathogenic`, and its cell
/// holds four pairs. Nine variants have two or three records of different
/// pairs, the two records of `rs9930000011` carry one pair, written once, and
/// `rs9930000012` has no record and no cell.
#[test]
fn clinvar_pairs_stand_in_the_order_of_their_records_identifiers() {
    let fixture = fixture_dir();
    let inputs = tempfile::tempdir().unwrap();
    write_indexed_vcf(
        inputs.path(),
        "homo_sapiens-chr21.vcf.gz",
        &["21"],
        "21\t1000\trs9930000001\tG\tA,C\t.\t.\tdbSNP_156;TSA=SNV\n\
         21\t1100\trs9930000002\tA\tG\t.\t.\tdbSNP_156;TSA=SNV\n\
         21\t1200\trs9930000003\tC\tT\t.\t.\tdbSNP_156;TSA=SNV\n\
         21\t1300\trs9930000004\tT\tC,G\t.\t.\tdbSNP_156;TSA=SNV\n\
         21\t1400\trs9930000005\tG\tA\t.\t.\tdbSNP_156;TSA=SNV\n\
         21\t1500\trs9930000006\tA\tC\t.\t.\tdbSNP_156;TSA=SNV\n\
         21\t1600\trs9930000007\tC\tA,G\t.\t.\tdbSNP_156;TSA=SNV\n\
         21\t1700\trs9930000008\tT\tA\t.\t.\tdbSNP_156;TSA=SNV\n\
         21\t1800\trs9930000009\tG\tT\t.\t.\tdbSNP_156;TSA=SNV\n\
         21\t1900\trs9930000010\tA\tT\t.\t.\tdbSNP_156;TSA=SNV\n\
         21\t2000\trs9930000011\tC\tG\t.\t.\tdbSNP_156;TSA=SNV\n\
         21\t2100\trs9930000012\tT\tC\t.\t.\tdbSNP_156;TSA=SNV\n",
    );
    // `write_dumps` writes every table the build reads; the two that hold the
    // ClinVar records are replaced by the rows of `CLINVAR_RECORDS`, in its
    // order, each record a germline feature with its significance and allele.
    write_dumps(inputs.path());
    let features: String = CLINVAR_RECORDS
        .iter()
        .map(|(id, name, position, _, _)| {
            format!(
                "{id}\t7\t32\t\\N\tVariation\t{name}\t1\t131537\t{position}\t{position}\t1\t\
                 Germline\n"
            )
        })
        .collect();
    let attributes: String = CLINVAR_RECORDS
        .iter()
        .map(|(id, _, _, allele, term)| format!("{id}\t480\t{term}\n{id}\t483\t{allele}\n"))
        .collect();
    std::fs::write(inputs.path().join("phenotype_feature.txt"), features).unwrap();
    std::fs::write(
        inputs.path().join("phenotype_feature_attrib.txt"),
        attributes,
    )
    .unwrap();
    let out = tempfile::tempdir().unwrap();

    let output = Command::new(env!("CARGO_BIN_EXE_vep-cache-builder"))
        .args([
            "--release",
            "116",
            "--assembly",
            "GRCh38",
            "--include-variations",
        ])
        .arg("--gff3")
        .arg(fixture.join("slice.gff3"))
        .arg("--fasta")
        .arg(fixture.join("peptides.fa"))
        .arg("--output-dir")
        .arg(out.path())
        .arg("--variation-vcf-dir")
        .arg(inputs.path())
        .arg("--variation-dumps")
        .arg(inputs.path())
        .output()
        .expect("vep-cache-builder runs");
    let log = String::from_utf8_lossy(&output.stderr);
    assert!(output.status.success(), "{log}");

    let info: Value =
        serde_json::from_str(&std::fs::read_to_string(out.path().join("info.json")).unwrap())
            .unwrap();
    let column = |name: &str| -> usize {
        info["variation_cols"]
            .as_array()
            .expect("variation_cols is written with the records")
            .iter()
            .position(|c| c == name)
            .unwrap_or_else(|| panic!("info.json names no {name} column"))
    };
    let (name, cell) = (column("variation_name"), column("clin_sig_allele"));
    let mut text = String::new();
    bgzf::io::Reader::new(
        std::fs::File::open(out.path().join("variations/21/all_vars.gz")).unwrap(),
    )
    .read_to_string(&mut text)
    .unwrap();
    let cells: Vec<(&str, &str)> = text
        .lines()
        .map(|line| {
            let fields: Vec<&str> = line.split('\t').collect();
            (fields[name], fields[cell])
        })
        .collect();
    assert_eq!(
        cells,
        [
            (
                "rs9930000001",
                "A:pathogenic;A:likely_pathogenic;C:uncertain_significance;C:likely_pathogenic"
            ),
            (
                "rs9930000002",
                "G:likely_benign;G:benign;G:uncertain_significance"
            ),
            (
                "rs9930000003",
                "T:uncertain_significance;T:likely_pathogenic"
            ),
            ("rs9930000004", "G:pathogenic;C:benign;G:likely_pathogenic"),
            ("rs9930000005", "A:risk_factor;A:pathogenic;A:drug_response"),
            ("rs9930000006", "C:benign;C:likely_benign"),
            (
                "rs9930000007",
                "A:likely_benign;G:uncertain_significance;A:benign"
            ),
            (
                "rs9930000008",
                "A:pathogenic;A:uncertain_significance;A:likely_pathogenic"
            ),
            ("rs9930000009", "T:likely_pathogenic;T:pathogenic"),
            ("rs9930000010", "T:not_provided;T:benign;T:likely_benign"),
            ("rs9930000011", "G:pathogenic"),
            ("rs9930000012", "."),
        ]
    );
}

/// The germline ClinVar records of five made-up variants of chromosome 21, in
/// the form of [`CLINVAR_RECORDS`], by rising identifier. A significance that
/// holds more than one classification joins them with a semicolon and a
/// blank, as the `clinvar_clin_sig` attribute of the variation database does.
const CLINVAR_RECORDS_OF_SEVERAL_CLASSIFICATIONS: &[(u32, &str, u32, &str, &str)] = &[
    (
        3001,
        "rs9940000001",
        1000,
        "C",
        "uncertain significance; risk factor",
    ),
    (3002, "rs9940000001", 1000, "C", "likely pathogenic"),
    (3003, "rs9940000001", 1000, "C", "not provided"),
    (3004, "rs9940000001", 1000, "C", "uncertain significance"),
    (3005, "rs9940000001", 1000, "C", "likely pathogenic"),
    (3006, "rs9940000001", 1000, "C", "likely benign"),
    (3007, "rs9940000001", 1000, "C", "uncertain significance"),
    (3008, "rs9940000001", 1000, "C", "uncertain significance"),
    (3011, "rs9940000002", 1100, "T", "risk factor"),
    (3012, "rs9940000002", 1100, "T", "benign/likely benign"),
    (
        3013,
        "rs9940000002",
        1100,
        "T",
        "benign/likely benign; other",
    ),
    (3014, "rs9940000002", 1100, "T", "benign/likely benign"),
    (3015, "rs9940000002", 1100, "T", "likely benign"),
    (
        3016,
        "rs9940000002",
        1100,
        "T",
        "conflicting classifications of pathogenicity",
    ),
    (
        3017,
        "rs9940000002",
        1100,
        "T",
        "likely benign; risk factor",
    ),
    (3018, "rs9940000002", 1100, "T", "benign"),
    (
        3019,
        "rs9940000002",
        1100,
        "T",
        "likely benign; risk factor",
    ),
    (3020, "rs9940000002", 1100, "T", "uncertain significance"),
    (
        3021,
        "rs9940000003",
        1200,
        "A",
        "conflicting classifications of pathogenicity; association; risk factor",
    ),
    (3031, "rs9940000004", 1300, "G", "pathogenic; other"),
    (3032, "rs9940000004", 1300, "T", "benign; other"),
    (3033, "rs9940000004", 1300, "G", "pathogenic; drug response"),
    (
        3041,
        "rs9940000005",
        1400,
        "A",
        "pathogenic/pathogenic, low penetrance",
    ),
];

/// A ClinVar significance that holds more than one classification gives
/// `clin_sig_allele` one `allele:term` pair for each, at the place of the
/// first record that carries it, so that every piece of the `;`-joined cell is
/// a pair and no pair stands twice. Ensembl VEP's adaptor cuts such a value at
/// its first semicolon (`PhenotypeFeatureAdaptor.pm:1273-1276`), so its cell
/// holds a pair for the first classification alone.
///
/// `rs9940000001` has eight records and `rs9940000002` ten: the first gets
/// `risk factor` from one value alone, after a semicolon; the second gets
/// `other` so, while `risk factor` and `likely benign` of its later values
/// are pairs an earlier record wrote. The one value of `rs9940000003` holds
/// three classifications. The two alleles of `rs9940000004` each get the
/// `other` of their own record. A value without a semicolon is one
/// classification, a comma and a slash inside it included (`rs9940000005`).
#[test]
fn a_clinvar_significance_of_several_classifications_gives_a_pair_for_each() {
    let fixture = fixture_dir();
    let inputs = tempfile::tempdir().unwrap();
    write_indexed_vcf(
        inputs.path(),
        "homo_sapiens-chr21.vcf.gz",
        &["21"],
        "21\t1000\trs9940000001\tT\tC\t.\t.\tdbSNP_156;TSA=SNV\n\
         21\t1100\trs9940000002\tC\tT\t.\t.\tdbSNP_156;TSA=SNV\n\
         21\t1200\trs9940000003\tG\tA\t.\t.\tdbSNP_156;TSA=SNV\n\
         21\t1300\trs9940000004\tA\tG,T\t.\t.\tdbSNP_156;TSA=SNV\n\
         21\t1400\trs9940000005\tC\tA\t.\t.\tdbSNP_156;TSA=SNV\n",
    );
    // `write_dumps` writes every table the build reads; the two that hold the
    // ClinVar records take the rows of the list above.
    write_dumps(inputs.path());
    let features: String = CLINVAR_RECORDS_OF_SEVERAL_CLASSIFICATIONS
        .iter()
        .map(|(id, name, position, _, _)| {
            format!(
                "{id}\t7\t32\t\\N\tVariation\t{name}\t1\t131537\t{position}\t{position}\t1\t\
                 Germline\n"
            )
        })
        .collect();
    let attributes: String = CLINVAR_RECORDS_OF_SEVERAL_CLASSIFICATIONS
        .iter()
        .map(|(id, _, _, allele, term)| format!("{id}\t480\t{term}\n{id}\t483\t{allele}\n"))
        .collect();
    std::fs::write(inputs.path().join("phenotype_feature.txt"), features).unwrap();
    std::fs::write(
        inputs.path().join("phenotype_feature_attrib.txt"),
        attributes,
    )
    .unwrap();
    let out = tempfile::tempdir().unwrap();

    let output = Command::new(env!("CARGO_BIN_EXE_vep-cache-builder"))
        .args([
            "--release",
            "116",
            "--assembly",
            "GRCh38",
            "--include-variations",
        ])
        .arg("--gff3")
        .arg(fixture.join("slice.gff3"))
        .arg("--fasta")
        .arg(fixture.join("peptides.fa"))
        .arg("--output-dir")
        .arg(out.path())
        .arg("--variation-vcf-dir")
        .arg(inputs.path())
        .arg("--variation-dumps")
        .arg(inputs.path())
        .output()
        .expect("vep-cache-builder runs");
    let log = String::from_utf8_lossy(&output.stderr);
    assert!(output.status.success(), "{log}");

    let info: Value =
        serde_json::from_str(&std::fs::read_to_string(out.path().join("info.json")).unwrap())
            .unwrap();
    let column = |name: &str| -> usize {
        info["variation_cols"]
            .as_array()
            .expect("variation_cols is written with the records")
            .iter()
            .position(|c| c == name)
            .unwrap_or_else(|| panic!("info.json names no {name} column"))
    };
    let (name, cell) = (column("variation_name"), column("clin_sig_allele"));
    let mut text = String::new();
    bgzf::io::Reader::new(
        std::fs::File::open(out.path().join("variations/21/all_vars.gz")).unwrap(),
    )
    .read_to_string(&mut text)
    .unwrap();
    let cells: Vec<(&str, &str)> = text
        .lines()
        .map(|line| {
            let fields: Vec<&str> = line.split('\t').collect();
            (fields[name], fields[cell])
        })
        .collect();
    assert_eq!(
        cells,
        [
            (
                "rs9940000001",
                "C:uncertain_significance;C:risk_factor;C:likely_pathogenic;C:not_provided;\
                 C:likely_benign"
            ),
            (
                "rs9940000002",
                "T:risk_factor;T:benign/likely_benign;T:other;T:likely_benign;\
                 T:conflicting_classifications_of_pathogenicity;T:benign;T:uncertain_significance"
            ),
            (
                "rs9940000003",
                "A:conflicting_classifications_of_pathogenicity;A:association;A:risk_factor"
            ),
            (
                "rs9940000004",
                "G:pathogenic;G:other;T:benign;T:other;G:drug_response"
            ),
            ("rs9940000005", "A:pathogenic/pathogenic,_low_penetrance"),
        ]
    );
    // Every piece of every cell is one `allele:term` pair.
    for (variant, cell) in &cells {
        for piece in cell.split(';') {
            let (allele, term) = piece
                .split_once(':')
                .unwrap_or_else(|| panic!("{variant}: {piece} is no pair"));
            assert!(
                !allele.is_empty() && !term.is_empty() && !term.starts_with('_'),
                "{variant}: {piece}"
            );
        }
    }
}

/// A build named one input at a time ends on the warning an installation
/// ends on for a sequence it read no variation record on, naming the build,
/// between its own lines and the line of the cells. The one file of another
/// species holds a record on chromosome 22 alone and the model's transcripts
/// stand on chromosome 21, so the warning names 21. The build of
/// `variation_files_carry_ensembl_form_and_read_back_through_the_runtime_source`,
/// whose chromosome holds records, prints no such warning: its closing
/// warnings are held whole.
#[test]
fn a_build_names_the_sequence_its_variation_file_holds_no_record_on() {
    let fixture = fixture_dir();
    let inputs = tempfile::tempdir().unwrap();
    write_indexed_vcf(
        inputs.path(),
        "danio_rerio.vcf.gz",
        &["22"],
        "22\t100\trs7\tA\tG\t.\t.\tdbSNP_156;TSA=SNV\n",
    );
    let out = tempfile::tempdir().unwrap();
    let output = Command::new(env!("CARGO_BIN_EXE_vep-cache-builder"))
        .args(["--release", "116", "--assembly", "GRCz11"])
        .args(["--species", "danio_rerio", "--include-variations"])
        .arg("--gff3")
        .arg(fixture.join("slice.gff3"))
        .arg("--fasta")
        .arg(fixture.join("peptides.fa"))
        .arg("--output-dir")
        .arg(out.path())
        .arg("--variation-vcf-dir")
        .arg(inputs.path())
        .output()
        .expect("vep-cache-builder runs");
    let log = String::from_utf8_lossy(&output.stderr);
    assert!(output.status.success(), "{log}");
    assert!(!out.path().join("variations/21").exists());
    let closing: Vec<&str> = log
        .lines()
        .filter_map(|line| line.split_once(" WARN vep_cache_builder: "))
        .map(|(_, warning)| warning)
        .collect();
    assert_eq!(
        closing,
        [
            "No SYMBOL_SOURCE, HGNC_ID, MANE_SELECT, MANE_PLUS_CLINICAL, APPRIS, miRNA, RefSeq, \
             SWISSPROT, TREMBL, UNIPARC, UNIPROT_ISOFORM, DOMAINS or GENE_PHENO value in the \
             cache of danio_rerio GRCz11, and no mature_miRNA_variant from it: the build read no \
             core table dumps; --core-dumps reads the release's (docs/cache-setup.md)",
            "No known variants on sequence 21 in the cache of danio_rerio GRCz11: the build read \
             no variation record there, in a variation VCF or in the variation_feature table; a \
             conversion of Ensembl VEP's cache carries any known variant that cache holds there \
             (docs/cache-conversion.md)",
            "No PubMed identifiers, variant synonyms or failed-variant flags in the cache of \
             danio_rerio GRCz11: the build reads them from no file; a conversion of Ensembl \
             VEP's cache carries them (docs/cache-conversion.md)",
        ],
        "{log}"
    );
}

/// A species other than `homo_sapiens` is built from the one genome-wide
/// variation file Ensembl publishes for it (`<species>.vcf.gz`), read through
/// its index: with no `--chromosomes`, each sequence the file holds records on
/// gets an `all_vars.gz` of its own records, and the genome-wide files Ensembl
/// publishes for human alone are not warned about, where a human build that
/// lacks them is. A directory holding neither a per-chromosome file nor the
/// genome-wide one ends the build naming both names looked for.
#[test]
fn another_species_is_built_from_its_one_genome_wide_variation_file() {
    let fixture = fixture_dir();
    let build = |species: &str, assembly: &str, vcf_dir: &Path, out: &Path| {
        Command::new(env!("CARGO_BIN_EXE_vep-cache-builder"))
            .args(["--release", "116", "--assembly", assembly])
            .args(["--species", species, "--include-variations"])
            .arg("--gff3")
            .arg(fixture.join("slice.gff3"))
            .arg("--fasta")
            .arg(fixture.join("peptides.fa"))
            .arg("--output-dir")
            .arg(out)
            .arg("--variation-vcf-dir")
            .arg(vcf_dir)
            .output()
            .expect("vep-cache-builder runs")
    };
    let human_only_file_warnings = |stderr: &str| {
        stderr
            .lines()
            .filter(|line| line.contains("WARN"))
            .filter(|line| {
                ["1000GENOMES", "somatic", "_associated"]
                    .iter()
                    .any(|name| line.contains(name))
            })
            .count()
    };
    let names = |cache: &Path, chromosome: &str| -> Vec<String> {
        let mut text = String::new();
        let all_vars = cache.join(format!("variations/{chromosome}/all_vars.gz"));
        bgzf::io::Reader::new(std::fs::File::open(&all_vars).unwrap())
            .read_to_string(&mut text)
            .unwrap();
        text.lines()
            .map(|line| line.split('\t').take(2).collect::<Vec<_>>().join(":"))
            .collect()
    };

    let inputs = tempfile::tempdir().unwrap();
    write_indexed_vcf(
        inputs.path(),
        "danio_rerio.vcf.gz",
        &["21", "22"],
        "21\t300\trs2\tA\tG\t.\t.\tdbSNP_156;TSA=SNV\n\
         22\t100\trs7\tA\tG\t.\t.\tdbSNP_156;TSA=SNV\n\
         22\t500\trs8\tA\tG\t.\t.\tdbSNP_156;TSA=SNV\n",
    );
    let out = tempfile::tempdir().unwrap();
    let output = build("danio_rerio", "GRCz11", inputs.path(), out.path());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(output.status.success(), "{stderr}");
    let mut written: Vec<String> = std::fs::read_dir(out.path().join("variations"))
        .unwrap()
        .map(|entry| entry.unwrap().file_name().into_string().unwrap())
        .collect();
    written.sort();
    assert_eq!(written, ["21", "22"]);
    assert_eq!(names(out.path(), "21"), ["21:rs2"]);
    assert_eq!(names(out.path(), "22"), ["22:rs7", "22:rs8"]);
    assert_eq!(human_only_file_warnings(&stderr), 0, "{stderr}");

    let human_inputs = tempfile::tempdir().unwrap();
    write_indexed_vcf(
        human_inputs.path(),
        "homo_sapiens-chr21.vcf.gz",
        &["21"],
        "21\t300\trs2\tA\tG\t.\t.\tdbSNP_156;TSA=SNV\n",
    );
    let human_out = tempfile::tempdir().unwrap();
    let human = build(
        "homo_sapiens",
        "GRCh38",
        human_inputs.path(),
        human_out.path(),
    );
    let human_stderr = String::from_utf8_lossy(&human.stderr);
    assert!(human.status.success(), "{human_stderr}");
    assert_eq!(names(human_out.path(), "21"), ["21:rs2"]);
    assert_eq!(human_only_file_warnings(&human_stderr), 4, "{human_stderr}");

    let neither = tempfile::tempdir().unwrap();
    let refused_out = tempfile::tempdir().unwrap();
    let refused = build(
        "danio_rerio",
        "GRCz11",
        neither.path(),
        &refused_out.path().join("cache"),
    );
    assert_eq!(refused.status.code(), Some(1));
    let refusal = String::from_utf8_lossy(&refused.stderr);
    let errors: Vec<&str> = refusal
        .lines()
        .filter(|line| line.starts_with("ERROR: ["))
        .collect();
    assert_eq!(
        errors,
        [format!(
            "ERROR: [cache-builder] Failed to resolve variation inputs: No variation VCF of \
             danio_rerio in {}: looked for *-chr*.vcf.gz, one file a chromosome, and \
             danio_rerio.vcf.gz, the one genome-wide file Ensembl publishes for a species other \
             than homo_sapiens",
            neither.path().display()
        )]
    );
}

/// The build's variation lines name a sequence as the cache spells it, with
/// nothing put before the name. The genome-wide file of another species with
/// records on `mitochondrion_genome`, then on `2L`, is logged one `Reading the
/// variation records of chromosome <name>: <file>` line a sequence, in name
/// order, and a human chromosome file one line naming `21`. A genome-wide file
/// whose `2L` and `mitochondrion_genome` records are both out of position
/// order ends the build on one error line, the first sequence's in name order,
/// that names the sequence the same way.
#[test]
fn variation_lines_name_a_sequence_as_the_cache_spells_it() {
    let fixture = fixture_dir();
    let build = |species: &str, assembly: &str, vcf_dir: &Path, out: &Path| {
        Command::new(env!("CARGO_BIN_EXE_vep-cache-builder"))
            .args(["--release", "116", "--assembly", assembly])
            .args(["--species", species, "--include-variations"])
            .arg("--gff3")
            .arg(fixture.join("slice.gff3"))
            .arg("--fasta")
            .arg(fixture.join("peptides.fa"))
            .arg("--output-dir")
            .arg(out)
            .arg("--variation-vcf-dir")
            .arg(vcf_dir)
            .output()
            .expect("vep-cache-builder runs")
    };
    // A log line is a timestamp, one space and the rest; the rest is compared.
    let reading_lines = |stderr: &str| -> Vec<String> {
        stderr
            .lines()
            .filter(|line| line.contains("Reading the variation records"))
            .filter_map(|line| line.split_once(' '))
            .map(|(_, rest)| rest.to_string())
            .collect()
    };
    let cache_names = |cache: &Path| -> Vec<String> {
        let mut names: Vec<String> = std::fs::read_dir(cache.join("variations"))
            .unwrap()
            .map(|entry| entry.unwrap().file_name().into_string().unwrap())
            .collect();
        names.sort();
        names
    };

    let inputs = tempfile::tempdir().unwrap();
    write_indexed_vcf(
        inputs.path(),
        "drosophila_melanogaster.vcf.gz",
        &["mitochondrion_genome", "2L"],
        "mitochondrion_genome\t100\trs7\tA\tG\t.\t.\tdbSNP_156;TSA=SNV\n\
         2L\t300\trs2\tA\tG\t.\t.\tdbSNP_156;TSA=SNV\n",
    );
    let genome_vcf = inputs.path().join("drosophila_melanogaster.vcf.gz");
    let out = tempfile::tempdir().unwrap();
    let output = build(
        "drosophila_melanogaster",
        "BDGP6.54",
        inputs.path(),
        out.path(),
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(output.status.success(), "{stderr}");
    assert_eq!(
        reading_lines(&stderr),
        [
            format!(
                " INFO vep_cache_builder: Reading the variation records of chromosome 2L: {}",
                genome_vcf.display()
            ),
            format!(
                " INFO vep_cache_builder: Reading the variation records of chromosome \
                 mitochondrion_genome: {}",
                genome_vcf.display()
            ),
        ],
        "{stderr}"
    );
    assert_eq!(cache_names(out.path()), ["2L", "mitochondrion_genome"]);

    let human_inputs = tempfile::tempdir().unwrap();
    write_indexed_vcf(
        human_inputs.path(),
        "homo_sapiens-chr21.vcf.gz",
        &["21"],
        "21\t300\trs2\tA\tG\t.\t.\tdbSNP_156;TSA=SNV\n",
    );
    let human_out = tempfile::tempdir().unwrap();
    let human = build(
        "homo_sapiens",
        "GRCh38",
        human_inputs.path(),
        human_out.path(),
    );
    let human_stderr = String::from_utf8_lossy(&human.stderr);
    assert!(human.status.success(), "{human_stderr}");
    assert_eq!(
        reading_lines(&human_stderr),
        [format!(
            " INFO vep_cache_builder: Reading the variation records of chromosome 21: {}",
            human_inputs
                .path()
                .join("homo_sapiens-chr21.vcf.gz")
                .display()
        )],
        "{human_stderr}"
    );
    assert_eq!(cache_names(human_out.path()), ["21"]);

    let unordered = tempfile::tempdir().unwrap();
    let mut writer = bgzf::io::Writer::new(
        std::fs::File::create(unordered.path().join("drosophila_melanogaster.vcf.gz")).unwrap(),
    );
    writer.write_all(VCF_HEADER.as_bytes()).unwrap();
    writer
        .write_all(
            b"2L\t500\trs5\tA\tG\t.\t.\tdbSNP_156;TSA=SNV\n\
              2L\t600\trs6\tA\tG\t.\t.\tdbSNP_156;TSA=SNV\n\
              2L\t100\trs1\tA\tG\t.\t.\tdbSNP_156;TSA=SNV\n\
              mitochondrion_genome\t500\trs8\tA\tG\t.\t.\tdbSNP_156;TSA=SNV\n\
              mitochondrion_genome\t100\trs9\tA\tG\t.\t.\tdbSNP_156;TSA=SNV\n",
        )
        .unwrap();
    writer.try_finish().unwrap();
    let refused_out = tempfile::tempdir().unwrap();
    let refused = build(
        "drosophila_melanogaster",
        "BDGP6.54",
        unordered.path(),
        &refused_out.path().join("cache"),
    );
    assert_eq!(refused.status.code(), Some(1));
    let refusal = String::from_utf8_lossy(&refused.stderr);
    let errors: Vec<&str> = refusal
        .lines()
        .filter(|line| line.starts_with("ERROR: ["))
        .collect();
    assert_eq!(
        errors,
        [
            "ERROR: [cache-builder] Failed to write the variation file of chromosome 2L: rs1 \
             starts at 100, before the record written before it (500): the variation records of \
             chr 2L are not in position order"
        ],
        "{refusal}"
    );
}

/// The variation lines count as the rest of the log does: a count of one
/// takes its noun in the singular, any other the plural. A human chromosome
/// file of one record is logged as one chromosome, one variation, one
/// variation record and one record; a genome-wide file of two records on `2L`
/// and one on `mitochondrion_genome` as two chromosomes and three records,
/// each sequence's lines by its own count.
#[test]
fn variation_lines_count_one_in_the_singular() {
    let fixture = fixture_dir();
    // The lines of a build that count variations, records or chromosomes,
    // each without its timestamp.
    let counting_lines = |species: &str, assembly: &str, vcf_dir: &Path, out: &Path| {
        let output = Command::new(env!("CARGO_BIN_EXE_vep-cache-builder"))
            .args(["--release", "116", "--assembly", assembly])
            .args(["--species", species, "--include-variations"])
            .arg("--gff3")
            .arg(fixture.join("slice.gff3"))
            .arg("--fasta")
            .arg(fixture.join("peptides.fa"))
            .arg("--output-dir")
            .arg(out)
            .arg("--variation-vcf-dir")
            .arg(vcf_dir)
            .output()
            .expect("vep-cache-builder runs");
        let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
        assert!(output.status.success(), "{stderr}");
        stderr
            .lines()
            .filter_map(|line| line.split_once(' '))
            .map(|(_, rest)| rest.to_string())
            .filter(|rest| {
                rest.starts_with(" INFO vep_cache_builder::variation: ")
                    || rest.contains("Resolved the variation VCFs")
                    || rest.contains("Variation files written")
            })
            .collect::<Vec<String>>()
    };

    let human_inputs = tempfile::tempdir().unwrap();
    write_indexed_vcf(
        human_inputs.path(),
        "homo_sapiens-chr21.vcf.gz",
        &["21"],
        "21\t300\trs2\tA\tG\t.\t.\tdbSNP_156;TSA=SNV\n",
    );
    let human_vcf = human_inputs.path().join("homo_sapiens-chr21.vcf.gz");
    let human_out = tempfile::tempdir().unwrap();
    assert_eq!(
        counting_lines(
            "homo_sapiens",
            "GRCh38",
            human_inputs.path(),
            human_out.path()
        ),
        [
            " INFO vep_cache_builder: Resolved the variation VCFs of 1 chromosome".to_string(),
            format!(
                " INFO vep_cache_builder::variation: Parsed 1 variation from {} (0 skipped)",
                human_vcf.display()
            ),
            " INFO vep_cache_builder::variation: chr 21: 1 variation record".to_string(),
            format!(
                " INFO vep_cache_builder::variation: chr 21: wrote {} (1 record)",
                human_out.path().join("variations/21/all_vars.gz").display()
            ),
            " INFO vep_cache_builder: Variation files written: 1 record across 1 chromosome"
                .to_string(),
        ]
    );

    let inputs = tempfile::tempdir().unwrap();
    write_indexed_vcf(
        inputs.path(),
        "drosophila_melanogaster.vcf.gz",
        &["mitochondrion_genome", "2L"],
        "mitochondrion_genome\t100\trs7\tA\tG\t.\t.\tdbSNP_156;TSA=SNV\n\
         2L\t300\trs2\tA\tG\t.\t.\tdbSNP_156;TSA=SNV\n\
         2L\t400\trs3\tA\tG\t.\t.\tdbSNP_156;TSA=SNV\n",
    );
    let genome_vcf = inputs.path().join("drosophila_melanogaster.vcf.gz");
    let out = tempfile::tempdir().unwrap();
    let wrote = |name: &str, records: &str| {
        format!(
            " INFO vep_cache_builder::variation: chr {name}: wrote {} ({records})",
            out.path()
                .join("variations")
                .join(name)
                .join("all_vars.gz")
                .display()
        )
    };
    // Two chromosomes are written at once, so their lines are compared as a set.
    let mut lines = counting_lines(
        "drosophila_melanogaster",
        "BDGP6.54",
        inputs.path(),
        out.path(),
    );
    lines.sort();
    let mut expected = [
        " INFO vep_cache_builder: Resolved the variation VCFs of 2 chromosomes".to_string(),
        format!(
            " INFO vep_cache_builder::variation: Parsed 2 variations from {} (0 skipped)",
            genome_vcf.display()
        ),
        " INFO vep_cache_builder::variation: chr 2L: 2 variation records".to_string(),
        wrote("2L", "2 records"),
        format!(
            " INFO vep_cache_builder::variation: Parsed 1 variation from {} (0 skipped)",
            genome_vcf.display()
        ),
        " INFO vep_cache_builder::variation: chr mitochondrion_genome: 1 variation record"
            .to_string(),
        wrote("mitochondrion_genome", "1 record"),
        " INFO vep_cache_builder: Variation files written: 3 records across 2 chromosomes"
            .to_string(),
    ];
    expected.sort();
    assert_eq!(lines, expected);
}

/// Ensembl VEP's dumper reads a VCF through ensembl-io's VCF parser
/// (`Dumper/Variation.pm:307`), which ends the run on a `##fileformat` line
/// naming a major version other than 4 (`BaseVCF4.pm:125`): a chromosome file
/// headed `VCFv3.3`, read through, and a frequency file headed so, read through
/// its index, each end the build with that text and no variation file; the
/// files headed `VCFv4.5` build.
#[test]
fn a_variation_vcf_whose_fileformat_is_not_vcf_4_ends_the_build() {
    let build = |chromosome: &str, frequency: &str| {
        build_with_headers(&headed(chromosome), &headed(frequency))
    };
    let refused = |stderr: &str, file: &str| {
        stderr.contains("Failed to read VCF header: ")
            && stderr.contains(&format!(
                "{file}: The VCF file format version 3.3 is not compatible with the parser \
                 version (VCF v4.2)"
            ))
    };

    let (ok, stderr, written) = build("VCFv3.3", "VCFv4.2");
    assert!(!ok, "{stderr}");
    assert!(refused(&stderr, "homo_sapiens-chr21.vcf.gz"), "{stderr}");
    assert!(written.is_none(), "a variation file was written: {stderr}");

    let (ok, stderr, written) = build("VCFv4.2", "VCFv3.3");
    assert!(!ok, "{stderr}");
    assert!(refused(&stderr, "1000GENOMES-phase_3.vcf.gz"), "{stderr}");
    assert!(written.is_none(), "a variation file was written: {stderr}");

    let (ok, stderr, written) = build("VCFv4.5", "VCFv4.5");
    assert!(ok, "{stderr}");
    assert_eq!(written.as_deref(), Some(JOINED_RECORD), "{stderr}");
}

/// The dumper's parser refuses no header line but a `##fileformat` line
/// (`BaseVCF4.pm:116-130`), reads the version as the line's first
/// `<digits>.<digits>` and prints it as the line spelt it: a header without
/// the line builds with its records read, as do `VCFv4.2.1` (4.2) and
/// `VCFv04.2`; `VCFv4` ends the build with the unreadable-version text and
/// `VCFv05.0` with the incompatible-version text naming `05.0`; each shape on
/// the chromosome file, read through, and on the frequency file, read through
/// its index.
#[test]
fn a_variation_vcf_header_reads_as_ensembl_vep_reads_it_on_every_fileformat_shape() {
    let no_line = VCF_HEADER.replacen("##fileformat=VCFv4.2\n", "", 1);
    for header in [no_line, headed("VCFv4.2.1"), headed("VCFv04.2")] {
        for (chromosome, frequency) in
            [(header.as_str(), VCF_HEADER), (VCF_HEADER, header.as_str())]
        {
            let (ok, stderr, written) = build_with_headers(chromosome, frequency);
            assert!(ok, "{header:?}: {stderr}");
            assert_eq!(
                written.as_deref(),
                Some(JOINED_RECORD),
                "{header:?}: {stderr}"
            );
        }
    }
    for (version, reason) in [
        (
            "VCFv4",
            "The script can't read the VCF file format version of 'fileformat'",
        ),
        (
            "VCFv05.0",
            "The VCF file format version 05.0 is not compatible with the parser version \
             (VCF v4.2)",
        ),
    ] {
        let header = headed(version);
        for (chromosome, frequency, file) in [
            (header.as_str(), VCF_HEADER, "homo_sapiens-chr21.vcf.gz"),
            (VCF_HEADER, header.as_str(), "1000GENOMES-phase_3.vcf.gz"),
        ] {
            let (ok, stderr, written) = build_with_headers(chromosome, frequency);
            assert!(!ok, "{version}: {stderr}");
            assert!(
                stderr.contains("Failed to read VCF header: ")
                    && stderr.contains(&format!("{file}: {reason}")),
                "{version} on {file}: {stderr}"
            );
            assert!(written.is_none(), "{version}: a variation file was written");
        }
    }
}

/// [`VCF_HEADER`] with its `##fileformat` line naming `version`.
fn headed(version: &str) -> String {
    VCF_HEADER.replacen("VCFv4.2", version, 1)
}

/// The one record [`build_with_headers`] writes when both files are read: the
/// chromosome file's SNV with the frequency file's populations joined.
const JOINED_RECORD: &str =
    "21\trs2\t.\t.\t300\t.\tA/G\t.\t.\t.\t.\t.\t.\t.\t.\tG:0.0001997\tG:0.0008\tG:0\tG:0\tG:0\tG:0\n";

/// Builds with a one-record chromosome file headed `chromosome`, read through,
/// and the matching frequency file headed `frequency`, read through its index:
/// whether the build succeeded, its standard error, and the chromosome's
/// variation file when one was written.
fn build_with_headers(chromosome: &str, frequency: &str) -> (bool, String, Option<String>) {
    let fixture = fixture_dir();
    let record = "21\t300\trs2\tA\tG\t.\t.\tdbSNP_156;TSA=SNV\n";
    let frequencies =
        "21\t300\trs2\tA\tG\t.\t.\tdbSNP_156;TSA=SNV;AFR=0.0008;AMR=0;EAS=0;EUR=0;SAS=0\n";
    let inputs = tempfile::tempdir().unwrap();
    for (name, header, records) in [
        ("homo_sapiens-chr21.vcf.gz", chromosome, record),
        ("1000GENOMES-phase_3.vcf.gz", frequency, frequencies),
    ] {
        write_indexed_vcf_with_header(inputs.path(), name, header, &["21"], records);
    }
    let out = tempfile::tempdir().unwrap();
    let output = Command::new(env!("CARGO_BIN_EXE_vep-cache-builder"))
        .args([
            "--release",
            "116",
            "--assembly",
            "GRCh38",
            "--include-variations",
        ])
        .arg("--gff3")
        .arg(fixture.join("slice.gff3"))
        .arg("--fasta")
        .arg(fixture.join("peptides.fa"))
        .arg("--output-dir")
        .arg(out.path())
        .arg("--variation-vcf-dir")
        .arg(inputs.path())
        .output()
        .expect("vep-cache-builder runs");
    let all_vars = out.path().join("variations/21/all_vars.gz");
    let written = all_vars.exists().then(|| {
        let mut text = String::new();
        bgzf::io::Reader::new(std::fs::File::open(&all_vars).unwrap())
            .read_to_string(&mut text)
            .unwrap();
        text
    });
    (
        output.status.success(),
        String::from_utf8_lossy(&output.stderr).into_owned(),
        written,
    )
}
