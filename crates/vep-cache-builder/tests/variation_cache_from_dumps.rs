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
    let path = dir.join(name);
    let mut writer = bgzf::io::Writer::new(std::fs::File::create(&path).unwrap());
    writer.write_all(VCF_HEADER.as_bytes()).unwrap();
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
        "1\t100\tCOSV9900001\tT\t.\t.\t.\tCOSMIC_102;TSA=SNV;E_Phenotype_or_Disease\n\
         21\t200\tCOSV1\tC\t.\t.\t.\tCOSMIC_102;TSA=SNV;E_Phenotype_or_Disease\n",
    );
    write_dumps(inputs.path());
    let out = tempfile::tempdir().unwrap();

    let status = Command::new(env!("CARGO_BIN_EXE_vep-cache-builder"))
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
        .status()
        .expect("vep-cache-builder runs");
    assert!(status.success(), "vep-cache-builder exited {status}");

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
