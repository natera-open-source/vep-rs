// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Unit tests of `classes`: the variant class of a row, the family classes, the genotype
//! classes of a sample and of a record, and the regulatory record classes.

use super::*;

fn set(names: &[&'static str]) -> BTreeSet<&'static str> {
    names.iter().copied().collect()
}

fn refusal(message: &str) -> Failure {
    failure(message)
}

#[test]
fn variant_class_inference() {
    let class = |location: &str, allele: &str| variant_class(location, allele, false);
    assert_eq!(class("21:100", "G"), "snv");
    assert_eq!(class("21:101", "-"), "deletion");
    assert_eq!(class("21:100-101", "C"), "insertion");
    assert_eq!(class("21:100-105", "-"), "deletion");
    assert_eq!(class("21:100-105", "deletion"), "symbolic_deletion");
    assert_eq!(
        class("21:100-105", "Alu_insertion"),
        "symbolic_Alu_insertion"
    );
    assert_eq!(class("21:100-101", "TG"), "insertion");
    assert_eq!(class("21:100", "TG"), "mnv_or_complex");
    assert_eq!(class("21:100", "A[21:200["), "breakend");
    // The mate's contig name is not inserted sequence.
    assert_eq!(class("21:100", "]MT:200]G"), "breakend");
    assert_eq!(class("21:100", "[21:200[TTGA"), "breakend_insertion");
    assert_eq!(class("21:100", ".G"), "single_breakend");
}

#[test]
fn variant_class_reads_a_substitution_and_coordinates_of_any_length() {
    assert_eq!(variant_class("21:100-101", "TG", true), "mnv_or_complex");
    assert_eq!(variant_class("21:100-102", "TGA", true), "mnv_or_complex");
    assert_eq!(variant_class("21:100", "G", true), "snv");
    let class = |location: &str, allele: &str| variant_class(location, allele, false);
    // A span is read after the first `:`, or whole without one; numbers are compared as
    // numbers.
    assert_eq!(class("100-101", "TG"), "insertion");
    assert_eq!(class("HLA:x:100-101", "TG"), "mnv_or_complex");
    assert_eq!(class("21:0100-101", "TG"), "insertion");
    assert_eq!(class("21:99-100", "TG"), "insertion");
    assert_eq!(class("21:0-1", "TG"), "insertion");
    assert_eq!(class("21:100-100", "TG"), "mnv_or_complex");
    assert_eq!(class("21:101-100", "TG"), "mnv_or_complex");
    assert_eq!(class("21:100-101-102", "TG"), "mnv_or_complex");
    assert_eq!(class("21:-101", "TG"), "mnv_or_complex");
    assert_eq!(
        class(
            "21:99999999999999999999999999999-100000000000000000000000000000",
            "TG"
        ),
        "insertion"
    );
    assert_eq!(class("21:100", ""), "mnv_or_complex");
    assert_eq!(class("21:100-101", ""), "insertion");
    assert_eq!(class("21:100", "G."), "single_breakend");
    assert_eq!(class("21:100", "*"), "snv");
    assert_eq!(class("21:100", "acgtn"), "mnv_or_complex");
    assert_eq!(class("21:100", "A]1:5]TT]1:9]"), "breakend_insertion");
    assert_eq!(class("21:100", "]1:5]A]1:9]"), "breakend");
    assert_eq!(class("21:100", "AT]1:5]"), "breakend_insertion");
    // One bracket closes nothing: every character stays.
    assert_eq!(class("21:100", "A["), "breakend_insertion");
    assert_eq!(class("21:100", "["), "breakend");
    assert_eq!(class("21:100", "\u{e9}"), "symbolic_\u{e9}");
}

#[test]
fn family_columns_are_names_values_custom_matches_and_agreements() {
    let specs = |arg: &str| parse_family_columns(Some(arg));
    assert_eq!(
        specs("HGVSc,custom_match=CE/CO"),
        Ok(vec!["HGVSc".to_string(), "custom_match=CE/CO".to_string()])
    );
    assert_eq!(
        specs("REFSEQ_MATCH=,SOURCE=,BAM_EDIT=,REFSEQ_OFFSET,GIVEN_REF/USED_REF").unwrap(),
        [
            "REFSEQ_MATCH=",
            "SOURCE=",
            "BAM_EDIT=",
            "REFSEQ_OFFSET",
            "GIVEN_REF/USED_REF"
        ]
    );
    assert_eq!(specs(" HGVSc , ,IMPACT= ").unwrap(), ["HGVSc", "IMPACT="]);
    assert_eq!(parse_family_columns(None), Ok(Vec::new()));
    assert_eq!(specs(""), Ok(Vec::new()));
    for bad in [
        "custom_match=CE",
        "custom_match=CE/CO/X",
        "custom_match=/CO",
        "GIVEN_REF/",
        "/USED_REF",
        "A/B/C",
        "GIVEN_REF/USED_REF=",
        "A B",
        "A==",
        "=",
    ] {
        assert_eq!(
            specs(bad),
            Err(refusal(&format!(
                "--family-columns wants NAME, NAME=, custom_match=EXACT/OVERLAP or A/B, got '{bad}'"
            ))),
            "{bad}"
        );
    }
    assert_eq!(
        specs("HGVSc,it's"),
        Err(refusal(
            "--family-columns wants NAME, NAME=, custom_match=EXACT/OVERLAP or A/B, got \"it's\""
        ))
    );
}

#[test]
fn custom_match_classes_follow_the_two_sources() {
    let class = |extra: &str, location: &str, allele: &str| {
        custom_match_class(&parse_extra(extra), "CE", "CO", location, allele)
    };
    let snv = |extra: &str| class(extra, "21:100", "G");
    assert_eq!(snv("CE=exact_1;CO=exact_1;CE_AF=0.10;CO_AF=0.10"), "exact");
    assert_eq!(snv("CO=overlap_1;CO_AF=0.20"), "overlap");
    assert_eq!(snv("IMPACT=MODIFIER"), "miss");
    assert_eq!(
        snv("CE=multi_1;CO=multi_1;CE_AF=0.30;CO_AF=0.30,0.40"),
        "multi_allelic"
    );
    assert_eq!(
        class(
            "CO=exact_1,overlap_1;CO_AF=0.10,0.20",
            "21:500-900",
            "deletion"
        ),
        "sv_span"
    );
    assert_eq!(class("CO=sv_1", "21:500", "A[21:900["), "sv_span");
    // An exact structural match is exact.
    assert_eq!(class("CE=sv_1;CO=sv_1", "21:500-900", "deletion"), "exact");
    // Two records, one value each.
    assert_eq!(
        snv("CE=exact_2;CO=exact_2,overlap_2;CE_AF=0.5;CO_AF=0.5,0.6"),
        "exact"
    );
    // FILTER lists filters, not alleles.
    assert_eq!(
        snv("CE=exact_3;CO=exact_3;CE_FILTER=q10,s50;CO_FILTER=q10,s50"),
        "exact"
    );
    // A comma the exact source lists whole.
    assert_eq!(
        snv("CE=exact_4;CO=exact_4;CE_CLNSIG=a,b;CO_CLNSIG=a,b"),
        "exact"
    );
    // A dash is an absent value.
    assert_eq!(snv("CE=-;CO=overlap_5"), "overlap");
    // A field the column gives twice counts with its last value.
    assert_eq!(snv("CE=m;CO=m;CO_AF=1,2;CE_AF=1;CO_AF=1"), "exact");
    assert_eq!(
        CUSTOM_MATCH_CLASSES,
        ["exact", "multi_allelic", "overlap", "sv_span", "miss"]
    );
}

#[test]
fn family_classes_name_one_class_a_spec() {
    assert_eq!(
        family_classes(
            &parse_extra("CO=o;IMPACT=HIGH"),
            &["custom_match=CE/CO", "IMPACT="],
            "21:100",
            "G"
        ),
        ["custom_match=overlap", "IMPACT=HIGH"]
    );
    let agree = |extra: &str| agreement_class(&parse_extra(extra), "GIVEN_REF", "USED_REF");
    assert_eq!(agree("GIVEN_REF=A;USED_REF=A"), "same");
    assert_eq!(agree("GIVEN_REF=A;USED_REF=G"), "differ");
    assert_eq!(agree("IMPACT=MODIFIER"), "-");
    assert_eq!(agree("GIVEN_REF=A"), "-");
    assert_eq!(agree("GIVEN_REF=A;USED_REF=-"), "-");
    assert_eq!(AGREEMENT_CLASSES, ["same", "differ", "-"]);
    let specs = [
        "REFSEQ_MATCH=",
        "SOURCE=",
        "BAM_EDIT=",
        "REFSEQ_OFFSET",
        "GIVEN_REF/USED_REF",
    ];
    let raw = "REFSEQ_MATCH=rseq_mrna_nonmatch,rseq_5p_mismatch;BAM_EDIT=OK;GIVEN_REF=C;USED_REF=T";
    let expected = [
        "REFSEQ_MATCH=rseq_mrna_nonmatch,rseq_5p_mismatch",
        "SOURCE=-",
        "BAM_EDIT=OK",
        "REFSEQ_OFFSET-",
        "GIVEN_REF/USED_REF=differ",
    ];
    assert_eq!(family_classes(&parse_extra(raw), &specs, "", ""), expected);
    let read: Vec<FamilySpec> = specs.iter().map(|spec| FamilySpec::read(spec)).collect();
    assert_eq!(row_family_classes(raw, &read, "", ""), expected);
    // An empty value and a dash: absent for presence, `-` for the value of an empty one.
    assert_eq!(
        family_classes(
            &parse_extra("A=;B=-;C"),
            &["A", "A=", "B", "B=", "C", "C=", "D", "D="],
            "",
            ""
        ),
        ["A-", "A=-", "B-", "B=-", "C-", "C=-", "D-", "D=-"]
    );
    // Any text is a spec: one of no other form asks for a column of its name.
    assert_eq!(
        family_classes(&parse_extra("A/B/C=1"), &["A/B/C", "x y="], "", ""),
        ["A/B/C+", "x y=-"]
    );
    assert_eq!(
        FamilySpec::read("custom_match/x"),
        FamilySpec::Agreement {
            spec: "custom_match/x".to_string(),
            a: "custom_match".to_string(),
            b: "x".to_string()
        }
    );
}

#[test]
fn genotype_labels_follow_vep_genotype_reading() {
    // Each GT's classes follow how Ensembl VEP reads a genotype for `--individual`: the bits
    // split on `|` when the GT has one, else on `/`, `.` bits dropped before translation,
    // so `./1` is a one-allele genotype and `./.` yields nothing, an allele containing `*`
    // outside the non-reference set, a haploid `1` homozygous, and an index read as a
    // number (`00` is `0`).
    let labels = |gt: &str, alts: &[&str]| genotype_labels(gt, alts, "r").unwrap();
    assert_eq!(labels("0|1", &["G"]), set(&["het"]));
    assert_eq!(labels("1/0", &["G"]), set(&["het"]));
    assert_eq!(labels("1|1", &["G"]), set(&["hom_alt"]));
    assert_eq!(labels("0/0", &["G"]), set(&["hom_ref"]));
    assert_eq!(labels("./.", &["G"]), set(&["missing"]));
    assert_eq!(labels(".", &["G"]), set(&["missing"]));
    assert_eq!(
        labels("./1", &["G"]),
        set(&["partially_missing", "hom_alt"])
    );
    assert_eq!(
        labels(".|0", &["G"]),
        set(&["partially_missing", "hom_ref"])
    );
    assert_eq!(labels("1", &["G"]), set(&["haploid", "hom_alt"]));
    assert_eq!(labels("0", &["G"]), set(&["haploid", "hom_ref"]));
    assert_eq!(labels("1/2", &["G", "T"]), set(&["multi_allelic"]));
    assert_eq!(labels("0/1/2", &["G", "T"]), set(&["multi_allelic"]));
    assert_eq!(labels("1/1/1", &["G", "T"]), set(&["hom_alt"]));
    assert_eq!(labels("0/00", &["G"]), set(&["hom_ref"]));
    assert_eq!(labels("01/1", &["G", "T"]), set(&["hom_alt"]));
    assert_eq!(labels("0/2", &["G", "*"]), set(&["star", "hom_ref"]));
    assert_eq!(labels("2/2", &["G", "*"]), set(&["star", "hom_ref"]));
    assert_eq!(labels("1/2", &["G", "*"]), set(&["star", "het"]));
    assert_eq!(
        labels("2", &["G", "*"]),
        set(&["haploid", "star", "hom_ref"])
    );
    assert_eq!(
        labels("./2", &["G", "*"]),
        set(&["partially_missing", "star", "hom_ref"])
    );
    assert_eq!(labels("0/2", &["G", "<*>"]), set(&["star", "hom_ref"]));
    for (gt, repr) in [
        ("0/3", "'0/3'"),
        ("a/b", "'a/b'"),
        ("0/-1", "'0/-1'"),
        ("1/", "'1/'"),
        ("", "''"),
        ("1\\2", "'1\\\\2'"),
        ("0/1|2", "'0/1|2'"),
        ("0/+1", "'0/+1'"),
        (
            "0/99999999999999999999999999",
            "'0/99999999999999999999999999'",
        ),
    ] {
        assert_eq!(
            genotype_labels(gt, &["G", "T"], "21 100 rs1 A G,T"),
            Err(refusal(&format!(
                "record 21 100 rs1 A G,T has a GT {repr} outside its alleles"
            ))),
            "{gt}"
        );
    }
}

#[test]
fn genotype_labels_refuse_a_non_ascii_digit() {
    // A GT bit is an ASCII allele index: a superscript two, a circled one and an
    // Arabic-Indic three are each refused as a bit outside the record's alleles, through
    // `genotype_labels` and through a record's kept sample.
    let record = "21 100 rs1 A G,T,C";
    for gt in ["0/\u{b2}", "\u{2460}", "\u{663}/0"] {
        assert_eq!(
            genotype_labels(gt, &["G", "T", "C"], record),
            Err(refusal(&format!(
                "record {record} has a GT '{gt}' outside its alleles"
            ))),
            "{gt}"
        );
    }
    let fields = [
        "21", "100", "rs1", "A", "G,T,C", ".", "PASS", ".", "GT", "0/0", "0/\u{b2}",
    ];
    assert_eq!(
        record_genotype_classes(&fields, 2),
        Err(refusal(&format!(
            "record {record} has a GT '0/\u{b2}' outside its alleles"
        )))
    );
}

/// The columns of a record: the eight fixed ones of `21 100 rs1 A G`, then `rest`.
fn with_fixed(rest: &[&'static str]) -> Vec<&'static str> {
    let mut fields = vec!["21", "100", "rs1", "A", "G", ".", "PASS", "."];
    fields.extend(rest);
    fields
}

#[test]
fn record_genotype_classes_add_the_record_shapes() {
    // A record shows every class any kept sample shows, plus the shapes of the whole kept
    // set: every genotype missing, no non-reference carrier, and one reference-homozygous
    // sample with the rest missing. GT is read at its FORMAT index, a bare `.` sample is
    // missing, and columns past the kept ones are not read. A record Ensembl VEP reads as a
    // structural variant shows no class; `<*>` is no such allele, and a bare `SVTYPE` flag
    // or `SVTYPE=` is no such value.
    let classes = |fields: &[&str], keep: usize| record_genotype_classes(fields, keep);
    let structural: [&[&str]; 5] = [
        &[
            "21", "100", "sv1", "A", "<CN0>", ".", "PASS", ".", "GT", "0|1", "1|1",
        ],
        &[
            "21",
            "100",
            "sv2",
            "N",
            "<INS:ME:ALU>",
            ".",
            "PASS",
            "SVTYPE=INS",
            "GT",
            "0/1",
            "0/0",
        ],
        &[
            "21",
            "100",
            "sv3",
            "A",
            "N",
            ".",
            "PASS",
            "SVTYPE=DEL;END=5000",
            "GT",
            "0/1",
            "0/0",
        ],
        &[
            "21",
            "100",
            "b1",
            "A",
            "A[21:200[",
            ".",
            "PASS",
            ".",
            "GT",
            "0/1",
            "0/0",
        ],
        &[
            "21", "100", "b2", "A", ".G", ".", "PASS", ".", "GT", "0/1", "0/0",
        ],
    ];
    for fields in structural {
        assert_eq!(classes(fields, 2), Ok(BTreeSet::new()), "{fields:?}");
    }
    // Plain bases: no structural variant.
    assert_eq!(
        classes(
            &[
                "21",
                "100",
                "sv5",
                "A",
                "G",
                ".",
                "PASS",
                "SVTYPE=DEL",
                "GT",
                "0/1",
                "0/0"
            ],
            2
        ),
        Ok(set(&["het", "hom_ref"]))
    );
    // An empty SVTYPE is none.
    assert_eq!(
        classes(
            &["21", "100", "sv6", "N", "NC", ".", "PASS", "SVTYPE=", "GT", "0/1", "0/0"],
            2
        ),
        Ok(set(&["het", "hom_ref"]))
    );
    // A bare SVTYPE flag is none.
    assert_eq!(
        classes(
            &["21", "100", "sv4", "N", "NC", ".", "PASS", "SVTYPE", "GT", "0/1", "0/0"],
            2
        ),
        Ok(set(&["het", "hom_ref"]))
    );
    assert_eq!(
        classes(
            &["21", "100", "g1", "A", "G,<*>", ".", "PASS", ".", "GT", "0/2", "1/2"],
            2
        ),
        Ok(set(&["star", "hom_ref", "het"]))
    );
    assert_eq!(
        classes(&with_fixed(&["GT", "./.", "./.", "./."]), 3),
        Ok(set(&["missing", "all_missing"]))
    );
    assert_eq!(
        classes(&with_fixed(&["GT", "0/0", "./.", "./."]), 3),
        Ok(set(&[
            "hom_ref",
            "missing",
            "no_carrier",
            "one_hom_ref_rest_missing"
        ]))
    );
    assert_eq!(
        classes(&with_fixed(&["GT", "0/0", "0/0", "./."]), 3),
        Ok(set(&["hom_ref", "missing", "no_carrier"]))
    );
    assert_eq!(
        classes(&with_fixed(&["GT", "0/0", "0/1", "./."]), 3),
        Ok(set(&["hom_ref", "het", "missing"]))
    );
    // The third sample is not kept.
    assert_eq!(
        classes(&with_fixed(&["GT:DP", "0/0:9", "0/0:9", "0/1:9"]), 2),
        Ok(set(&["hom_ref", "no_carrier"]))
    );
    assert_eq!(
        classes(&with_fixed(&["DP:GT", "9:1", ".:.", ".", "9:0"]), 4),
        Ok(set(&["haploid", "hom_alt", "hom_ref", "missing"]))
    );
    assert_eq!(
        classes(&with_fixed(&["GT", "0/0", "0/0"]), 2),
        Ok(set(&["hom_ref", "no_carrier"]))
    );
    assert_eq!(
        classes(&with_fixed(&["DP", "9", "9"]), 2),
        Err(refusal("record 21 100 rs1 A G has no GT in its FORMAT"))
    );
    assert_eq!(
        classes(&with_fixed(&["GT", "0/0"]), 2),
        Err(refusal(
            "record 21 100 rs1 A G has 1 sample columns, 2 are kept"
        ))
    );
    assert_eq!(
        classes(&with_fixed(&[]), 1),
        Err(refusal(
            "record 21 100 rs1 A G has 0 sample columns, 1 are kept"
        ))
    );
}

#[test]
fn record_genotype_classes_of_no_kept_sample_and_of_a_short_record() {
    // No kept sample: a record with a FORMAT column shows `all_missing`, one without is
    // refused; a record of fewer than five columns is named by the ones it has.
    assert_eq!(
        record_genotype_classes(&with_fixed(&["GT", "0/1"]), 0),
        Ok(set(&["all_missing"]))
    );
    assert_eq!(
        record_genotype_classes(&with_fixed(&[]), 0),
        Err(refusal(
            "record 21 100 rs1 A G has 0 sample columns, 0 are kept"
        ))
    );
    assert_eq!(
        record_genotype_classes(&["21", "100"], 1),
        Err(refusal("record 21 100 has 0 sample columns, 1 are kept"))
    );
    let owned: Vec<String> = with_fixed(&["GT", "1/1", "0/1"])
        .into_iter()
        .map(str::to_string)
        .collect();
    assert_eq!(
        record_genotype_classes(&owned, 2),
        Ok(set(&["hom_alt", "het"]))
    );
}

#[test]
fn structural_record_reads_svtype_as_ensembl_vep_does() {
    // `SVTYPE` makes a record structural only with a value Perl reads as true: a bare flag
    // and `SVTYPE=` store nothing, `SVTYPE=0` stores Perl's false `0`, and of a repeated key
    // the last one counts; `SVTYPE=0.0` is true. Plain `ACGT` alleles are never structural,
    // and the allele test runs on the comma-joined ALT string, so a single breakend counts
    // at that string's ends only.
    let structural =
        |ref_allele: &str, alts: &[&str], info: &str| is_structural_record(ref_allele, alts, info);
    for info in [
        "SVTYPE",
        "END=5;SVTYPE",
        "SVTYPE=",
        "SVTYPE==",
        "SVTYPE=0",
        "SVTYPE=0=X",
        "SVTYPE=DEL;SVTYPE=0",
        ".",
    ] {
        assert!(!structural("N", &["NC"], info), "{info}");
    }
    for info in [
        "SVTYPE=DEL",
        "SVTYPE=DEL=X",
        "SVTYPE=0.0",
        "SVTYPE=0;SVTYPE=DEL",
        "END=5;SVTYPE=INS",
    ] {
        assert!(structural("N", &["NC"], info), "{info}");
    }
    assert!(!structural("A", &["G"], "SVTYPE=DEL"));
    assert!(structural("A", &["<CN0>"], ".") && !structural("A", &["<*>"], "."));
    assert!(structural("A", &[".T", "G"], ".") && !structural("A", &["G", ".T"], "."));
    assert!(structural("A", &["G", "T."], ".") && !structural("A", &["T.", "G"], "."));
}

#[test]
fn structural_alleles_are_read_on_the_joined_alt() {
    let structural = |alts: &[&str]| is_structural_record("N", alts, ".");
    assert!(structural(&["N[21:5["]) && structural(&["]21:5]N"]));
    // A bracket pair needs a character between its ends, and none of them a `*`.
    assert!(!structural(&["<>"]) && !structural(&["<a*b>"]));
    assert!(structural(&["<a", "b>"]));
    assert!(structural(&["._"]) && structural(&["\u{e9}."]) && structural(&[".5"]));
    assert!(!structural(&["."]) && !structural(&[".-"]) && !structural(&["-."]));
    assert!(structural(&["A.\n"]) && !structural(&["A.\n\n"]));
    // Without bases the alleles are not plain bases; lower case is not plain either.
    assert!(!is_structural_record("", &[""], "."));
    assert!(is_structural_record("", &[""], "SVTYPE=DEL"));
    assert!(is_structural_record("a", &["g"], "SVTYPE=DEL"));
    assert!(!is_structural_record("A", &["", "G"], "SVTYPE=DEL"));
}

#[test]
fn class_names_of_a_command_line_are_held_to_the_lists() {
    let names = |given: &[&str]| {
        genotype_class_names(
            &given
                .iter()
                .map(|name| name.to_string())
                .collect::<Vec<_>>(),
        )
    };
    assert_eq!(names(&["het", "star", "het"]).unwrap(), ["het", "star"]);
    assert_eq!(names(&[]), Ok(Vec::new()));
    assert_eq!(
        names(&["het", "homozygous"]),
        Err(refusal(
            "--require-genotype-class wants one of het, hom_alt, hom_ref, missing, \
             partially_missing, haploid, star, multi_allelic, all_missing, no_carrier, \
             one_hom_ref_rest_missing, got 'homozygous'"
        ))
    );
    assert_eq!(record_class_key(None), Ok(None));
    assert_eq!(
        record_class_key(Some("regulatory_class")),
        Ok(Some("regulatory_class"))
    );
    assert_eq!(
        record_class_key(Some("splice_class")),
        Err(refusal(
            "--record-classes wants one of regulatory_class, got 'splice_class'"
        ))
    );
    assert_eq!(RECORD_CLASS_RESERVOIR, 16);
    let keys: Vec<&str> = RECORD_CLASS_FUNCTIONS.iter().map(|(key, _)| *key).collect();
    assert_eq!(keys, RECORD_CLASS_KEYS);
    assert!(record_class_function("regulatory_class").is_some());
    assert!(record_class_function("splice_class").is_none());
}

fn motif_row(location: &str, allele: &str) -> Row {
    Row {
        dataset: "reg".to_string(),
        uploaded_variation: "x".to_string(),
        location: location.to_string(),
        allele: allele.to_string(),
        feature: "ENSM4".to_string(),
        feature_type: "MotifFeature".to_string(),
        consequence_set: "TF_binding_site_variant".to_string(),
        strand: "1".to_string(),
        flags: String::new(),
        positions: ["-", "-", "-"].map(str::to_string),
        family: Vec::new(),
        substitution: false,
    }
}

fn pairs(alleles: &[(&str, &str)]) -> Vec<(String, String)> {
    alleles
        .iter()
        .map(|(ref_allele, alt)| (ref_allele.to_string(), alt.to_string()))
        .collect()
}

#[test]
fn allele_shape_reads_the_input_pair_the_allele_came_from() {
    // A two-base substitution and a two-base insertion print one row shape: the input's
    // alleles tell them apart, and a row matching no pair reads as the row alone does.
    let mnv_row = motif_row("21:5000-5001", "GC");
    assert_eq!(allele_shape(&mnv_row, &pairs(&[("AT", "GC")])), "mnv");
    assert_eq!(allele_shape(&mnv_row, &pairs(&[("A", "AGC")])), "indel");
    assert_eq!(allele_shape(&mnv_row, &[]), "indel");
    assert_eq!(
        allele_shape(&mnv_row, &pairs(&[("A", "<DEL>"), ("AT", "GC")])),
        "mnv"
    );
    let del_row = motif_row("21:6001", "-");
    assert_eq!(allele_shape(&del_row, &pairs(&[("AT", "A")])), "indel");
    // A multi-allelic record with an indel loses the first base of every allele and nothing
    // more.
    let multi = pairs(&[("AT", "A"), ("AT", "ATT")]);
    assert_eq!(allele_shape(&motif_row("21:1001", "TT"), &multi), "indel");
    assert_eq!(allele_shape(&motif_row("21:1001", "-"), &multi), "indel");
    let cg_row = motif_row("21:1001-1002", "CG");
    assert_eq!(
        allele_shape(&cg_row, &pairs(&[("ATG", "AG"), ("ATG", "ACG")])),
        "mnv"
    );
    let g_row = motif_row("21:1000", "G");
    assert_eq!(
        allele_shape(&g_row, &pairs(&[("A", "AT"), ("A", "G")])),
        "snv"
    );
    let star_row = motif_row("21:1001", "*");
    assert_eq!(
        allele_shape(&star_row, &pairs(&[("AT", "A"), ("AT", "*")])),
        "indel"
    );
    // The motif allele shape reads a row its record's pairs do not match from the row alone.
    let mut trimmed = motif_row("21:401-402", "GG");
    trimmed.substitution = true;
    assert_eq!(trimmed.klass(), "mnv_or_complex");
    assert_eq!(allele_shape(&trimmed, &pairs(&[("AAAG", "AGGG")])), "indel");
    // Letters compare upper-cased; a row read alone is `snv`, `mnv` or `indel` by its class.
    assert_eq!(
        allele_shape(&motif_row("21:5000-5001", "gc"), &pairs(&[("at", "Gc")])),
        "mnv"
    );
    assert_eq!(allele_shape(&motif_row("21:7", "G"), &[]), "snv");
    assert_eq!(allele_shape(&motif_row("21:7", "GC"), &[]), "mnv");
    assert_eq!(allele_shape(&motif_row("21:7", "<DEL>"), &[]), "indel");
}

#[test]
fn regulatory_classes_follow_the_feature_types_of_a_record() {
    let row = |feature_type: &str, allele: &str, strand: &str| Row {
        feature_type: feature_type.to_string(),
        strand: strand.to_string(),
        ..motif_row("21:3000", allele)
    };
    assert_eq!(regulatory_classes(&[], &[]), BTreeSet::new());
    assert_eq!(
        regulatory_classes(&[row("Transcript", "G", "1")], &[]),
        BTreeSet::new()
    );
    let snv = pairs(&[("A", "G")]);
    assert_eq!(
        regulatory_classes(
            &[
                row("Transcript", "G", "1"),
                row("RegulatoryFeature", "G", "")
            ],
            &snv
        ),
        set(&["rf_only"])
    );
    assert_eq!(
        regulatory_classes(
            &[
                row("Transcript", "G", "1"),
                row("RegulatoryFeature", "G", ""),
                row("MotifFeature", "G", "-1"),
            ],
            &snv
        ),
        set(&["motif_reverse", "motif_snv", "rf_and_mf"])
    );
    assert_eq!(
        regulatory_classes(
            &[row("RegulatoryFeature", "G", ""), row("-", "G", "")],
            &snv
        ),
        set(&["regulatory_only", "rf_only"])
    );
    // A motif row with a symbolic allele shows no shape, and still its strand.
    assert_eq!(
        regulatory_classes(&[row("MotifFeature", "deletion", "-1")], &[]),
        set(&["mf_only", "motif_reverse", "regulatory_only"])
    );
    assert_eq!(
        regulatory_classes(&[row("MotifFeature", "-", "1")], &pairs(&[("AT", "A")])),
        set(&["mf_only", "motif_indel", "regulatory_only"])
    );
}
