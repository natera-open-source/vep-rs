// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The ordered map, the reservoir step, the pass that samples the rows of a reference on
//! any count of threads, and the choice of exemplars among sampled rows.

use super::*;
use crate::TempDir;
use std::fs;

fn scratch() -> TempDir {
    TempDir::create(&std::env::temp_dir(), "vep-golden-corpus-sampling-").unwrap()
}

fn row(name: &str, location: &str, allele: &str, consequence: &str, extra: &str) -> Row {
    Row {
        dataset: "s".to_string(),
        uploaded_variation: name.to_string(),
        location: location.to_string(),
        allele: allele.to_string(),
        feature: "ENST1".to_string(),
        feature_type: "Transcript".to_string(),
        consequence_set: consequence.to_string(),
        strand: extra.to_string(),
        flags: String::new(),
        positions: ["-", "-", "-"].map(str::to_string),
        family: Vec::new(),
        substitution: false,
    }
}

fn names(rows: &[Row]) -> Vec<&str> {
    rows.iter()
        .map(|row| row.uploaded_variation.as_str())
        .collect()
}

fn reservoir_of(combinations: &[(&str, Vec<Row>)]) -> OrderedMap<String, Sample> {
    let mut reservoir = OrderedMap::new();
    for (combination, rows) in combinations {
        let sample = Sample {
            seen: rows.len(),
            rows: rows.clone(),
        };
        reservoir.insert(combination.to_string(), sample);
    }
    reservoir
}

#[test]
fn an_ordered_map_keeps_its_keys_where_they_were_first_set() {
    let mut map: OrderedMap<String, Vec<u32>> = OrderedMap::new();
    assert!(map.is_empty());
    map.entry("b").push(1);
    map.entry("a").push(2);
    map.entry("b").push(3);
    map.insert("c".to_string(), vec![4]);
    map.insert("a".to_string(), vec![5]);
    let entries: Vec<(&str, &[u32])> = map
        .iter()
        .map(|(key, value)| (key.as_str(), value.as_slice()))
        .collect();
    assert_eq!(
        entries,
        [("b", &[1, 3][..]), ("a", &[5][..]), ("c", &[4][..])]
    );
    assert_eq!(
        (map.len(), map.place("a"), map.place("z")),
        (3, Some(1), None)
    );
    assert_eq!(map.place_or_insert_with("c", Vec::new), 2);
    assert_eq!(map.place_or_insert_with("d", || vec![6]), 3);
    assert_eq!(map.at(3), (&"d".to_string(), &vec![6]));
    map.retain_places(|place| place != 1);
    assert_eq!(map.keys().collect::<Vec<_>>(), ["b", "c", "d"]);
    assert_eq!((map.place("d"), map.get("a")), (Some(2), None));
    assert!(map.contains_key("c") && !map.contains_key("a"));
}

#[test]
fn a_reservoir_step_fills_the_bucket_then_draws_a_place() {
    let mut rng = Random::new(1);
    let mut bucket = Vec::new();
    for offered in 1..=3 {
        reservoir_step(&mut rng, &mut bucket, offered, offered, 3);
    }
    assert_eq!(bucket, [1, 2, 3]);
    // No draw is made while the bucket fills.
    assert_eq!(rng.getrandbits(32), Random::new(1).getrandbits(32));
    // A capacity of zero or less makes the draw and keeps nothing.
    let (mut drawing, mut plain) = (Random::new(5), Random::new(5));
    let mut empty: Vec<u32> = Vec::new();
    for capacity in [0, -1] {
        reservoir_step_with(&mut drawing, &mut empty, 10, capacity, || 7);
        plain.randrange(10);
    }
    assert!(empty.is_empty());
    assert_eq!(drawing.getrandbits(32), plain.getrandbits(32));
    // The item is asked for only when it is kept.
    let mut asked = 0;
    let mut bucket = vec![0u32];
    let mut rng = Random::new(9);
    for offered in 2..200 {
        reservoir_step_with(&mut rng, &mut bucket, offered, 1, || {
            asked += 1;
            offered as u32
        });
    }
    assert!(asked > 0 && asked < 20, "{asked}");
}

/// A reference of `rows` rows over a few consequence sets and Extra values, some lines of
/// it no rows, and an input whose records hold a two-base substitution every tenth record.
fn reference(directory: &TempDir, rows: usize, line_end: &str) -> (String, String) {
    let sets = [
        "missense_variant",
        "intron_variant,non_coding_transcript_variant",
        "non_coding_transcript_variant,intron_variant",
        "synonymous_variant",
        "upstream_gene_variant",
        " intron_variant ",
    ];
    let mut text = format!("## a header{line_end}#Uploaded_variation\tLocation{line_end}");
    let mut input =
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n".to_string();
    let mut state = 7u32;
    for place in 0..rows {
        state = state.wrapping_mul(1_103_515_245).wrapping_add(12345);
        let set = sets[(state >> 16) as usize % sets.len()];
        let name = format!("rs{}", place / 3);
        let position = 1000 + place / 3 * 10;
        let (location, allele) = if place / 3 % 10 == 0 {
            (format!("21:{position}-{}", position + 1), "TG")
        } else {
            (format!("21:{position}"), "G")
        };
        let impact = ["HIGH", "LOW", "MODIFIER"][(state >> 8) as usize % 3];
        text += &format!(
            "{name}\t{location}\t{allele}\tENSG1\tENST{place}\tTranscript\t{set}\t-\t-\t-\t-\t-\t-\tIMPACT={impact};STRAND=1{line_end}"
        );
        if place % 3 == 0 {
            let (ref_allele, alt) = if place / 3 % 10 == 0 {
                ("AC", "TG")
            } else {
                ("A", "G")
            };
            input += &format!("21\t{position}\t{name}\t{ref_allele}\t{alt}\t.\t.\t.\n");
        }
        if place % 97 == 0 {
            text += &format!("too\tfew\tcolumns{line_end}{line_end}");
        }
    }
    let reference = directory.path.join("ref.txt");
    fs::write(&reference, text).unwrap();
    let input_path = directory.path.join("in.vcf");
    fs::write(&input_path, input).unwrap();
    (
        reference.to_str().unwrap().to_string(),
        input_path.to_str().unwrap().to_string(),
    )
}

/// What a pass holds, as text.
fn pass_text(pass: &ReservoirPass) -> String {
    let sample = |sample: &Sample| format!("{} {:?}", sample.seen, sample.rows);
    let mut text = String::new();
    for (set, held) in pass.reservoir.iter() {
        text += &format!("set {set}: {}\n", sample(held));
    }
    for (stratum, held) in pass.families.iter() {
        text += &format!("stratum {stratum}: {}\n", sample(held));
    }
    for (key, rows) in pass.named_rows.iter() {
        text += &format!("named {key:?}: {rows:?}\n");
    }
    let mut datasets: Vec<String> = pass
        .per_dataset
        .iter()
        .map(|(sid, count)| format!("{sid} {} {}", count.rows, count.combinations.len()))
        .collect();
    datasets.sort();
    text + &datasets.join("\n")
}

#[test]
fn a_pass_samples_the_same_rows_on_every_count_of_threads() {
    let scratch = scratch();
    // More rows than one block of the reader holds, so that blocks are read side by side.
    let (reference_path, input) = reference(&scratch, 60_000, "\n");
    let datasets = vec![
        ("a".to_string(), reference_path.clone(), input.clone()),
        ("b".to_string(), reference_path, input),
    ];
    let named: HashMap<String, HashSet<String>> = HashMap::from([(
        "b".to_string(),
        HashSet::from([
            "rs10".to_string(),
            "rs19999".to_string(),
            "absent".to_string(),
        ]),
    )]);
    let specs = ["IMPACT=".to_string()];
    let one = reservoir_pass(&datasets, 3, &specs, &named, 1).unwrap();
    let one_text = pass_text(&one);
    for threads in [2, 5, 64] {
        let pass = reservoir_pass(&datasets, 3, &specs, &named, threads).unwrap();
        assert!(pass_text(&pass) == one_text, "{threads} threads");
    }
    // Six ways to write five sets: the third is the second with its terms in another order.
    let sets: Vec<&String> = one.reservoir.keys().collect();
    assert_eq!(sets.len(), 5);
    assert!(one
        .reservoir
        .values()
        .all(|sample| sample.rows.len() == RESERVOIR));
    assert_eq!(
        one.reservoir
            .values()
            .map(|sample| sample.seen)
            .sum::<usize>(),
        120_000
    );
    assert_eq!(one.families.len(), 6);
    assert!(one
        .families
        .values()
        .all(|sample| sample.rows.len() == FAMILY_RESERVOIR));
    assert_eq!(
        (one.per_dataset["a"].rows, one.per_dataset["b"].rows),
        (60_000, 60_000)
    );
    assert_eq!(one.per_dataset["a"].combinations.len(), 5);
    let keys: Vec<&(String, String)> = one.named_rows.keys().collect();
    assert_eq!(
        keys,
        [
            &("b".to_string(), "rs10".to_string()),
            &("b".to_string(), "rs19999".to_string())
        ]
    );
    assert!(one.named_rows.values().all(|rows| rows.len() == 3));
    // A row of a two-base substitution is read with the input, in a sample and in a
    // family stratum.
    let substituted = one
        .named_rows
        .get(&("b".to_string(), "rs10".to_string()))
        .unwrap();
    assert!(substituted
        .iter()
        .all(|row| row.substitution && row.klass() == "mnv_or_complex"));
    assert!(one
        .families
        .keys()
        .any(|stratum| stratum.ends_with(";mnv_or_complex")));
    // Another seed draws other rows; without family specs no stratum is sampled.
    let other = reservoir_pass(&datasets, 4, &[], &HashMap::new(), 3).unwrap();
    assert!(other.families.is_empty() && other.named_rows.is_empty());
    assert_ne!(
        format!(
            "{:?}",
            other.reservoir.values().next().unwrap().rows[..8].to_vec()
        ),
        format!(
            "{:?}",
            one.reservoir.values().next().unwrap().rows[..8].to_vec()
        )
    );
}

#[test]
fn a_pass_refuses_an_input_before_its_reference_and_a_dataset_after_the_ones_before_it() {
    let scratch = scratch();
    let (reference_path, input) = reference(&scratch, 30, "\r\n");
    let absent = scratch.path.join("absent");
    let absent = absent.to_str().unwrap();
    let dataset =
        |reference: &str, input: &str| ("s".to_string(), reference.to_string(), input.to_string());
    let refusal = |datasets: &[Dataset]| {
        let pass = reservoir_pass(datasets, 1, &[], &HashMap::new(), 2);
        pass.unwrap_err().line
    };
    let unreadable = |path: &str| {
        format!(
            "ERROR: [build_golden_corpus] {path} is not a readable text file ([Errno 2] No such file or directory: '{path}')"
        )
    };
    assert_eq!(refusal(&[dataset(absent, absent)]), unreadable(absent));
    let other = format!("{absent}.txt");
    assert_eq!(refusal(&[dataset(&other, &input)]), unreadable(&other));
    assert_eq!(
        refusal(&[dataset(&reference_path, &input), dataset(&other, absent)]),
        unreadable(absent)
    );
    // Rows are read whatever the line end.
    let pass = reservoir_pass(
        &[dataset(&reference_path, &input)],
        1,
        &[],
        &HashMap::new(),
        2,
    );
    assert_eq!(pass.unwrap().per_dataset["s"].rows, 30);
}

#[test]
fn exemplars_are_the_first_rows_of_distinct_strata_in_their_order() {
    let wide = row("wide", "21:100-5000", "deletion", "c", "1");
    let flagged = Row {
        flags: "cds_end_NF".to_string(),
        ..row("flagged", "21:300", "G", "c", "1")
    };
    let reverse = row("reverse", "21:400", "G", "c", "-1");
    let plain = row("plain", "21:500", "G", "c", "1");
    let rows = vec![
        wide.clone(),
        plain.clone(),
        reverse.clone(),
        flagged.clone(),
    ];
    let reservoir = reservoir_of(&[("c", rows.clone()), ("a,focus", rows.clone())]);
    let choose = |k: i64, max_span: i64, focus: &[&str], focus_k: Option<i64>| {
        let focus: Vec<String> = focus.iter().map(|term| term.to_string()).collect();
        let mut rng = Random::new(1);
        choose_exemplars(&reservoir, k, &mut rng, max_span, &focus, focus_k, None).unwrap()
    };
    // The combinations in sorted order; the highest preference first, a row wider than the
    // bound last, and a second stratum before a second row of the first.
    let chosen = choose(2, 1000, &[], None);
    assert_eq!(chosen.keys().collect::<Vec<_>>(), ["a,focus", "c"]);
    assert_eq!(names(chosen.get("c").unwrap()), ["flagged", "reverse"]);
    let three = choose(3, 1000, &[], None);
    assert_eq!(names(&three.get("c").unwrap()[..2]), ["flagged", "reverse"]);
    assert_eq!(names(choose(4, 1000, &[], None).get("c").unwrap()).len(), 4);
    assert_eq!(
        names(choose(4, 1000, &[], None).get("c").unwrap())[2],
        "wide"
    );
    // Without a bound the wide row is one more stratum among the others.
    assert_eq!(
        names(choose(1, 0, &[], None).get("c").unwrap()),
        ["flagged"]
    );
    // A count of zero or less takes every candidate; a focus term takes its own count.
    for k in [0, -3] {
        assert_eq!(choose(k, 1000, &[], None).get("c").unwrap().len(), 4);
    }
    let focused = choose(1, 1000, &["focus", "other"], Some(3));
    assert_eq!(focused.get("a,focus").unwrap().len(), 3);
    assert_eq!(focused.get("c").unwrap().len(), 1);
    assert_eq!(
        choose(1, 1000, &["focus"], None)
            .get("a,focus")
            .unwrap()
            .len(),
        1
    );
    // A row the reservoir holds twice is chosen once.
    let twice = reservoir_of(&[("c", vec![plain.clone(), plain.clone()])]);
    let chosen = choose_exemplars(&twice, 5, &mut Random::new(1), 0, &[], None, None).unwrap();
    assert_eq!(names(chosen.get("c").unwrap()), ["plain"]);
    assert_eq!(exemplar_order(&wide, 1000, None), Ok((1, 0, 0)));
    assert_eq!(exemplar_order(&wide, 0, None), Ok((0, 0, 0)));
    assert_eq!(exemplar_order(&plain, -1, None), Ok((1, 0, 0)));
    assert_eq!(exemplar_order(&flagged, 1000, None), Ok((0, 0, -4)));
}

#[test]
fn every_class_and_every_family_stratum_keeps_an_exemplar() {
    let snv = row("snv", "21:100", "G", "b", "1");
    let deletion = row("deletion", "21:200-205", "-", "b", "1");
    let later_deletion = Row {
        flags: "cds_start_NF".to_string(),
        ..row("flagged_deletion", "21:300-305", "-", "a", "1")
    };
    let insertion = row("insertion", "21:400-401", "TT", "a", "1");
    let reservoir = reservoir_of(&[
        ("b", vec![snv.clone(), deletion.clone()]),
        ("a", vec![insertion.clone(), later_deletion.clone()]),
    ]);
    let mut chosen = Chosen::new();
    chosen.insert("a".to_string(), vec![insertion.clone()]);
    chosen.insert("b".to_string(), vec![]);
    ensure_class_exemplars(&mut chosen, &reservoir, 0, None).unwrap();
    // The classes no chosen row has, in the order the reservoir first shows them: the snv,
    // then the deletion that sorts first, each in the combination it was sampled for.
    assert_eq!(names(chosen.get("b").unwrap()), ["snv"]);
    assert_eq!(
        names(chosen.get("a").unwrap()),
        ["insertion", "flagged_deletion"]
    );

    let family = |name: &str, classes: &[&str], consequence: &str| Row {
        family: classes.iter().map(|class| class.to_string()).collect(),
        ..row(name, "21:100", "G", consequence, "1")
    };
    let mut families: OrderedMap<String, Sample> = OrderedMap::new();
    let high = vec![
        family("high1", &["IMPACT=HIGH"], "b"),
        family("high2", &["IMPACT=HIGH"], "b"),
    ];
    let low = vec![family("low1", &["IMPACT=LOW"], "new")];
    let held = vec![family("held", &["IMPACT=MODIFIER"], "a")];
    for (stratum, rows) in [
        ("IMPACT=LOW;snv", low),
        ("IMPACT=HIGH;snv", high),
        ("IMPACT=MODIFIER;snv", held.clone()),
    ] {
        families.insert(
            stratum.to_string(),
            Sample {
                seen: rows.len(),
                rows,
            },
        );
    }
    let mut chosen = Chosen::new();
    chosen.insert("a".to_string(), held);
    ensure_family_exemplars(&mut chosen, &families, 0, None).unwrap();
    // The strata in sorted order, each the first row of its sample that sorts first; a
    // combination no row was chosen for goes behind the others.
    assert_eq!(chosen.keys().collect::<Vec<_>>(), ["a", "b", "new"]);
    assert_eq!(names(chosen.get("a").unwrap()), ["held"]);
    assert_eq!(names(chosen.get("b").unwrap()), ["high1"]);
    assert_eq!(names(chosen.get("new").unwrap()), ["low1"]);
}
