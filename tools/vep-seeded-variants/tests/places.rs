// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Every record lies where its ID says: the place is worked out again from the cache
//! alone, with the geometry of `common`, for each pilot slice and for the made-up cache.

mod common;

use common::{
    as_read, bacterium, category_of, coding_index, codon_at, complete, exon_of, exons_low_to_high,
    introns, mitochondrion, mitochondrion_under, pilot_cache, records_of, spliced, substituted,
    Genome, Scratch, PILOTS,
};
use vep_core::codon::{is_start_codon, translate_codon_with_table};
use vep_seeded_variants::cache::{Strand, Transcript};
use vep_seeded_variants::generate::{Generated, Records};
use vep_seeded_variants::vcf::Record;

fn residue(codon: &[u8; 3], code: u8) -> u8 {
    translate_codon_with_table(codon, code)
}

/// Whether the positions from `low` to `high` lie in the coding part of one exon and, in
/// a transcript of several exons, six bases or more from that exon's ends.
fn inside_one_coding_exon(transcript: &Transcript, low: u64, high: u64) -> bool {
    let Some((start, end)) = exon_of(transcript, low) else {
        return false;
    };
    let clear = transcript.exons.len() == 1 || (start + 6 <= low && high + 6 <= end);
    high <= end && clear && (low..=high).all(|at| coding_index(transcript, at).is_some())
}

/// The lowest and the highest sequence position of the codon that holds `at`, where its
/// three bases are neighbours on the sequence.
fn codon_span(transcript: &Transcript, at: u64) -> Option<(u64, u64)> {
    let (number, offset, _) = codon_at(transcript, at)?;
    let offset = offset as u64;
    let (low, high) = match transcript.strand {
        Strand::Forward => (at.checked_sub(offset)?, at + 2 - offset),
        Strand::Reverse => ((at + offset).checked_sub(2)?, at + offset),
    };
    let same = |end: u64| codon_at(transcript, end).is_some_and(|(other, ..)| other == number);
    (same(low) && same(high)).then_some((low, high))
}

/// The check of one code category: `kind` is the category's name after its code.
fn code_record(transcript: &Transcript, code: u8, kind: &str, record: &Record) -> bool {
    if transcript.code() != code || !complete(transcript) {
        return false;
    }
    let Some((number, before, after)) = substituted(
        transcript,
        record.position,
        record.reference.as_bytes(),
        record.alt.as_bytes(),
    ) else {
        return false;
    };
    let codons = transcript
        .coding
        .as_ref()
        .map_or(0, |c| c.sequence.len() as u64 / 3);
    match kind {
        "start_codon" => number == 0,
        "stop_gained" => residue(&before, code) != b'*' && residue(&after, code) == b'*',
        "stop_lost" => {
            number + 1 == codons && residue(&before, code) == b'*' && residue(&after, code) != b'*'
        }
        "changed_codon" => {
            let differs = |codon: &[u8; 3]| is_start_codon(codon, code) != is_start_codon(codon, 1);
            residue(&before, code) != residue(&before, 1)
                || (number == 0 && (differs(&before) || differs(&after)))
        }
        "made_codon_synonymous" | "made_codon_missense" | "made_codon_stop" => {
            let unlike = |codon: &[u8; 3]| residue(codon, code) != residue(codon, 1);
            let (was, is) = (residue(&before, code), residue(&after, code));
            let made = !unlike(&before) && unlike(&after);
            // The two categories of a sense codon stand as the missense and the
            // synonymous records do: in a codon that is neither the first nor the
            // last, clear of its exon's ends.
            let inside = 0 < number
                && number + 1 < codons
                && codon_span(transcript, record.position)
                    .is_some_and(|(low, high)| inside_one_coding_exon(transcript, low, high));
            let whole = transcript
                .coding
                .as_ref()
                .is_some_and(|coding| coding.sequence.len() % 3 == 0);
            // The category of the closing stop codon stands as the stop-lost record
            // does: in a last codon whose three bases lie in one exon.
            let closing = whole
                && codons >= 2
                && number + 1 == codons
                && codon_span(transcript, record.position).is_some();
            match kind {
                "made_codon_synonymous" => made && inside && was == is && was != b'*',
                "made_codon_missense" => made && inside && was != is && was != b'*' && is != b'*',
                _ => made && closing && was == b'*',
            }
        }
        _ => false,
    }
}

/// Why a record does not lie where its ID says, or `Ok`.
fn lies_where_it_says(genome: &Genome, record: &Record) -> Result<(), String> {
    let category = category_of(&record.id);
    let sequence = genome
        .resolve(&record.sequence)
        .ok_or("the cache holds no such sequence")?;
    if category != "synonym" && sequence != record.sequence {
        return Err("the sequence is not written by the cache's own name".to_owned());
    }
    let transcripts = genome.on(&sequence);
    let at = record.position;
    let reference = record.reference.as_bytes();
    let alt = record.alt.as_bytes();
    let any = |test: &dyn Fn(&Transcript) -> bool| transcripts.iter().any(test);
    let inside_any = transcripts.iter().any(|t| t.start <= at && at <= t.end);
    let sequences = genome.cache.sequences();
    let holds = if let Some(rest) = category.strip_prefix("code") {
        let (code, kind) = rest.split_once('_').ok_or("no code in the ID")?;
        let code: u8 = code.parse().map_err(|_| "no code in the ID")?;
        any(&|t| code_record(t, code, kind, record))
    } else {
        match category {
            "first_sequence" => Some(&sequence) == sequences.first(),
            "last_sequence" => Some(&sequence) == sequences.last() && sequences.len() > 1,
            "organelle" => genome.organelles().contains(&sequence),
            "synonym" => {
                let unlinked = vep_core::synonyms::ChromosomeSynonyms::new();
                !genome.holds(&record.sequence)
                    && unlinked.source_chr_name(&record.sequence, |name| genome.holds(name))
                        != sequence
            }
            "missense_forward" | "synonymous_forward" => any(&|t| {
                t.strand == Strand::Forward && complete(t) && coding_index(t, at).is_some()
            }),
            "missense_reverse" | "synonymous_reverse" => any(&|t| {
                t.strand == Strand::Reverse && complete(t) && coding_index(t, at).is_some()
            }),
            "splice_donor" | "splice_acceptor" => any(&|t| {
                let at_low_end = (category == "splice_donor") == (t.strand == Strand::Forward);
                introns(t).iter().any(|(low, high)| {
                    high - low + 1 >= 20
                        && if at_low_end {
                            at == *low || at == low + 1
                        } else {
                            at == *high || at + 1 == *high
                        }
                })
            }),
            "five_prime_utr" => any(&|t| {
                complete(t)
                    && spliced(t, at)
                        .zip(t.coding.as_ref())
                        .is_some_and(|(place, coding)| place < coding.cdna_start)
            }),
            "three_prime_utr" => any(&|t| {
                complete(t)
                    && spliced(t, at)
                        .zip(t.coding.as_ref())
                        .is_some_and(|(place, coding)| place > coding.cdna_end)
            }),
            "intron" => any(&|t| {
                introns(t)
                    .iter()
                    .any(|(low, high)| low + 20 <= at && at + 20 <= *high)
            }),
            "upstream" | "downstream" => {
                !inside_any
                    && any(&|t| {
                        let before_start =
                            (category == "upstream") == (t.strand == Strand::Forward);
                        if before_start {
                            at < t.start && t.start - at <= 5000
                        } else {
                            at > t.end && at - t.end <= 5000
                        }
                    })
            }
            "intergenic" => transcripts
                .iter()
                .all(|t| at + 5000 < t.start || at > t.end + 5000),
            "frameshift_insertion" => {
                reference.len() == 1
                    && matches!(alt.len(), 2 | 3)
                    && alt.starts_with(reference)
                    && any(&|t| complete(t) && inside_one_coding_exon(t, at, at + 1))
            }
            "frameshift_deletion" => {
                matches!(reference.len(), 2 | 3)
                    && alt == &reference[..1]
                    && any(&|t| {
                        complete(t)
                            && inside_one_coding_exon(t, at, at + reference.len() as u64 - 1)
                    })
            }
            "inframe_deletion" => {
                reference.len() == 4
                    && alt == &reference[..1]
                    && any(&|t| {
                        let first = (at + 1..=at + 3).filter_map(|p| coding_index(t, p)).min();
                        complete(t)
                            && inside_one_coding_exon(t, at, at + 3)
                            && first.is_some_and(|first| first % 3 == 0)
                    })
            }
            "first_base" => at == 1,
            "last_base" => genome.length(&sequence) == Some(at),
            "origin" => {
                let region = genome.cache.sequence_region(&sequence);
                region
                    .is_some_and(|region| region.circular == Some(true) && at + 1 == region.length)
                    && reference.len() == 4
                    && alt == &reference[..1]
            }
            _ => return Err("the ID names no category".to_owned()),
        }
    };
    if holds {
        Ok(())
    } else {
        Err(format!("{}:{at} is not such a place", record.sequence))
    }
}

fn check_all(what: &str, genome: &Genome, generated: &Generated) {
    assert!(!generated.records.is_empty(), "{what}: no record");
    for record in &generated.records {
        if let Err(why) = lies_where_it_says(genome, record) {
            panic!("{what}: {}: {why}", record.id);
        }
    }
}

#[test]
fn each_record_lies_where_its_category_says() {
    for (corpus, _) in PILOTS {
        let dir = pilot_cache(corpus);
        let genome = Genome::load(&dir);
        for seed in [1, 2, 3, 20] {
            let light = records_of(&dir, seed, Records::Light);
            check_all(&format!("{corpus}, seed {seed}"), &genome, &light);
        }
        let wide = records_of(&dir, 1, Records::Wide(12));
        check_all(&format!("{corpus}, wide"), &genome, &wide);
    }
    for (what, dir) in [
        ("the made-up bacterium", bacterium()),
        ("the made-up cache under code 4", mitochondrion()),
    ] {
        let genome = Genome::load(&dir);
        for seed in 1..=10 {
            let light = records_of(&dir, seed, Records::Light);
            check_all(&format!("{what}, seed {seed}"), &genome, &light);
        }
    }
}

#[test]
fn each_pilot_slice_gets_its_codes() {
    for (corpus, code) in PILOTS {
        let generated = records_of(&pilot_cache(corpus), 1, Records::Light);
        let ids: Vec<&str> = generated.records.iter().map(|r| r.id.as_str()).collect();
        for kind in ["changed_codon", "start_codon", "stop_gained", "stop_lost"] {
            let id = format!("code{code}_{kind}_1");
            assert!(ids.contains(&id.as_str()), "{corpus}: no {id}");
        }
        for kind in ["start_codon", "stop_gained", "stop_lost"] {
            let id = format!("code1_{kind}_1");
            assert!(ids.contains(&id.as_str()), "{corpus}: no {id}");
        }
        assert!(
            !ids.iter().any(|id| id.starts_with("code1_changed_codon")),
            "{corpus}: the standard code has no codon it reads differently from itself"
        );
        let codes: Vec<&str> = ids
            .iter()
            .filter_map(|id| id.strip_prefix("code"))
            .filter_map(|rest| rest.split_once('_'))
            .map(|(code, _)| code)
            .collect();
        assert!(
            codes
                .iter()
                .all(|named| *named == "1" || *named == code.to_string()),
            "{corpus}: {codes:?}"
        );
    }
}

/// Whether a transcript of the genome on `code` holds a substitution of the category
/// `kind`, by a scan written here: every base of every exon against each of the three
/// bases it is not.
fn holds_a_substitution(genome: &Genome, code: u8, kind: &str) -> bool {
    genome.transcripts.iter().any(|(sequence, transcripts)| {
        transcripts
            .iter()
            .filter(|transcript| transcript.code() == code)
            .any(|transcript| {
                exons_low_to_high(transcript)
                    .into_iter()
                    .flat_map(|(start, end)| start..=end)
                    .any(|at| {
                        let Some((_, offset, codon)) = codon_at(transcript, at) else {
                            return false;
                        };
                        let reference = as_read(transcript, codon[offset]);
                        b"ACGT".iter().any(|alt| {
                            let record = Record {
                                sequence: sequence.clone(),
                                position: at,
                                id: String::new(),
                                reference: char::from(reference).to_string(),
                                alt: char::from(*alt).to_string(),
                            };
                            *alt != reference && code_record(transcript, code, kind, &record)
                        })
                    })
            })
    })
}

#[test]
fn a_made_codon_record_is_written_wherever_a_transcript_holds_its_substitution() {
    // Each cache with its code and, for the three categories in order, whether the
    // cache must hold a substitution; `None` leaves it to the scan.
    type Held = [Option<bool>; 3];
    const ALL: Held = [Some(true); 3];
    const NONE: Held = [Some(false); 3];
    let mut caches: Vec<(String, std::path::PathBuf, u8, Held)> = PILOTS
        .iter()
        .map(|(corpus, code)| {
            // The plastid's code differs from the standard code in start codons alone.
            let held = if *code == 11 { NONE } else { ALL };
            ((*corpus).to_owned(), pilot_cache(corpus), *code, held)
        })
        .collect();
    caches.push((
        "the made-up cache under code 4".to_owned(),
        mitochondrion(),
        4,
        ALL,
    ));
    // The same coding sequence under two codes that read one codon differently. Under
    // code 12 no substitution keeps a residue or starts from a stop codon; under code
    // 23 the codon is a stop codon, which the closing stop codon's category alone takes.
    let scratch = Scratch::new("made-other-codes");
    for (code, held) in [
        (12, [Some(false), None, Some(false)]),
        (23, [Some(false), Some(false), Some(true)]),
    ] {
        caches.push((
            format!("the made-up cache under code {code}"),
            mitochondrion_under(&scratch, code),
            code,
            held,
        ));
    }
    for (what, dir, code, expected) in caches {
        let genome = Genome::load(&dir);
        for (kind, expected) in [
            "made_codon_synonymous",
            "made_codon_missense",
            "made_codon_stop",
        ]
        .into_iter()
        .zip(expected)
        {
            let held = holds_a_substitution(&genome, code, kind);
            assert!(
                expected.is_none_or(|expected| expected == held),
                "{what}: {kind} is held: {held}"
            );
            // A code that translates every codon as the standard code does makes none.
            let translates_alike = (0..64u8).all(|number| {
                let codon =
                    [number >> 4, number >> 2, number].map(|two| b"ACGT"[(two & 3) as usize]);
                residue(&codon, code) == residue(&codon, 1)
            });
            assert!(!(translates_alike && held), "{what}: {kind}");
            for seed in [1, 2, 3, 20] {
                let light = records_of(&dir, seed, Records::Light);
                let id = format!("code{code}_{kind}_1");
                let written = light.records.iter().any(|record| record.id == id);
                assert_eq!(written, held, "{what}, seed {seed}: {id}");
            }
        }
    }
}

#[test]
fn the_organelle_records_follow_the_stated_rule() {
    let expected = [
        ("GRCm39-mouse", "MT"),
        ("BDGP6.54-fly", "mitochondrion_genome"),
        ("TAIR10-arabidopsis", "Pt"),
        ("R64-1-1-yeast", "Mito"),
    ];
    for (corpus, organelle) in expected {
        let dir = pilot_cache(corpus);
        assert_eq!(Genome::load(&dir).organelles(), [organelle], "{corpus}");
        let generated = records_of(&dir, 1, Records::Light);
        let on: Vec<&str> = generated
            .records
            .iter()
            .filter(|record| category_of(&record.id) == "organelle")
            .map(|record| record.sequence.as_str())
            .collect();
        assert_eq!(on, [organelle], "{corpus}");
    }
    // Every transcript of the made-up cache names one code and the cache holds two
    // sequences, so the plasmid, the shorter one, is the organelle.
    assert_eq!(Genome::load(&bacterium()).organelles(), ["pMADE"]);
    let generated = records_of(&bacterium(), 1, Records::Light);
    let on: Vec<&str> = generated
        .records
        .iter()
        .filter(|record| category_of(&record.id) == "organelle")
        .map(|record| record.sequence.as_str())
        .collect();
    assert_eq!(on, ["pMADE"]);
}

#[test]
fn missense_and_synonymous_are_what_they_claim() {
    let mut checked = 0;
    for (corpus, _) in PILOTS {
        let dir = pilot_cache(corpus);
        let genome = Genome::load(&dir);
        let mut records = records_of(&dir, 1, Records::Wide(12)).records;
        for seed in 1..=5 {
            records.extend(records_of(&dir, seed, Records::Light).records);
        }
        for record in &records {
            let category = category_of(&record.id);
            let (claim, strand) = match category {
                "missense_forward" => ("missense", Strand::Forward),
                "missense_reverse" => ("missense", Strand::Reverse),
                "synonymous_forward" => ("synonymous", Strand::Forward),
                "synonymous_reverse" => ("synonymous", Strand::Reverse),
                _ => continue,
            };
            let holds = genome.on(&record.sequence).iter().any(|transcript| {
                if transcript.strand != strand || !complete(transcript) {
                    return false;
                }
                let Some((_, before, after)) = substituted(
                    transcript,
                    record.position,
                    record.reference.as_bytes(),
                    record.alt.as_bytes(),
                ) else {
                    return false;
                };
                let code = transcript.code();
                let (was, is) = (residue(&before, code), residue(&after, code));
                let sense = was != b'*' && is != b'*';
                match claim {
                    "missense" => sense && was != is,
                    _ => sense && was == is && before != after,
                }
            });
            assert!(holds, "{corpus}: {} is not {claim}", record.id);
            checked += 1;
        }
    }
    assert!(checked >= 4 * 4 * 6, "only {checked} records were checked");
}

#[test]
fn a_codon_record_changes_one_base_of_its_codon() {
    for (corpus, code) in PILOTS {
        let dir = pilot_cache(corpus);
        let genome = Genome::load(&dir);
        for record in &records_of(&dir, 1, Records::Light).records {
            if !record.id.starts_with(&format!("code{code}_")) {
                continue;
            }
            let found = genome.on(&record.sequence).iter().any(|transcript| {
                transcript.code() == code
                    && codon_at(transcript, record.position).is_some()
                    && record.reference.len() == 1
                    && record.alt.len() == 1
                    && record.reference != record.alt
            });
            assert!(found, "{corpus}: {}", record.id);
        }
    }
}
