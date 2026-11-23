// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The class `allele_number_repeated_or_reference_equal_alt`: the allele number of a row of
//! a record whose written allele list repeats an alternate or lists an alternate equal to
//! the reference allele. The shape is read from the cells' texts alone;
//! `tools/vep-golden-corpus` compiles this file by path.

/// The fields that hold a row's allele number: `ALLELE_NUM` in the default, tab and VCF
/// formats, `allele_num` in the JSON format.
pub const ALLELE_NUMBER_FIELDS: [&str; 2] = ["ALLELE_NUM", "allele_num"];

/// The fields under which a row carries its record's written allele list:
/// `UPLOADED_ALLELE` in the default, tab and VCF formats, `uploaded_allele` in the JSON
/// format.
pub const WRITTEN_ALLELE_LIST_FIELDS: [&str; 2] = ["UPLOADED_ALLELE", "uploaded_allele"];

/// The written alternates as Ensembl VEP's VCF parser spells them on the rows of an indel
/// record (Parser/VCF.pm:294-335 at release/116.3), `None` for a list the parser leaves as
/// written. The parser chops the first character where both of these hold: an alternate's
/// length differs from the reference allele's, and the reference allele and every alternate
/// that holds no `*` begin with one character. The reference allele and each alternate that
/// holds no `*` are then the text after that character, `-` where nothing is left, and an
/// alternate that holds a `*` keeps its written text.
fn chopped_alternates<'a>(reference_allele: &str, alternates: &[&'a str]) -> Option<Vec<&'a str>> {
    let length = |allele: &str| allele.chars().count();
    let indel = alternates
        .iter()
        .any(|written| length(written) != length(reference_allele));
    let anchor = reference_allele.chars().next()?;
    let chopped_ones = || alternates.iter().filter(|written| !written.contains('*'));
    if !indel || !chopped_ones().all(|written| written.starts_with(anchor)) {
        return None;
    }
    let chopped = |written: &'a str| {
        if written.contains('*') {
            return written;
        }
        match &written[anchor.len_utf8()..] {
            "" => "-",
            rest => rest,
        }
    };
    Some(alternates.iter().copied().map(chopped).collect())
}

/// Whether the allele numbers of one row, `ev` the reference's and `av` the candidate's,
/// have the shape of `allele_number_repeated_or_reference_equal_alt`. `allele` is the row's
/// allele, and `reference_list` and `candidate_list` are the record's written allele list
/// as each side's row carries it, `None` for a row that carries none.
///
/// VariationFeatureOverlap.pm:459-467 at release/116.3 drop from the list an allele string
/// seen before, the reference allele's being the first seen, and :477-491 number what is
/// left: the reference allele's own entry is passed over and every other entry takes the
/// next number, `*` and `<DEL:*>` among them. The candidate numbers an allele by its place
/// in the written list.
///
/// The shape holds where all of this does, and nowhere else:
///
/// 1. Both rows carry the list and the two lists are one text.
/// 2. The list has three or more parts between its `/`: the reference allele, then the
///    written alternates in order. An empty part reads as `-`.
/// 3. `allele` is a written alternate, or an alternate as Ensembl VEP's VCF parser spells
///    it on the rows of an indel record (`chopped_alternates`). Its place is that of the
///    first such alternate among the written alternates, counted from 1. A row does not
///    say which of the two spellings its allele is: where the allele is a written
///    alternate at one place and the parser's spelling of an alternate at another, the
///    pair is not of the class. Nor is it where that alternate equals the reference
///    allele, an entry Ensembl VEP drops and writes no row for.
/// 4. Ensembl VEP's count for that place is the number of written alternates up to and
///    including it that equal neither the reference allele nor an alternate written before
///    them.
/// 5. `av` is the place and `ev` the count, each as its decimal text, and the two differ.
pub fn allele_number_place_pair(
    reference_list: Option<&str>,
    candidate_list: Option<&str>,
    allele: &str,
    ev: &str,
    av: &str,
) -> bool {
    let (Some(list), Some(candidate_list)) = (reference_list, candidate_list) else {
        return false;
    };
    if list != candidate_list {
        return false;
    }
    // Ensembl VEP reads an empty allele of the list as `-` where it drops the repeated
    // ones (VariationFeatureOverlap.pm:459-467).
    let parts: Vec<&str> = list
        .split('/')
        .map(|part| if part.is_empty() { "-" } else { part })
        .collect();
    if parts.len() < 3 {
        return false;
    }
    let (reference_allele, alternates) = (parts[0], &parts[1..]);
    // The place of the first alternate of `spelt` that is the allele, with Ensembl VEP's
    // count for that place, which is taken over the written alternates.
    let located = |spelt: &[&str]| -> Option<(usize, usize)> {
        let place = spelt.iter().position(|alternate| *alternate == allele)? + 1;
        let count = alternates[..place]
            .iter()
            .enumerate()
            .filter(|(at, written)| {
                **written != reference_allele && !alternates[..*at].contains(written)
            })
            .count();
        Some((place, count))
    };
    let written = located(alternates);
    let chopped =
        chopped_alternates(reference_allele, alternates).and_then(|chopped| located(&chopped));
    // The row spells its allele as its record writes it or as the parser chops it, and
    // does not say which: an allele both spellings place differently decides nothing.
    let (place, count) = match (written, chopped) {
        (Some(written), Some(chopped)) if written != chopped => return false,
        (Some(found), _) | (None, Some(found)) => found,
        (None, None) => return false,
    };
    // An alternate equal to the reference allele is dropped without a number
    // (VariationFeatureOverlap.pm:463) and writes no row.
    place != count
        && alternates[place - 1] != reference_allele
        && av == place.to_string()
        && ev == count.to_string()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The rule over one list that both rows carry.
    fn pair(list: &str, allele: &str, ev: &str, av: &str) -> bool {
        allele_number_place_pair(Some(list), Some(list), allele, ev, av)
    }

    #[test]
    fn an_alternate_written_after_a_repeated_one_is_taken() {
        assert!(pair("C/G/G/T", "T", "2", "3"));
        // Each repeat before the allele lowers Ensembl VEP's number by one more.
        assert!(pair("C/G/G/G/T", "T", "2", "4"));
        assert!(pair("C/G/T/G/T/A", "A", "3", "5"));
    }

    #[test]
    fn an_alternate_written_after_one_equal_to_the_reference_allele_is_taken() {
        assert!(pair("A/A/C", "C", "1", "2"));
        assert!(pair("A/C/A/G", "G", "2", "3"));
        assert!(pair("A/A/C/C/G", "G", "2", "4"));
    }

    #[test]
    fn a_star_and_a_symbolic_deletion_take_a_number_as_any_alternate_does() {
        assert!(pair("C/*/G/G/T", "T", "3", "4"));
        assert!(pair("C/<DEL:*>/C/T", "T", "2", "3"));
        assert!(pair("C/*/*/G", "G", "2", "3"));
    }

    #[test]
    fn an_empty_part_of_the_list_reads_as_a_dash() {
        assert!(pair("A//G/G/T", "T", "3", "4"));
        assert!(pair("A/A//", "-", "1", "2"));
        assert!(pair("/-/T", "T", "1", "2"));
    }

    #[test]
    fn a_row_without_the_list_or_with_another_list_stays_a_difference() {
        let list = Some("C/G/G/T");
        assert!(allele_number_place_pair(list, list, "T", "2", "3"));
        assert!(!allele_number_place_pair(None, list, "T", "2", "3"));
        assert!(!allele_number_place_pair(list, None, "T", "2", "3"));
        assert!(!allele_number_place_pair(None, None, "T", "2", "3"));
        assert!(!allele_number_place_pair(
            list,
            Some("C/G/T"),
            "T",
            "2",
            "3"
        ));
        // A list of fewer than three parts, as a row without the field reads in the tab
        // format, holds no entry to drop.
        for short in ["-", "", "C", "C/T", "T/T"] {
            assert!(!pair(short, "T", "2", "3"), "{short}");
            assert!(!pair(short, "T", "0", "1"), "{short}");
        }
    }

    #[test]
    fn an_allele_the_list_does_not_hold_stays_a_difference() {
        assert!(!pair("C/G/G/T", "A", "2", "3"));
        // The reference allele's own entry is no written alternate.
        assert!(!pair("C/G/G/T", "C", "0", "1"));
    }

    #[test]
    fn a_row_of_an_indel_record_is_matched_by_the_allele_the_parser_spells() {
        // The parser spells `CA/C/C/CAT` as `A/-/-/AT` on the record's rows.
        assert!(pair("CA/C/C/CAT", "AT", "2", "3"));
        // The repeated deletion has its row under its first occurrence, where the place and
        // the count are one number.
        assert!(!pair("CA/C/C/CAT", "-", "1", "2"));
        assert!(!pair("CA/C/C/CAT", "-", "1", "1"));
        // An alternate equal to the reference allele is dropped under either spelling.
        assert!(pair("TA/TA/T", "-", "1", "2"));
        assert!(pair("C/C/CG", "G", "1", "2"));
        // An alternate that holds a `*` keeps its text, and its first character is not
        // held to the others'.
        assert!(pair("CA/*/C/C/CAT", "AT", "3", "4"));
        assert!(pair("CA/<DEL:*>/CA/CAT", "AT", "2", "3"));
        assert!(pair("CA/C/C/*", "*", "2", "3"));
    }

    #[test]
    fn an_allele_of_one_spelling_alone_is_placed_by_that_spelling() {
        // The alternates do not share the reference allele's first base, so the parser
        // leaves the list as written and a row carries the written allele.
        assert!(pair("CA/C/C/GAT", "GAT", "2", "3"));
        assert!(!pair("CA/C/C/GAT", "AT", "2", "3"));
        // An alternate that holds a `*` is spelt alike as written and as chopped: both
        // spellings give one place and one count, which decide.
        assert!(pair("CA/C/C/*", "*", "2", "3"));
    }

    #[test]
    fn an_allele_the_two_spellings_place_differently_stays_a_difference() {
        // `GG/G/GGG/GG` is chopped to `-`, `GG` and `G`. The allele `GG` is the third
        // written alternate, where the count is 2, and the chopped second, where the count
        // is 2 as well: a row of the record prints 2 in both engines, and 2 and 3 are the
        // count and the place of the written spelling. The row does not say which it is.
        assert!(!pair("GG/G/GGG/GG", "GG", "2", "3"));
        // `C/CT/CT/CCT` is chopped to `T`, `T` and `CT`. The allele `CT` is the first
        // written alternate, numbered 1 by both engines, and the chopped third, numbered 2
        // and 3.
        assert!(!pair("C/CT/CT/CCT", "CT", "2", "3"));
        // The allele `G` of the first list is its first written alternate and its chopped
        // third, where the count is 2.
        assert!(!pair("GG/G/GGG/GG", "G", "2", "3"));
        // `GG/G/GGG/GG/GGGG` is chopped to `-`, `GG`, `G` and `GGG`: the allele `GGG` is
        // the second written alternate and the chopped fourth, where the count is 3.
        assert!(!pair("GG/G/GGG/GG/GGGG", "GGG", "3", "4"));
        // An allele of the same lists that one spelling alone holds is placed by it: the
        // written `GGGG` and `CCT`, and the chopped `T`, the repeated alternate's first
        // occurrence, where the place and the count are one number.
        assert!(pair("GG/G/GGG/GG/GGGG", "GGGG", "3", "4"));
        assert!(pair("C/CT/CT/CCT", "CCT", "2", "3"));
        assert!(!pair("C/CT/CT/CCT", "T", "1", "2"));
    }

    #[test]
    fn an_allele_of_neither_spelling_stays_a_difference() {
        assert!(!pair("CA/C/C/CAT", "T", "2", "3"));
        // The parser's spelling of the reference allele is no alternate.
        assert!(!pair("CA/C/C/CAT", "A", "2", "3"));
        // Alleles of one length make no indel record, and nothing is chopped.
        assert!(!pair("CA/CG/CG/CT", "T", "2", "3"));
        // Nor is anything chopped where an alternate begins with another character, a base
        // in the other case included.
        assert!(!pair("CA/C/C/cAT", "AT", "2", "3"));
        assert!(!pair("CA/G/G/CAT", "AT", "2", "3"));
    }

    #[test]
    fn an_alternate_equal_to_the_reference_allele_stays_a_difference() {
        // VariationFeatureOverlap.pm:463 drops an alternate equal to the reference allele,
        // so no row carries it: its place and the count of the alternates kept before it
        // are no pair of the class.
        assert!(!pair("C/G/C/T", "C", "1", "2"));
        assert!(!pair("A/A/C", "A", "0", "1"));
        assert!(!pair("A/*/A", "A", "1", "2"));
        // The parser spells `CA/C/CA/CAT` as `A/-/A/AT`: the second alternate is the
        // reference allele under either spelling.
        assert!(!pair("CA/C/CA/CAT", "CA", "1", "2"));
        assert!(!pair("CA/C/CA/CAT", "A", "1", "2"));
        // The alternates written after it keep their pairs.
        assert!(pair("C/G/C/T", "T", "2", "3"));
        assert!(pair("A/A/C", "C", "1", "2"));
        assert!(pair("CA/C/CA/CAT", "AT", "2", "3"));
    }

    #[test]
    fn numbers_that_are_not_the_place_and_the_count_stay_a_difference() {
        // The candidate numbered 3 where the allele's place is 2.
        assert!(!pair("A/A/C", "C", "1", "3"));
        // The two sides exchanged, and a reference number that is not the count.
        assert!(!pair("C/G/G/T", "T", "3", "2"));
        assert!(!pair("C/G/G/T", "T", "1", "3"));
        // A repeated alternate has rows under its first occurrence alone, where the place
        // and the count are one number.
        assert!(!pair("C/G/G/T", "G", "1", "2"));
        assert!(!pair("C/G/G/T", "G", "2", "2"));
        // A list that drops no entry before the allele numbers it alike on both sides.
        assert!(!pair("C/G/T/T", "T", "2", "3"));
        assert!(!pair("C/G/A/T", "T", "2", "3"));
        // A number in any other text than its decimal one, and no number.
        for (ev, av) in [
            ("02", "3"),
            ("2", "03"),
            ("2", "+3"),
            ("2.0", "3"),
            ("-", "3"),
        ] {
            assert!(!pair("C/G/G/T", "T", ev, av), "{ev} {av}");
        }
        assert!(!pair("C/G/G/T", "T", "2", "-"));
    }
}
