// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! How a comparison that stops leaves the report directory. Every stop sets the discordance
//! report aside under its `.partial` name. Only the refusal of keys that do not ascend keeps
//! the directory of the key files, and its message names the output file and the lines of it
//! that first write the key it stopped on and the key before it, and neither key; any other
//! error of the merge, carried as the same error type, removes the directory.

use super::*;

/// The names of the entries of `dir`, sorted.
fn entries(dir: &Path) -> Vec<String> {
    let mut names: Vec<String> = std::fs::read_dir(dir)
        .unwrap()
        .map(|entry| entry.unwrap().file_name().to_string_lossy().into_owned())
        .collect();
    names.sort_unstable();
    names
}

/// A comparison under way: the report directory, its discordance report open with the header
/// written, and the directory of the key files with one key file in it.
struct UnderWay {
    report: TempDir,
    discordant_path: PathBuf,
    discordant_out: BufWriter<File>,
    key_files: TempDir,
}

fn under_way() -> UnderWay {
    let report = TempDir::create(&std::env::temp_dir(), "compare_stop_test_").unwrap();
    let discordant_path = report.path.join("discordant.tsv");
    let mut discordant_out = BufWriter::new(File::create(&discordant_path).unwrap());
    discordant_out
        .write_all(b"file_name\tsource\tlocation\tallele\tfeature\tfeature_type\tconsequence_set\n")
        .unwrap();
    let key_files = TempDir::create(&report.path, "compare_vep_tmp_").unwrap();
    std::fs::write(
        key_files.path.join("output.txt.perl.keys.sorted.txt"),
        "b\na\n",
    )
    .unwrap();
    UnderWay {
        report,
        discordant_path,
        discordant_out,
        key_files,
    }
}

/// The refusal of the second of the reference's `keys`, which does not sort after the first.
fn refusal_of_key_order(keys: [&str; 2]) -> io::Error {
    let keys = SortedKeys::in_the_order_written(format!("{}\n{}\n", keys[0], keys[1]).into());
    let none = SortedKeys::in_the_order_written(Vec::new());
    merge_keys(
        &mut MemoryKeys::new(&keys),
        &mut MemoryKeys::new(&none),
        "output.txt",
        None::<&mut Vec<u8>>,
        None,
        false,
        Collector::None,
    )
    .map(drop)
    .unwrap_err()
}

const NEEDS: &str = "; the comparison needs each side's keys in order with every key once";

#[test]
fn an_error_other_than_keys_out_of_order_removes_the_key_files_and_sets_the_report_aside() {
    let run = under_way();
    let absent = run.report.path.join("absent.txt");
    let partial = run.report.path.join("discordant.tsv.partial");
    let stopped = Stopped::on(
        "output.txt",
        io::Error::other("another error"),
        &run.key_files,
        [&absent, &absent],
        None,
    );
    drop(run.key_files);
    assert_eq!(
        stopped.reported(run.discordant_out, &run.discordant_path),
        format!(
            "output.txt: another error; the rows written before the stop are kept in {}",
            partial.display()
        )
    );
    assert_eq!(entries(&run.report.path), ["discordant.tsv.partial"]);
    assert!(std::fs::read_to_string(&partial)
        .unwrap()
        .starts_with("file_name\tsource\t"));

    // The same steps under the refusal keep the key files, so their removal above is the
    // error's doing. No file yields the two keys, so the message names the side alone.
    let run = under_way();
    let partial = run.report.path.join("discordant.tsv.partial");
    let stopped = Stopped::on(
        "output.txt",
        refusal_of_key_order(["b", "a"]),
        &run.key_files,
        [&absent, &absent],
        None,
    );
    let kept = run.key_files.path.clone();
    drop(run.key_files);
    assert_eq!(
        stopped.reported(run.discordant_out, &run.discordant_path),
        format!(
            "output.txt: the sorted keys of the --perl-dir side do not ascend: a key of the \
             side does not sort after the key before it{NEEDS}; the key files are kept in {} \
             and the rows written before the stop in {}",
            kept.display(),
            partial.display()
        )
    );
    assert_eq!(entries(&run.report.path).len(), 2);
    assert!(kept.join("output.txt.perl.keys.sorted.txt").is_file());
    assert!(partial.is_file());
}

/// A stop that names no paired file, as the run makes of a failure outside the comparison
/// of one, sets the report aside as well.
#[test]
fn a_stop_made_of_a_message_sets_the_report_aside() {
    let run = under_way();
    let stopped = Stopped::from("the report cannot be written".to_string());
    assert_eq!(
        stopped.reported(run.discordant_out, &run.discordant_path),
        format!(
            "the report cannot be written; the rows written before the stop are kept in {}",
            run.report.path.join("discordant.tsv.partial").display()
        )
    );
    let kept = run.key_files.path.clone();
    drop(run.key_files);
    assert!(!kept.exists());
    assert_eq!(entries(&run.report.path), ["discordant.tsv.partial"]);
}

/// A report that cannot be renamed stays under its name, which the message gives with the
/// reason.
#[test]
fn a_report_that_cannot_be_renamed_is_named_as_it_stands() {
    let run = under_way();
    std::fs::create_dir(run.report.path.join("discordant.tsv.partial")).unwrap();
    std::fs::write(
        run.report.path.join("discordant.tsv.partial").join("held"),
        "",
    )
    .unwrap();
    let stopped = Stopped::from("another error".to_string());
    let message = stopped.reported(run.discordant_out, &run.discordant_path);
    let opening = format!(
        "another error; the rows written before the stop are kept in {}, which was not renamed: ",
        run.discordant_path.display()
    );
    assert!(message.starts_with(&opening), "{message}");
    assert!(run.discordant_path.is_file());
}

/// Three rows of an output, the last two of one tuple key, behind two `#` lines, an empty
/// line and a short row, the lines ended in three ways, and the two keys the rows write.
const OUTPUT: &str = "## header\r\n#Uploaded_variation\tLocation\n\nshort\trow\r\
                      v1\t1:200\tk\tg1\tt1\tTranscript\tintron_variant,NMD_transcript_variant\n\
                      v2\t1:100\tk\tg1\tt1\tTranscript\tintron_variant\n\
                      v3\t1:100\tk\tg1\tt1\tTranscript\tintron_variant\n";
const KEY_OF_LINE_5: &str = "1:200\tk\tt1\tTranscript\tNMD_transcript_variant,intron_variant";
const KEY_OF_LINE_6: &str = "1:100\tk\tt1\tTranscript\tintron_variant";

#[test]
fn the_line_of_a_key_is_counted_over_every_line_of_the_file() {
    let scratch = TempDir::create(&std::env::temp_dir(), "compare_stop_test_").unwrap();
    let output = scratch.path.join("output.txt");
    std::fs::write(&output, OUTPUT).unwrap();
    let lines = |keys: [&str; 2], collapsed: bool, contigs: Option<&HashSet<String>>| {
        first_lines_writing(&output, keys.map(str::as_bytes), collapsed, contigs).unwrap()
    };
    assert_eq!(
        lines([KEY_OF_LINE_6, KEY_OF_LINE_5], false, None),
        [Some(6), Some(5)],
        "the first of the two lines that write the key of line 6"
    );
    assert_eq!(
        lines([KEY_OF_LINE_6, KEY_OF_LINE_6], false, None),
        [Some(6), Some(6)],
        "one key asked for twice"
    );
    let collapsed_key_of_line_5 = "1:200\tk\tt1\tTranscript\tintron_variant";
    assert_eq!(
        lines([collapsed_key_of_line_5, KEY_OF_LINE_6], true, None),
        [Some(5), Some(6)]
    );
    assert_eq!(
        lines([collapsed_key_of_line_5, KEY_OF_LINE_5], false, None),
        [None, Some(5)],
        "a key no row writes beside one the fifth line writes"
    );
    let other_contig: HashSet<String> = HashSet::from(["2".to_string()]);
    assert_eq!(
        lines([KEY_OF_LINE_6, KEY_OF_LINE_5], false, Some(&other_contig)),
        [None, None]
    );
    let absent = scratch.path.join("absent.txt");
    let keys = [KEY_OF_LINE_6, KEY_OF_LINE_5].map(str::as_bytes);
    assert!(first_lines_writing(&absent, keys, false, None).is_err());
}

/// A file that yields one of the two keys alone gives the message no line: the side's
/// output writes the key the merge stopped on, and its fifth line another Location than the
/// key before it holds.
#[test]
fn the_refusal_names_no_line_where_the_output_yields_one_of_its_two_keys_alone() {
    let run = under_way();
    let (perl, rust) = (
        run.report.path.join("perl.txt"),
        run.report.path.join("rust.txt"),
    );
    std::fs::write(&perl, OUTPUT.replace("1:200", "1:201")).unwrap();
    std::fs::write(&rust, "").unwrap();
    let stopped = Stopped::on(
        "output.txt",
        refusal_of_key_order([KEY_OF_LINE_5, KEY_OF_LINE_6]),
        &run.key_files,
        [&perl, &rust],
        None,
    );
    assert_eq!(
        stopped.message,
        format!(
            "output.txt: the sorted keys of the --perl-dir side do not ascend: a key of the \
             side does not sort after the key before it{NEEDS}"
        )
    );
}

#[test]
fn the_refusal_names_the_lines_of_its_sides_output_that_first_write_its_two_keys() {
    let run = under_way();
    let (perl, rust) = (
        run.report.path.join("perl.txt"),
        run.report.path.join("rust.txt"),
    );
    std::fs::write(&perl, OUTPUT).unwrap();
    std::fs::write(&rust, "").unwrap();
    let stopped = Stopped::on(
        "output.txt",
        refusal_of_key_order([KEY_OF_LINE_5, KEY_OF_LINE_6]),
        &run.key_files,
        [&perl, &rust],
        None,
    );
    assert_eq!(
        stopped.message,
        format!(
            "output.txt: the sorted keys of the --perl-dir side do not ascend: the key of \
             line 6 of {} does not sort after the key of line 5, which stands before it{NEEDS}",
            perl.display()
        )
    );
    assert_eq!(
        stopped.key_files.as_deref(),
        Some(run.key_files.path.as_path())
    );
}
