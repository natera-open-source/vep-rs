// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The reference release option and the two options that name an input VCF.

use super::*;

fn parsed(arguments: &[&str]) -> Result<Args, String> {
    parse_args(arguments.iter().map(OsString::from).collect())
}

/// The release the command line selects, or the message it is refused with.
fn release_of(arguments: &[&str]) -> Result<String, String> {
    parsed(arguments).map(|args| args.reference_release)
}

#[test]
fn the_reference_release_defaults_to_116_3_and_is_one_of_the_three_releases() {
    assert_eq!(DEFAULT_REFERENCE_RELEASE, "116.3");
    assert_eq!(REFERENCE_RELEASES, ["115.2", "116.2", "116.3"]);
    assert_eq!(release_of(&[]).as_deref(), Ok("116.3"));
    for release in REFERENCE_RELEASES {
        assert_eq!(
            release_of(&["--reference-release", release]).as_deref(),
            Ok(release)
        );
    }
    assert_eq!(
        release_of(&["--reference-release", "116.1"]),
        Err(
            "argument --reference-release: invalid choice: '116.1' (choose from 115.2, 116.2, 116.3)"
                .to_string()
        )
    );
}

#[test]
fn an_input_directory_and_an_input_file_are_kept_as_given() {
    let scratch = TempDir::create(&std::env::temp_dir(), "input_option_test_").unwrap();
    let inputs = scratch.path.join("inputs");
    std::fs::create_dir(&inputs).unwrap();
    let vcf = inputs.join("s.vcf");
    std::fs::write(&vcf, "").unwrap();
    let directories = ["--perl-dir", "p", "--rust-dir", "r", "--report-dir", "o"];
    let with = |option: &str, value: &Path| {
        let named = [option, value.to_str().unwrap()];
        let arguments = [&directories[..], &named[..]].concat();
        let args = parsed(&arguments).unwrap_or_else(|message| panic!("{message}"));
        (args.input_dir, args.input_vcf)
    };
    assert_eq!(with("--input-dir", &inputs), (Some(inputs.clone()), None));
    assert_eq!(with("--input-vcf", &vcf), (None, Some(vcf.clone())));
}
