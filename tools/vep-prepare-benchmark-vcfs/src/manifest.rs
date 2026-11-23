// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `manifest.json`: the run and its files as `json.dumps(indent=2)` writes them, with no
//! newline after the closing brace.

use crate::lines::Counts;
use crate::pyfmt::{push_json_os_str, push_json_str};
use std::fmt::Write;
use std::time::{SystemTime, UNIX_EPOCH};

/// One prepared file.
pub struct FileSummary {
    pub input_file: Vec<u8>,
    pub output_file: Vec<u8>,
    pub input_style: &'static str,
    pub cache_style: &'static str,
    pub mode: &'static str,
    pub counts: Counts,
    pub sha256: String,
}

/// The run: every value in the bytes the command line or the file system gave.
pub struct Manifest<'a> {
    pub generated_at_utc: String,
    pub input_dir: &'a [u8],
    pub output_dir: &'a [u8],
    pub requested_chrom_mode: &'static str,
    pub cache_dir: Option<&'a [u8]>,
    pub max_variants: Option<&'a str>,
    pub assembly: Option<&'a [u8]>,
    pub policy: &'static str,
    pub contigs: Option<&'a [Vec<u8>]>,
    pub files: &'a [FileSummary],
}

fn key(out: &mut String, indent: usize, name: &str) {
    for _ in 0..indent {
        out.push(' ');
    }
    push_json_str(out, name);
    out.push_str(": ");
}

fn os_str_or_null(out: &mut String, value: Option<&[u8]>) {
    match value {
        Some(value) => push_json_os_str(out, value),
        None => out.push_str("null"),
    }
}

impl Manifest<'_> {
    pub fn dumps(&self) -> String {
        let mut out = String::from("{\n");
        key(&mut out, 2, "generated_at_utc");
        push_json_str(&mut out, &self.generated_at_utc);
        out.push_str(",\n");
        key(&mut out, 2, "input_dir");
        push_json_os_str(&mut out, self.input_dir);
        out.push_str(",\n");
        key(&mut out, 2, "output_dir");
        push_json_os_str(&mut out, self.output_dir);
        out.push_str(",\n");
        key(&mut out, 2, "requested_chrom_mode");
        push_json_str(&mut out, self.requested_chrom_mode);
        out.push_str(",\n");
        key(&mut out, 2, "cache_dir");
        os_str_or_null(&mut out, self.cache_dir);
        out.push_str(",\n");
        key(&mut out, 2, "max_variants");
        out.push_str(self.max_variants.unwrap_or("null"));
        out.push_str(",\n");
        key(&mut out, 2, "assembly");
        os_str_or_null(&mut out, self.assembly);
        out.push_str(",\n");
        key(&mut out, 2, "canonical_contigs_policy");
        push_json_str(&mut out, self.policy);
        out.push_str(",\n");
        key(&mut out, 2, "canonical_contigs_filter_active");
        out.push_str(if self.contigs.is_some() {
            "true"
        } else {
            "false"
        });
        out.push_str(",\n");
        key(&mut out, 2, "canonical_contigs_set");
        match self.contigs {
            Some(names) => {
                out.push('[');
                for (index, name) in names.iter().enumerate() {
                    out.push_str(if index == 0 { "\n    " } else { ",\n    " });
                    push_json_os_str(&mut out, name);
                }
                out.push_str("\n  ]");
            }
            None => out.push_str("null"),
        }
        out.push_str(",\n");
        key(&mut out, 2, "files");
        out.push('[');
        for (index, file) in self.files.iter().enumerate() {
            out.push_str(if index == 0 {
                "\n    {\n"
            } else {
                ",\n    {\n"
            });
            file.dump_into(&mut out);
            out.push_str("\n    }");
        }
        out.push_str(if self.files.is_empty() {
            "]\n}"
        } else {
            "\n  ]\n}"
        });
        out
    }
}

impl FileSummary {
    fn dump_into(&self, out: &mut String) {
        key(out, 6, "input_file");
        push_json_os_str(out, &self.input_file);
        out.push_str(",\n");
        key(out, 6, "output_file");
        push_json_os_str(out, &self.output_file);
        out.push_str(",\n");
        key(out, 6, "detected_input_chrom_style");
        push_json_str(out, self.input_style);
        out.push_str(",\n");
        key(out, 6, "detected_cache_chrom_style");
        push_json_str(out, self.cache_style);
        out.push_str(",\n");
        key(out, 6, "normalization_mode");
        push_json_str(out, self.mode);
        out.push_str(",\n");
        key(out, 6, "headers_written");
        writeln!(out, "{},", self.counts.headers).unwrap();
        key(out, 6, "variants_written");
        writeln!(out, "{},", self.counts.variants).unwrap();
        key(out, 6, "variants_dropped_non_canonical");
        writeln!(out, "{},", self.counts.dropped).unwrap();
        key(out, 6, "output_sha256");
        push_json_str(out, &self.sha256);
        out.push_str(",\n");
        key(out, 6, "chr_stripped");
        out.push_str(if self.counts.chr_stripped {
            "true"
        } else {
            "false"
        });
    }
}

/// The civil date of a day counted from 1970-01-01, in the proleptic Gregorian calendar.
fn civil_from_days(days: i64) -> (i64, u32, u32) {
    let shifted = days + 719_468;
    let era = shifted.div_euclid(146_097);
    let day_of_era = shifted.rem_euclid(146_097);
    let year_of_era =
        (day_of_era - day_of_era / 1_460 + day_of_era / 36_524 - day_of_era / 146_096) / 365;
    let day_of_year = day_of_era - (365 * year_of_era + year_of_era / 4 - year_of_era / 100);
    let month_index = (5 * day_of_year + 2) / 153;
    let day = (day_of_year - (153 * month_index + 2) / 5 + 1) as u32;
    let month = if month_index < 10 {
        month_index + 3
    } else {
        month_index - 9
    } as u32;
    let year = year_of_era + era * 400 + i64::from(month <= 2);
    (year, month, day)
}

/// `datetime.isoformat()` of a UTC time: seconds, six decimals where the microseconds are
/// not zero, and the `+00:00` offset.
pub fn utc_isoformat(time: SystemTime) -> String {
    let since_epoch = time.duration_since(UNIX_EPOCH).unwrap_or_default();
    let seconds = since_epoch.as_secs() as i64;
    let micros = since_epoch.subsec_micros();
    let (year, month, day) = civil_from_days(seconds.div_euclid(86_400));
    let second_of_day = seconds.rem_euclid(86_400);
    let mut out = format!(
        "{year:04}-{month:02}-{day:02}T{:02}:{:02}:{:02}",
        second_of_day / 3_600,
        second_of_day % 3_600 / 60,
        second_of_day % 60
    );
    if micros != 0 {
        write!(out, ".{micros:06}").unwrap();
    }
    out.push_str("+00:00");
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::Duration;

    fn at(seconds: u64, micros: u32) -> String {
        utc_isoformat(UNIX_EPOCH + Duration::new(seconds, micros * 1_000))
    }

    #[test]
    fn a_time_is_written_as_isoformat_writes_it() {
        assert_eq!(at(0, 0), "1970-01-01T00:00:00+00:00");
        assert_eq!(at(0, 1), "1970-01-01T00:00:00.000001+00:00");
        assert_eq!(at(951_782_399, 999_999), "2000-02-28T23:59:59.999999+00:00");
        assert_eq!(at(951_782_400, 120_000), "2000-02-29T00:00:00.120000+00:00");
        assert_eq!(at(1_709_251_200, 0), "2024-03-01T00:00:00+00:00");
        assert_eq!(at(1_735_689_599, 5), "2024-12-31T23:59:59.000005+00:00");
        assert_eq!(at(4_102_444_800, 0), "2100-01-01T00:00:00+00:00");
        // Below a microsecond nothing is written.
        assert_eq!(
            utc_isoformat(UNIX_EPOCH + Duration::new(86_400, 999)),
            "1970-01-02T00:00:00+00:00"
        );
    }

    fn summary(name: &str, counts: Counts) -> FileSummary {
        FileSummary {
            input_file: format!("/in/{name}.gz").into_bytes(),
            output_file: format!("/out/{name}").into_bytes(),
            input_style: "chr",
            cache_style: "unknown",
            mode: "keep",
            counts,
            sha256: "00ff".to_string(),
        }
    }

    #[test]
    fn the_manifest_is_laid_out_as_json_dumps_lays_it_out() {
        let contigs = [b"1".to_vec(), b"sca\xfffold".to_vec()];
        let files = [
            summary(
                "a.vcf",
                Counts {
                    headers: 3,
                    variants: 2,
                    dropped: 1,
                    chr_stripped: true,
                },
            ),
            summary("b \"quoted\".vcf", Counts::default()),
        ];
        let manifest = Manifest {
            generated_at_utc: "2024-03-01T00:00:00+00:00".to_string(),
            input_dir: b"/in",
            output_dir: b"/out",
            requested_chrom_mode: "auto",
            cache_dir: Some(b"/cache"),
            max_variants: Some("-5"),
            assembly: Some(b"GRCh38"),
            policy: "named_set",
            contigs: Some(&contigs),
            files: &files,
        };
        let expected = r#"{
  "generated_at_utc": "2024-03-01T00:00:00+00:00",
  "input_dir": "/in",
  "output_dir": "/out",
  "requested_chrom_mode": "auto",
  "cache_dir": "/cache",
  "max_variants": -5,
  "assembly": "GRCh38",
  "canonical_contigs_policy": "named_set",
  "canonical_contigs_filter_active": true,
  "canonical_contigs_set": [
    "1",
    "sca\udcfffold"
  ],
  "files": [
    {
      "input_file": "/in/a.vcf.gz",
      "output_file": "/out/a.vcf",
      "detected_input_chrom_style": "chr",
      "detected_cache_chrom_style": "unknown",
      "normalization_mode": "keep",
      "headers_written": 3,
      "variants_written": 2,
      "variants_dropped_non_canonical": 1,
      "output_sha256": "00ff",
      "chr_stripped": true
    },
    {
      "input_file": "/in/b \"quoted\".vcf.gz",
      "output_file": "/out/b \"quoted\".vcf",
      "detected_input_chrom_style": "chr",
      "detected_cache_chrom_style": "unknown",
      "normalization_mode": "keep",
      "headers_written": 0,
      "variants_written": 0,
      "variants_dropped_non_canonical": 0,
      "output_sha256": "00ff",
      "chr_stripped": false
    }
  ]
}"#;
        assert_eq!(manifest.dumps(), expected);
    }

    #[test]
    fn what_a_run_does_not_have_is_null() {
        let files = [summary("a.vcf", Counts::default())];
        let manifest = Manifest {
            generated_at_utc: "t".to_string(),
            input_dir: b"/in",
            output_dir: b"/out",
            requested_chrom_mode: "keep",
            cache_dir: None,
            max_variants: None,
            assembly: None,
            policy: "disabled",
            contigs: None,
            files: &files,
        };
        let text = manifest.dumps();
        for line in [
            "  \"cache_dir\": null,\n",
            "  \"max_variants\": null,\n",
            "  \"assembly\": null,\n",
            "  \"canonical_contigs_filter_active\": false,\n",
            "  \"canonical_contigs_set\": null,\n",
        ] {
            assert!(text.contains(line), "{line}");
        }
        assert!(text.ends_with("\n    }\n  ]\n}"));
    }
}
