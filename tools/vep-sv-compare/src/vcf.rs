// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The input VCF: its records, each with the variant type the reports group by.

use std::io::Read;
use std::os::unix::ffi::OsStrExt;
use std::path::{Path, PathBuf};

use flate2::read::MultiGzDecoder;

use crate::failure::Failure;
use crate::normalize::{location_integer_digits, normalize_location};
use crate::perl::perl_pos;
use crate::pyfmt::py_strip;
use crate::rules::gvcf::perl_split;
use crate::text::{
    first_defect, line_count, lines, open_file, os_failure, py_isalpha, py_name, py_stem,
    py_suffix, read_blocks, reading_failure, LineDefect, RunDefect, UntilFailure, BLOCK_BYTES,
    INT_MAX_STR_DIGITS,
};
use crate::threads::BlockSource;

/// One record of an input VCF as the comparison reads it.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct VcfVariant {
    pub chrom: String,
    /// POS as Ensembl VEP numifies it.
    pub pos: i128,
    /// The first of the `;`-separated identifiers of the ID column: what the first column
    /// of a row holds where an engine names the record by its ID.
    pub vid: String,
    pub reference: String,
    pub alt: String,
    pub variant_type: String,
    /// `chrom:pos`, normalised as an output's Location is.
    pub location: String,
    /// The INFO column as written; empty when the record has fewer than eight columns.
    pub info: String,
    /// The POS column as written; `None` stands for the decimal spelling of `pos`.
    pub pos_text: Option<String>,
}

/// The name of an input file without its `.vcf` or `.vcf.gz` ending, as the bytes of the
/// file name: the name the two outputs of the file are looked up under.
pub fn vcf_basename_bytes(name: &[u8]) -> &[u8] {
    let stem = py_stem(name);
    stem.strip_suffix(b".vcf").unwrap_or(stem)
}

/// The name of an input file without its `.vcf` or `.vcf.gz` ending; of a path that ends
/// in `..`, that component.
pub fn vcf_basename(vcf_path: &Path) -> String {
    let name = py_name(vcf_path.as_os_str().as_bytes());
    String::from_utf8_lossy(vcf_basename_bytes(name)).into_owned()
}

/// The variant type of a record by its REF and ALT columns.
pub fn classify_variant(reference: &str, alt: &str) -> &'static str {
    // More than one ALT field, a trailing empty one being none, comes first: the other
    // rules read one allele.
    let alts = perl_split(alt, ',');
    if alts.len() > 1 {
        return "Multi_allelic";
    }
    let alt = alts.first().copied().unwrap_or(alt);
    if alt == "*" {
        return "Spanning";
    }
    if alt == "." {
        return "RefOnly";
    }
    if alt.starts_with("<NON_REF") || alt == "<*>" {
        return "NON_REF";
    }
    if alt.starts_with('<') {
        let upper = alt.to_uppercase();
        if upper.starts_with("<DEL") {
            return "Symbolic_DEL";
        }
        if upper.starts_with("<INS") {
            return if upper.contains("ME") {
                "Mobile_Element"
            } else {
                "Symbolic_INS"
            };
        }
        if upper.starts_with("<DUP") {
            return "Symbolic_DUP";
        }
        if upper.starts_with("<INV") {
            return "Symbolic_INV";
        }
        if upper.starts_with("<CN") {
            return "Symbolic_CNV";
        }
        if upper.starts_with("<CPX") {
            return "Complex_SV";
        }
        // A symbolic breakend with INFO/CHR2 and END2 is a breakend record like its
        // bracket-notation form.
        if upper.starts_with("<BND") {
            return "BND";
        }
        return "Other";
    }
    if alt.contains('[') || alt.contains(']') {
        return "BND";
    }
    // A single breakend: letters before or after one dot.
    if alt.strip_suffix('.').is_some_and(py_isalpha)
        || alt.strip_prefix('.').is_some_and(py_isalpha)
    {
        return "BND";
    }
    let bases = |allele: &str| {
        allele
            .bytes()
            .all(|b| matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T'))
    };
    if bases(reference) && bases(alt) {
        let (lr, la) = (reference.len(), alt.len());
        return if lr == 1 && la == 1 {
            "SNV"
        } else if lr == la && lr > 1 {
            "MNP"
        } else if lr > la {
            "Small_DEL"
        } else if lr < la {
            "Small_INS"
        } else {
            "Complex"
        };
    }
    "Other"
}

/// The digits of a POS column Python's integer conversion refuses: a decimal integer,
/// blanks around it and a sign before it allowed, of more than [`INT_MAX_STR_DIGITS`]
/// digits.
fn pos_integer_digits(pos_text: &str) -> Option<usize> {
    if pos_text.len() <= INT_MAX_STR_DIGITS {
        return None;
    }
    let stripped = py_strip(pos_text);
    let digits = stripped.strip_prefix(['+', '-']).unwrap_or(stripped);
    (digits.len() > INT_MAX_STR_DIGITS && digits.bytes().all(|b| b.is_ascii_digit()))
        .then_some(digits.len())
}

/// The first of the `;`-separated identifiers of an ID column, empty where the column
/// opens with a `;`; the whole column where it holds no `;`.
fn first_identifier(ids: &str) -> &str {
    ids.split_once(';').map_or(ids, |(first, _)| first)
}

/// The record a line of an input VCF holds: `None` for a header line and for a line of
/// fewer than five columns.
fn parse_line(line: &str) -> Option<Result<VcfVariant, LineDefect>> {
    if line.starts_with('#') {
        return None;
    }
    let mut columns = line.split('\t');
    let (chrom, pos_text, ids, reference, alt) = (
        columns.next()?,
        columns.next()?,
        columns.next()?,
        columns.next()?,
        columns.next()?,
    );
    let info = columns.nth(2).unwrap_or("");
    if let Some(digits) = pos_integer_digits(pos_text) {
        return Some(Err(LineDefect::IntegerDigits(digits)));
    }
    let pos = perl_pos(pos_text);
    let written = format!("{chrom}:{pos}");
    if let Some(digits) = location_integer_digits(&written) {
        return Some(Err(LineDefect::IntegerDigits(digits)));
    }
    Some(Ok(VcfVariant {
        chrom: chrom.to_string(),
        pos,
        vid: first_identifier(ids).to_string(),
        reference: reference.to_string(),
        alt: alt.to_string(),
        variant_type: classify_variant(reference, alt).to_string(),
        location: normalize_location(&written),
        info: info.to_string(),
        pos_text: Some(pos_text.to_string()),
    }))
}

/// What one run of whole lines of an input VCF holds.
#[derive(Debug, Default)]
pub struct VcfBlock {
    /// How many lines the run has.
    pub lines: u64,
    /// The records of the lines before the first defective one.
    pub variants: Vec<VcfVariant>,
    pub defect: Option<RunDefect>,
}

/// Reads the records of one run of whole lines, up to its first defective line.
pub fn parse_vcf_block(block: &[u8]) -> VcfBlock {
    let mut parsed = VcfBlock {
        lines: line_count(block),
        ..VcfBlock::default()
    };
    for (index, line) in lines(block).enumerate() {
        let record = match std::str::from_utf8(line) {
            Ok(line) => parse_line(line),
            Err(_) => Some(Err(LineDefect::NotUtf8)),
        };
        match record {
            None => {}
            Some(Ok(variant)) => parsed.variants.push(variant),
            Some(Err(defect)) => {
                parsed.defect = Some((index as u64 + 1, defect));
                break;
            }
        }
    }
    parsed
}

/// An input VCF as a source of runs of lines.
pub struct VcfSource {
    pub path: PathBuf,
    /// The bytes a run holds, about.
    pub block_bytes: usize,
}

impl VcfSource {
    pub fn new(path: PathBuf) -> VcfSource {
        VcfSource {
            path,
            block_bytes: BLOCK_BYTES,
        }
    }

    /// Whether the file is read through a gzip decoder: its name has the suffix `.gz`.
    fn gzipped(&self) -> bool {
        py_suffix(py_name(self.path.as_os_str().as_bytes())) == b".gz"
    }

    /// The bytes of the file as they are read: those of a gzipped file as `zcat` writes
    /// them to a reader that ignores how it ends, every member of the stream and nothing
    /// after the point where it cannot be decoded, in silence.
    fn open(&self) -> std::io::Result<Box<dyn Read>> {
        let file = open_file(&self.path)?;
        Ok(if self.gzipped() {
            Box::new(UntilFailure::new(MultiGzDecoder::new(file)))
        } else {
            Box::new(file)
        })
    }

    /// The records of the file from its parsed runs, or the failure its reading ends in
    /// when a line of it is defective.
    pub fn variants(
        &self,
        read: Result<(), Failure>,
        blocks: Vec<VcfBlock>,
    ) -> Result<Vec<VcfVariant>, Failure> {
        // A line is read, and can fail, before the file is known to end badly.
        let defect = first_defect(blocks.iter().map(|block| (block.lines, block.defect)));
        if let Some((line, defect)) = defect {
            return Err(reading_failure(&self.path, self.open(), line, defect));
        }
        read?;
        Ok(blocks
            .into_iter()
            .flat_map(|block| block.variants)
            .collect())
    }
}

impl BlockSource for VcfSource {
    type Read = Result<(), Failure>;
    type Parsed = VcfBlock;

    /// A gzipped file that cannot be opened holds nothing, in silence. A plain file that
    /// cannot be opened or read is a failure.
    fn read(&self, emit: &mut dyn FnMut(Vec<u8>)) -> Result<(), Failure> {
        match (self.open(), self.gzipped()) {
            (Ok(stream), _) => read_blocks(stream, self.block_bytes, emit)
                .map_err(|error| os_failure(&error, None)),
            (Err(_), true) => Ok(()),
            (Err(error), false) => Err(os_failure(&error, Some(&self.path))),
        }
    }

    fn parse(&self, _sequence: u64, block: Vec<u8>) -> VcfBlock {
        parse_vcf_block(&block)
    }
}

/// The records of an input VCF, plain or gzipped, each with its variant type. POS is read
/// as Ensembl VEP numifies it, so a record whose POS is not an integer is kept and
/// measured.
pub fn parse_vcf(path: &Path) -> Result<Vec<VcfVariant>, Failure> {
    let source = VcfSource::new(path.to_path_buf());
    let mut blocks = Vec::new();
    let read = source.read(&mut |block| blocks.push(parse_vcf_block(&block)));
    source.variants(read, blocks)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::failure::FailureKind;
    use crate::TempDir;
    use flate2::write::GzEncoder;
    use flate2::Compression;
    use std::io::Write;

    fn gzip(bytes: &[u8]) -> Vec<u8> {
        let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
        encoder.write_all(bytes).unwrap();
        encoder.finish().unwrap()
    }

    #[test]
    fn test_variant_classification() {
        assert_eq!(classify_variant("A", "T"), "SNV");
        assert_eq!(classify_variant("N", "<DEL>"), "Symbolic_DEL");
        assert_eq!(classify_variant("N", "<DUP:TANDEM>"), "Symbolic_DUP");
        assert_eq!(classify_variant("N", "<CN=3>"), "Symbolic_CNV");
        assert_eq!(classify_variant("N", "<INS:ME:ALU>"), "Mobile_Element");
        assert_eq!(classify_variant("N", "N[21:100["), "BND");
        assert_eq!(classify_variant("N", "<BND>"), "BND");
        assert_eq!(classify_variant("N", "<NON_REF>"), "NON_REF");
        assert_eq!(classify_variant("A", "*"), "Spanning");
    }

    #[test]
    fn a_record_is_typed_by_the_first_rule_its_alleles_meet() {
        let cases = [
            ("A", "", "Small_DEL"),
            ("", "", "Complex"),
            ("A", ",", "Other"),
            ("A", ",,", "Other"),
            ("A", "T,", "SNV"),
            ("A", "T,,C", "Multi_allelic"),
            ("A", ",T", "Multi_allelic"),
            ("A", "<in\u{17f}>", "Symbolic_INS"),
            ("A", "<ins:me>", "Mobile_Element"),
            ("A", "<INSME>", "Mobile_Element"),
            ("A", "<DEL>,", "Symbolic_DEL"),
            ("A", "A.", "BND"),
            ("A", ".A", "BND"),
            ("A", ".", "RefOnly"),
            ("A", "\u{2167}.", "Other"),
            ("A", "\u{e9}.", "BND"),
            ("A", "A1.", "Other"),
            ("A", "..A", "Other"),
            ("A", "a\u{301}.", "Other"),
            ("N", "<NON_REF>x", "NON_REF"),
            ("N", "<non_ref>", "Other"),
            ("N", "<*>", "NON_REF"),
            ("N", "<* >", "Other"),
            ("AC", "A", "Small_DEL"),
            ("A", "AC", "Small_INS"),
            ("AC", "GT", "MNP"),
            ("ac", "gt", "MNP"),
            ("N", "A", "Other"),
            ("A", "N", "Other"),
            ("", "A", "Small_INS"),
            ("A", "<CN0>", "Symbolic_CNV"),
            ("A", "<cn=2>", "Symbolic_CNV"),
            ("A", "<CPX>", "Complex_SV"),
            ("A", "<bnd>", "BND"),
            ("A", "<X>", "Other"),
            ("A", "]1:5]A", "BND"),
            ("A", "\u{ff41}.", "BND"),
            ("A", "A\tB", "Other"),
        ];
        for (reference, alt, expected) in cases {
            assert_eq!(
                classify_variant(reference, alt),
                expected,
                "{reference:?} {alt:?}"
            );
        }
    }

    #[test]
    fn a_file_name_loses_its_vcf_ending() {
        let cases = [
            ("d/01_snv.vcf.gz", "01_snv"),
            ("d/x.vcf", "x"),
            (".vcf.gz", ""),
            (".vcf", ""),
            ("x.vcf.vcf", "x"),
            ("a.b.vcf.gz", "a.b"),
            ("x..vcf.gz", "x."),
            ("...vcf", ".."),
            ("..vcf.gz", "."),
            ("x.VCF.gz", "x.VCF"),
            ("a/..", ".."),
            ("..", ".."),
            ("a/../..", ".."),
            ("a/.", "a"),
            ("a/", "a"),
            (".", ""),
            ("", ""),
            ("/", ""),
            ("x/...gz", "...gz"),
            ("x/..a.gz", "..a"),
            ("x/foo.", "foo"),
            ("x/a.vcf.vcf.gz", "a.vcf"),
        ];
        for (path, expected) in cases {
            assert_eq!(vcf_basename(Path::new(path)), expected, "{path}");
        }
    }

    #[test]
    fn a_plain_file_is_read_record_by_record() {
        let dir = TempDir::new("vcf_plain_").unwrap();
        let path = dir.path.join("in.vcf");
        std::fs::write(
            &path,
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\n\
             chr21\t100\tid1\tA\tT\t.\t.\tEND=5\tGT\n\
             chrM\t 3e2\t.\tN\t<DEL>\r\n\
             21\t7\tshort\tA\n\
             \n\
             22\tabc\tid3\tN\tN[chr1:5[\t.\t.\rX\t+12\tid4\tAC\tA\t.\tPASS\t.",
        )
        .unwrap();
        let variants = parse_vcf(&path).unwrap();
        let summary: Vec<String> = variants
            .iter()
            .map(|v| {
                format!(
                    "{}|{}|{}|{}|{}|{}|{:?}",
                    v.chrom, v.pos, v.vid, v.variant_type, v.location, v.info, v.pos_text
                )
            })
            .collect();
        assert_eq!(
            summary,
            [
                r#"chr21|100|id1|SNV|21:100|END=5|Some("100")"#,
                r#"chrM|300|.|Symbolic_DEL|MT:300||Some(" 3e2")"#,
                r#"22|0|id3|BND|22:0||Some("abc")"#,
                r#"X|12|id4|Small_DEL|X:12|.|Some("+12")"#,
            ]
        );
        assert_eq!(variants[0].reference, "A");
        assert_eq!(variants[2].alt, "N[chr1:5[");
    }

    #[test]
    fn a_record_is_named_by_the_first_identifier_of_its_id_column() {
        let dir = TempDir::new("vcf_ids_").unwrap();
        let path = dir.path.join("in.vcf");
        // Ensembl VEP and vep-rs split the column at `;` and name a record by the first
        // part alone; an empty first part and `.` name none.
        let cases = [
            ("one;two", "one"),
            ("one;two;three", "one"),
            ("one;", "one"),
            ("one", "one"),
            (";two", ""),
            (";", ""),
            ("", ""),
            (".", "."),
            (".;two", "."),
        ];
        let text: String = cases
            .iter()
            .enumerate()
            .map(|(n, (ids, _))| format!("21\t{}\t{ids}\tA\tT\n", n + 1))
            .collect();
        std::fs::write(&path, text).unwrap();
        let named: Vec<String> = parse_vcf(&path)
            .unwrap()
            .into_iter()
            .map(|v| v.vid)
            .collect();
        let expected: Vec<&str> = cases.iter().map(|(_, first)| *first).collect();
        assert_eq!(named, expected);
    }

    #[test]
    fn blocks_of_any_size_hold_the_records_of_the_file() {
        let dir = TempDir::new("vcf_blocks_").unwrap();
        let path = dir.path.join("in.vcf");
        let mut text = String::from("##header\n");
        for n in 0..200 {
            text.push_str(&format!("21\t{n}\tid{n}\tA\tT\t.\t.\tK={n}\r\n"));
        }
        std::fs::write(&path, &text).unwrap();
        let whole = parse_vcf(&path).unwrap();
        assert_eq!(whole.len(), 200);
        for block_bytes in [1, 7, 64, 1000] {
            let source = VcfSource {
                path: path.clone(),
                block_bytes,
            };
            for threads in [1, 4] {
                let (read, blocks) =
                    crate::threads::read_sources(threads, std::slice::from_ref(&source))
                        .pop()
                        .unwrap();
                assert_eq!(source.variants(read, blocks).unwrap(), whole);
            }
        }
    }

    #[test]
    fn a_gzipped_file_is_every_member_up_to_where_the_stream_fails() {
        let dir = TempDir::new("vcf_gz_").unwrap();
        let first = gzip(b"#h\n21\t1\ta\tA\tT\n21\t2\tb\tA\tT\n");
        let second = gzip(b"21\t3\tc\tA\tT\n21\t4\td\tA\tT\n");
        let ids = |name: &str, bytes: &[u8]| -> Vec<String> {
            let path = dir.path.join(name);
            std::fs::write(&path, bytes).unwrap();
            parse_vcf(&path)
                .unwrap()
                .into_iter()
                .map(|v| v.vid)
                .collect()
        };
        let both = [first.clone(), second.clone()].concat();
        assert_eq!(ids("both.vcf.gz", &both), ["a", "b", "c", "d"]);
        // Bytes after the last member that are no member: the members before stand.
        let trailing = [both.clone(), b"not a member".to_vec()].concat();
        assert_eq!(ids("trailing.vcf.gz", &trailing), ["a", "b", "c", "d"]);
        // A second member cut short inside its header: the first member stands.
        let cut = [first.clone(), second[..6].to_vec()].concat();
        assert_eq!(ids("cut.vcf.gz", &cut), ["a", "b"]);
        // A member whose checksum is wrong is read, and nothing after it.
        let mut wrong = first.clone();
        let checksum = wrong.len() - 8;
        wrong[checksum] ^= 0xff;
        let wrong = [wrong, second.clone()].concat();
        assert_eq!(ids("wrong.vcf.gz", &wrong), ["a", "b"]);
        assert!(ids("plain.vcf.gz", b"21\t1\ta\tA\tT\n").is_empty());
        assert!(ids("empty.vcf.gz", b"").is_empty());
        assert!(parse_vcf(&dir.path.join("absent.vcf.gz"))
            .unwrap()
            .is_empty());
        std::fs::create_dir(dir.path.join("folder.vcf.gz")).unwrap();
        assert!(parse_vcf(&dir.path.join("folder.vcf.gz"))
            .unwrap()
            .is_empty());
    }

    #[test]
    fn a_name_of_dots_before_gz_has_no_suffix_and_is_read_plain() {
        let dir = TempDir::new("vcf_dots_").unwrap();
        let record: &[u8] = b"21\t1\ta\tA\tT\n";
        let count = |name: &str, bytes: &[u8]| {
            let path = dir.path.join(name);
            std::fs::write(&path, bytes).unwrap();
            parse_vcf(&path).unwrap().len()
        };
        for plain in ["...gz", "..gz", ".gz"] {
            assert_eq!(count(plain, record), 1, "{plain}");
        }
        for gzipped in ["..a.gz", ".vcf.gz", "a..gz"] {
            assert_eq!(count(gzipped, record), 0, "{gzipped}");
            assert_eq!(count(gzipped, &gzip(record)), 1, "{gzipped}");
        }
    }

    #[test]
    fn a_plain_file_that_cannot_be_read_is_a_failure() {
        let dir = TempDir::new("vcf_fail_").unwrap();
        let path = dir.path.join("absent.vcf");
        let absent = parse_vcf(&path).unwrap_err();
        assert_eq!(absent.kind, FailureKind::Io);
        assert_eq!(
            absent.message,
            format!("[Errno 2] No such file or directory: '{}'", path.display())
        );
        let path = dir.path.join("folder.vcf");
        std::fs::create_dir(&path).unwrap();
        let folder = parse_vcf(&path).unwrap_err();
        assert_eq!(folder.kind, FailureKind::Io);
        assert_eq!(
            folder.message,
            format!("[Errno 21] Is a directory: '{}'", path.display())
        );
    }

    #[test]
    fn a_line_that_is_not_utf8_is_a_failure_in_either_kind_of_file() {
        let dir = TempDir::new("vcf_utf8_").unwrap();
        let bytes = b"#h\n21\t1\ta\tA\tT\n21\t2\t\xff\tA\tT\n21\t3\tc\tA\tT\n";
        let plain = dir.path.join("in.vcf");
        std::fs::write(&plain, bytes).unwrap();
        let expected = "'utf-8' codec can't decode byte 0xff in position 19: invalid start byte";
        let failure = parse_vcf(&plain).unwrap_err();
        assert_eq!(failure.kind, FailureKind::Decode);
        assert_eq!(failure.message, expected);
        let gz = dir.path.join("in.vcf.gz");
        std::fs::write(&gz, gzip(bytes)).unwrap();
        let failure = parse_vcf(&gz).unwrap_err();
        assert_eq!(failure.kind, FailureKind::Decode);
        assert_eq!(failure.message, expected);
    }

    #[test]
    fn a_pos_past_the_integer_limit_is_a_failure() {
        let dir = TempDir::new("vcf_int_").unwrap();
        let path = dir.path.join("in.vcf");
        let digits = "1".repeat(INT_MAX_STR_DIGITS + 1);
        std::fs::write(&path, format!("21\t {digits}\ta\tA\tT\n")).unwrap();
        let failure = parse_vcf(&path).unwrap_err();
        assert_eq!(failure.kind, FailureKind::Value);
        assert!(failure.message.contains("value has 4301 digits"));
        std::fs::write(&path, format!("21\t{}\ta\tA\tT\n", &digits[1..])).unwrap();
        assert_eq!(parse_vcf(&path).unwrap()[0].pos, i128::MAX);
    }
}
