// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The parameters of a `--plugin Name,param1,param2` argument, in the three
//! grammars Ensembl VEP's plugins read them by.
//!
//! `params_to_hash` (ensembl-variation Utils/BaseVepPlugin.pm:115-130) reads
//! every parameter as `key=value` once any parameter carries `=`, and one
//! without `=` beside one with it fails with `Failed to parse parameter <token>`
//! (`:126-127`); when none carries `=` the hash is empty (`:123`) and a plugin
//! that also reads by position takes the tokens in order. The grammars:
//!
//! - [`Grammar::PositionalOrKeyValue`]: that reading beside a positional one
//!   (REVEL.pm:105-110, CADD.pm:118-131, LoFtool.pm:72, pLI.pm:104-109). A token
//!   of the grammar's flag set is a flag of the positional reading, not a value:
//!   CADD's trailing `1` is `force_annotate` (CADD.pm:119-120) and pLI's `sigfig`
//!   anywhere selects three significant figures (pLI.pm:105-106).
//! - [`Grammar::KeyValueOnly`]: `params_to_hash` alone (AlphaMissense.pm:175-181,
//!   GWAS.pm:94-96, SpliceAI.pm:160-163). A list without `=` is the data file
//!   given without its key, which those modules refuse as a missing file.
//! - [`Grammar::PositionalAndKeyValue`]: positional tokens and `key=value` pairs
//!   side by side (dbNSFP.pm:192, 586).
//!
//! [`PluginParams::parse`] reads the first grammar with no flag token, the
//! reading most plugins share; [`PluginParams::parse_with`] takes the grammar.
//! Both add two refusals of their own: a key given twice and a key the plugin
//! does not declare, where Ensembl VEP keeps the later value (`:128`) and never
//! reads the unknown key, except dbNSFP, which dies on it (dbNSFP.pm:601);
//! [`PluginParams::parse_any_key`] drops the second refusal for a module that
//! opens every key it does not name as a data file (CADD.pm:124-129). A
//! `key=value` token is split as `perl_key_value` splits it, the reading of
//! `my ($key, $val) = split('=', $param)` (`:126`, NearestExonJB.pm:101,
//! dbNSFP.pm:589): the value is the second `=`-field alone, so `file=a=b` is
//! `file` and `a`, and `file==b` and `file=` carry the empty value.
//!
//! A positional token past the ones a module reads by index is never read
//! there (pLI.pm:108-109, AncestralAllele.pm:90, Conservation.pm:91,
//! MaxEntScan.pm:110-115; a module reading no parameter at all, Blosum62.pm,
//! Downstream.pm, NMD.pm, HGVSIntronOffset.pm): the run goes on without it.
//! `warn_unread_parameters` logs each such token so a misspelt flag is seen;
//! `parse_none` is the whole reading of a module that takes no parameter.

use crate::perl_format::{numify_arithmetic, perl_number};
use crate::traits::PluginError;

/// Perl's number of the parameter `key`'s `value`, as a numeric comparison or
/// arithmetic in the plugin `name`'s module reads it (`$CONFIG{intronic} == 1`,
/// NearestExonJB.pm:135; `$cutoff < 0`, SpliceAI.pm:167): the infinity or
/// not-a-number Perl's `grok_number` flags, else its leading number, 0 when
/// none leads it ([`numify_arithmetic`]). A value that is not a number through
/// to its end is logged, the mark Perl's `isn't numeric` warning gives it in a
/// module that enables warnings (UTRAnnotator.pm does not, and prints none).
pub(crate) fn perl_numeric(name: &str, key: &str, value: &str) -> f64 {
    let number = numify_arithmetic(value).unwrap_or(0.0);
    if !perl_whole_number(value) {
        tracing::warn!(
            "{name}: {key}='{value}' is not a number; it is read as {}, as Ensembl VEP reads it",
            perl_number(number)
        );
    }
    number
}

/// Whether `value` is wholly a number to Perl (`looks_like_number`, its
/// `grok_number`): blanks of Perl's `isSPACE` at either end, one sign, then
/// digits with an optional fraction and exponent, a fraction alone, `inf` or
/// `infinity`, `nan` with an optional `q` or `s` on either side and a
/// parenthesised payload of digits or a `0x` or `0b` number, or Visual C's
/// `1.#` before `inf` or `ind` with trailing zeros or before the `nan` forms,
/// each word in any case; or the ten bytes `0 but true`. Every other text is
/// numified under a warning there: `--1` (two signs), `0x1A`, `1e`,
/// `nan(abc)` and the empty value.
fn perl_whole_number(value: &str) -> bool {
    if value == "0 but true" {
        return true;
    }
    let perl_space = |c: char| matches!(c, ' ' | '\t' | '\n' | '\x0b' | '\x0c' | '\r');
    let text = value.trim_matches(perl_space);
    let body = text
        .strip_prefix(['+', '-'])
        .unwrap_or(text)
        .to_ascii_lowercase();
    let (visual_c, word) = match body.strip_prefix("1.#").or_else(|| body.strip_prefix("1#")) {
        Some(word) => (true, word),
        None => (false, body.as_str()),
    };
    if word == "inf" || word == "infinity" {
        return true;
    }
    if visual_c && !word.is_empty() && matches!(word.trim_end_matches('0'), "inf" | "ind") {
        return true;
    }
    let digits = |s: &str| !s.is_empty() && s.bytes().all(|b| b.is_ascii_digit());
    let quiet_or_signalling = ['q', 's'];
    let unmarked = word.strip_prefix(quiet_or_signalling).unwrap_or(word);
    if let Some(tail) = unmarked.strip_prefix("nan") {
        let tail = tail.strip_prefix(quiet_or_signalling).unwrap_or(tail);
        let Some(payload) = tail.strip_prefix('(').and_then(|t| t.strip_suffix(')')) else {
            return tail.is_empty();
        };
        return if let Some(hex) = payload.strip_prefix("0x") {
            !hex.is_empty() && hex.bytes().all(|b| b.is_ascii_hexdigit())
        } else if let Some(bits) = payload.strip_prefix("0b") {
            !bits.is_empty() && bits.bytes().all(|b| matches!(b, b'0' | b'1'))
        } else {
            digits(payload)
        };
    }
    if visual_c {
        return false;
    }
    let (mantissa, exponent) = match word.split_once('e') {
        Some((mantissa, exponent)) => (mantissa, Some(exponent)),
        None => (word, None),
    };
    let mantissa_ok = match mantissa.split_once('.') {
        None => digits(mantissa),
        Some((whole, fraction)) => {
            (digits(whole) && (fraction.is_empty() || digits(fraction)))
                || (whole.is_empty() && digits(fraction))
        }
    };
    let exponent_ok = exponent
        .is_none_or(|exponent| digits(exponent.strip_prefix(['+', '-']).unwrap_or(exponent)));
    mantissa_ok && exponent_ok
}

/// Splits a `key=value` token as `my ($key, $val) = split('=', $param)` splits
/// it (Utils/BaseVepPlugin.pm:126): the key before the first `=`, the value
/// between the first `=` and the second, every later field dropped. `file=a=b`
/// is `file` and `a`; `file==b` and `file=` are `file` and the empty value.
/// `None` for a token without `=`.
pub(crate) fn perl_key_value(param: &str) -> Option<(&str, &str)> {
    let (key, rest) = param.split_once('=')?;
    let value = rest.split_once('=').map_or(rest, |(value, _)| value);
    Some((key, value))
}

/// Logs each token of `unread` as a parameter the plugin `name` does not
/// read, naming what it reads instead; the plugin then runs as its Ensembl VEP
/// module runs, which never looks at the token.
pub(crate) fn warn_unread_parameters(name: &str, unread: &[String], reads: &str) {
    for token in unread {
        tracing::warn!("{name}: ignoring parameter '{token}'; {name} reads {reads}");
    }
}

/// Parses the parameters of the plugin `name`, whose Ensembl VEP module reads
/// none: a `key=value` token is refused as the undeclared key every built-in
/// refuses, and every other token is logged as unread.
pub(crate) fn parse_none(name: &str, params: &[String]) -> Result<(), PluginError> {
    if let Some((key, _)) = params.iter().find_map(|param| perl_key_value(param)) {
        return Err(PluginError::Init(format!(
            "{name}: unknown parameter '{key}'; {name} takes no parameter"
        )));
    }
    warn_unread_parameters(name, params, "no parameter");
    Ok(())
}

/// How a plugin's Ensembl VEP module reads its parameters.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum Grammar<'a> {
    /// Every parameter `key=value`, or every one by position; `flags` are the
    /// tokens of the positional reading that are flags, not values.
    PositionalOrKeyValue { flags: &'a [&'a str] },
    /// Every parameter `key=value`; a positional list is refused as the data
    /// file given without its key, `file_key`.
    KeyValueOnly { file_key: &'a str },
    /// Positional tokens beside `key=value` pairs; `flags` as above.
    PositionalAndKeyValue { flags: &'a [&'a str] },
}

/// The parameters after a plugin's name, read by the plugin's grammar.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct PluginParams {
    positional: Vec<String>,
    keyed: Vec<(String, String)>,
    flags: Vec<String>,
    trailing_flag: Option<String>,
}

impl PluginParams {
    /// Parses the parameters of the plugin `name` as every parameter `key=value`
    /// or every one by position, with no flag token: the reading most Ensembl
    /// VEP plugins share. `declared` are the keys the plugin accepts.
    ///
    /// Every refusal is [`PluginError::Init`] with a message opening `<name>: `
    /// and naming the parameter at fault: a parameter without `=` beside one
    /// with it, a key given twice, or a key outside `declared`. Positional
    /// parameters are not checked against `declared`.
    pub fn parse(name: &str, declared: &[&str], params: &[String]) -> Result<Self, PluginError> {
        Self::parse_with(
            name,
            Grammar::PositionalOrKeyValue { flags: &[] },
            declared,
            params,
        )
    }

    /// Parses the parameters of the plugin `name` under `grammar`; `declared`
    /// are the keys the plugin accepts.
    ///
    /// The refusals of [`Self::parse`], plus a positional list under
    /// [`Grammar::KeyValueOnly`].
    pub fn parse_with(
        name: &str,
        grammar: Grammar<'_>,
        declared: &[&str],
        params: &[String],
    ) -> Result<Self, PluginError> {
        Self::parse_keys(name, grammar, Some(declared), params)
    }

    /// Parses the parameters of the plugin `name` under `grammar` with every
    /// key accepted, for a module that reads each key it does not name as a
    /// data file (CADD.pm:124-129): every refusal of [`Self::parse_with`]
    /// except the undeclared key's.
    pub fn parse_any_key(
        name: &str,
        grammar: Grammar<'_>,
        params: &[String],
    ) -> Result<Self, PluginError> {
        Self::parse_keys(name, grammar, None, params)
    }

    /// The reading [`Self::parse_with`] and [`Self::parse_any_key`] share;
    /// `declared` is `None` when every key is accepted.
    fn parse_keys(
        name: &str,
        grammar: Grammar<'_>,
        declared: Option<&[&str]>,
        params: &[String],
    ) -> Result<Self, PluginError> {
        let keyed_example = params.iter().find(|p| p.contains('='));
        // A grammar that reads the whole list one way puts every parameter
        // under `params_to_hash` once one carries `=`: a token without `=` is
        // then its parse failure, and a flag token is no flag (a bare `1`
        // beside `snv=...` fails). `keyed_reading` holds the parameter that
        // put the list under that reading; the grammar that reads positional
        // tokens beside pairs never has one.
        let (keyed_reading, flags) = match grammar {
            Grammar::PositionalAndKeyValue { flags } => (None, flags),
            Grammar::PositionalOrKeyValue { flags } => (keyed_example, flags),
            Grammar::KeyValueOnly { file_key } => {
                if keyed_example.is_none() && !params.is_empty() {
                    return Err(PluginError::Init(format!(
                        "{name}: no file specified; use {file_key}=<path>"
                    )));
                }
                (keyed_example, &[][..])
            }
        };
        let flags = if keyed_reading.is_none() { flags } else { &[] };

        let mut parsed = Self::default();
        for param in params {
            if flags.contains(&param.as_str()) {
                parsed.flags.push(param.clone());
                continue;
            }
            let Some((key, value)) = perl_key_value(param) else {
                if let Some(example) = keyed_reading {
                    let alternative = match grammar {
                        Grammar::KeyValueOnly { .. } => "",
                        _ => ", or every one by position",
                    };
                    return Err(PluginError::Init(format!(
                        "{name}: parameter '{param}' is not key=value while '{example}' is; \
                         give every parameter as key=value{alternative}"
                    )));
                }
                parsed.positional.push(param.clone());
                continue;
            };
            if let Some(declared) = declared {
                if !declared.contains(&key) {
                    let accepted = if declared.is_empty() {
                        format!("{name} reads its parameters by position")
                    } else {
                        format!("{name} accepts {}", declared.join(", "))
                    };
                    return Err(PluginError::Init(format!(
                        "{name}: unknown parameter '{key}'; {accepted}"
                    )));
                }
            }
            if parsed.keyed.iter().any(|(seen, _)| seen == key) {
                return Err(PluginError::Init(format!(
                    "{name}: parameter '{key}' is given twice"
                )));
            }
            parsed.keyed.push((key.to_string(), value.to_string()));
        }
        parsed.trailing_flag = params
            .last()
            .filter(|last| flags.contains(&last.as_str()))
            .cloned();
        Ok(parsed)
    }

    /// The value of `key`; `None` when the parameters do not name it.
    pub fn get(&self, key: &str) -> Option<&str> {
        self.keyed
            .iter()
            .find(|(k, _)| k == key)
            .map(|(_, v)| v.as_str())
    }

    /// The positional parameters in their written order, flags left out.
    pub fn positional(&self) -> &[String] {
        &self.positional
    }

    /// The `key=value` pairs in their written order.
    pub fn keyed(&self) -> &[(String, String)] {
        &self.keyed
    }

    /// Whether the flag `name` was given anywhere in the list.
    pub fn flag(&self, name: &str) -> bool {
        self.flags.iter().any(|flag| flag == name)
    }

    /// The flag that closes the list, when its last token is one.
    pub fn trailing_flag(&self) -> Option<&str> {
        self.trailing_flag.as_deref()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn params(tokens: &[&str]) -> Vec<String> {
        tokens.iter().map(|t| t.to_string()).collect()
    }

    fn pairs(items: &[(&str, &str)]) -> Vec<(String, String)> {
        items
            .iter()
            .map(|(k, v)| (k.to_string(), v.to_string()))
            .collect()
    }

    fn init_message<T: std::fmt::Debug>(result: Result<T, PluginError>) -> String {
        match result {
            Err(PluginError::Init(message)) => message,
            other => panic!("expected PluginError::Init, got {other:?}"),
        }
    }

    #[test]
    fn keyed_params_parse_in_order() {
        let parsed =
            PluginParams::parse("Example", &["cols", "file"], &params(&["file=a", "cols=b"]))
                .unwrap();
        assert_eq!(parsed.keyed(), pairs(&[("file", "a"), ("cols", "b")]));
        assert_eq!(parsed.get("file"), Some("a"));
        assert_eq!(parsed.get("cols"), Some("b"));
        assert_eq!(parsed.get("other"), None);
        assert!(parsed.positional().is_empty());
    }

    #[test]
    fn mixed_keyed_and_positional_is_an_init_error_naming_the_token() {
        let message = init_message(PluginParams::parse(
            "Example",
            &["file"],
            &params(&["file=a", "x"]),
        ));
        assert_eq!(
            message,
            "Example: parameter 'x' is not key=value while 'file=a' is; \
             give every parameter as key=value, or every one by position"
        );
    }

    #[test]
    fn duplicate_key_is_an_init_error() {
        let message = init_message(PluginParams::parse(
            "Example",
            &["file"],
            &params(&["file=a", "file=b"]),
        ));
        assert_eq!(message, "Example: parameter 'file' is given twice");
    }

    #[test]
    fn positional_params_keep_order() {
        let parsed = PluginParams::parse("Example", &["file"], &params(&["b", "a", "c"])).unwrap();
        assert_eq!(parsed.positional(), params(&["b", "a", "c"]).as_slice());
        assert_eq!(parsed.get("file"), None);
        assert!(parsed.keyed().is_empty());
        assert_eq!(parsed.trailing_flag(), None);
    }

    /// `my ($key, $val) = split('=', $param)` (Utils/BaseVepPlugin.pm:126,
    /// NearestExonJB.pm:101, dbNSFP.pm:589) assigns the first two fields and
    /// drops the rest: `gff3=0=x` runs NearestExonJB on the cache model (`0`),
    /// `intronic=1=0` reports both exons of an intron (`1`), `max_range=50=100`
    /// applies 50 bases, and `raw=1=0` writes AVI's raw score. A `==` leaves
    /// the second field empty.
    #[test]
    fn a_value_is_the_second_field_alone_as_ensembl_vep_reads_it() {
        for (token, value) in [
            ("file=a=b", "a"),
            ("file=a=b=c", "a"),
            ("file=a=", "a"),
            ("file==b", ""),
            ("file===", ""),
        ] {
            let parsed = PluginParams::parse("Example", &["file"], &params(&[token])).unwrap();
            assert_eq!(parsed.keyed(), pairs(&[("file", value)]), "{token}");
            assert_eq!(parsed.get("file"), Some(value), "{token}");
        }
        assert_eq!(perl_key_value("=b=c"), Some(("", "b")));
        assert_eq!(perl_key_value("="), Some(("", "")));
        assert_eq!(perl_key_value("file"), None);
        assert_eq!(perl_key_value(""), None);
    }

    /// Perl numifies a value's leading number and reads a value with none as
    /// 0: `01` is 1, `50.0` is 50, `1e1` is 10, `1abc` is 1, and ``, `abc`,
    /// `0x1A` and `len` are 0; `inf` and `nan` are themselves. The whole-number
    /// flag follows `looks_like_number`: blanks around a number, a trailing
    /// dot, a leading dot and `nanq` are whole; a second sign (`--1`, which
    /// `strtod` reads as 1), a blank after the sign (`- 5`, which it reads as
    /// -5), an exponent without digits, `0x1A`, text after the number and the
    /// empty value are not.
    #[test]
    fn perl_numeric_reads_the_leading_number() {
        for (value, number, whole) in [
            ("1", 1.0, true),
            ("01", 1.0, true),
            ("50.0", 50.0, true),
            ("1.", 1.0, true),
            (".5", 0.5, true),
            ("1e1", 10.0, true),
            ("1E+1", 10.0, true),
            (" 1 ", 1.0, true),
            ("1\x0b", 1.0, true),
            ("\x0b\x0c1", 1.0, true),
            ("+1", 1.0, true),
            ("-5", -5.0, true),
            ("0 but true", 0.0, true),
            ("--1", 1.0, false),
            ("+-1", -1.0, false),
            ("- 5", -5.0, false),
            ("1e", 1.0, false),
            ("1abc", 1.0, false),
            ("", 0.0, false),
            (" ", 0.0, false),
            (".", 0.0, false),
            ("abc", 0.0, false),
            ("0x1A", 0.0, false),
            ("len", 0.0, false),
            ("inf", f64::INFINITY, true),
            ("-Infinity", f64::NEG_INFINITY, true),
            ("1.#INF00", f64::INFINITY, true),
            ("1.#Infinity", f64::INFINITY, true),
            ("infx", f64::INFINITY, false),
            ("inf00", f64::INFINITY, false),
            ("1.#INF0x", f64::INFINITY, false),
            ("1.#INFI", f64::INFINITY, false),
            ("1.#5", 1.0, false),
            ("1.#infinity00", f64::INFINITY, false),
        ] {
            assert_eq!(perl_numeric("Example", "key", value), number, "{value:?}");
            assert_eq!(perl_whole_number(value), whole, "{value:?}");
        }
        for (value, whole) in [
            ("nan", true),
            ("NaNq", true),
            ("qnan", true),
            ("snanq", true),
            ("-SNaN", true),
            ("1.#QNAN(", false),
            ("nan(123)", true),
            ("nan(0x1F)", true),
            ("nan(0b101)", true),
            ("1.#IND00", true),
            ("1.#QNAN", true),
            ("1.#NAN", true),
            ("1.#nanq", true),
            ("1.#QNAN(123)", true),
            ("1.#indq", false),
            ("1.#nan00", false),
            ("nanx", false),
            ("nan()", false),
            ("nan(abc)", false),
            ("nan(0x)", false),
            ("1.#QNAN0", false),
        ] {
            assert!(perl_numeric("Example", "key", value).is_nan(), "{value:?}");
            assert_eq!(perl_whole_number(value), whole, "{value:?}");
        }
    }

    #[test]
    fn undeclared_key_is_an_init_error() {
        let message = init_message(PluginParams::parse(
            "Example",
            &["file", "cols"],
            &params(&["file=a", "colour=b"]),
        ));
        assert_eq!(
            message,
            "Example: unknown parameter 'colour'; Example accepts file, cols"
        );
    }

    /// `parse_any_key` keeps every key, as a module that opens each key's value
    /// as a data file reads them (CADD.pm:124-129), and keeps the other
    /// refusals: a key given twice and a bare token beside a `key=value` one.
    #[test]
    fn parse_any_key_accepts_every_key_and_keeps_the_other_refusals() {
        let grammar = Grammar::PositionalOrKeyValue { flags: &["1"] };
        let parsed =
            PluginParams::parse_any_key("Example", grammar, &params(&["file=a", "colour=b"]))
                .unwrap();
        assert_eq!(parsed.keyed(), pairs(&[("file", "a"), ("colour", "b")]));
        assert_eq!(
            init_message(PluginParams::parse_any_key(
                "Example",
                grammar,
                &params(&["file=a", "file=b"]),
            )),
            "Example: parameter 'file' is given twice"
        );
        assert_eq!(
            init_message(PluginParams::parse_any_key(
                "Example",
                grammar,
                &params(&["file=a", "1"]),
            )),
            "Example: parameter '1' is not key=value while 'file=a' is; \
             give every parameter as key=value, or every one by position"
        );
    }

    #[test]
    fn an_empty_value_is_kept_as_ensembl_vep_keeps_it() {
        // `my ($key, $val) = split('=', 'file=')` leaves `$val` the empty
        // string (Utils/BaseVepPlugin.pm:126), so `file=` is a keyed parameter
        // whose value is empty, not a parse failure.
        let parsed = PluginParams::parse("Example", &["file"], &params(&["file="])).unwrap();
        assert_eq!(parsed.keyed(), pairs(&[("file", "")]));
        assert_eq!(parsed.get("file"), Some(""));
    }

    #[test]
    fn no_params_are_positional_and_empty() {
        let parsed = PluginParams::parse("Example", &["file"], &[]).unwrap();
        assert_eq!(parsed, PluginParams::default());
        assert!(parsed.positional().is_empty());
    }

    /// A plugin whose module reads no parameter logs a positional token and
    /// refuses a `key=value` one, naming the key.
    #[test]
    fn parse_none_logs_positional_tokens_and_refuses_a_key() {
        assert!(parse_none("NMD", &[]).is_ok());
        assert!(parse_none("NMD", &params(&["a", "b"])).is_ok());
        assert_eq!(
            init_message(parse_none("NMD", &params(&["a", "file=x"]))),
            "NMD: unknown parameter 'file'; NMD takes no parameter"
        );
    }

    /// AlphaMissense reads `params_to_hash` alone and dies without `file`
    /// (AlphaMissense.pm:175-181), so a file given by position is the file
    /// given without its key.
    #[test]
    fn key_value_only_refuses_positional_file() {
        let grammar = Grammar::KeyValueOnly { file_key: "file" };
        let message = init_message(PluginParams::parse_with(
            "AlphaMissense",
            grammar,
            &["file", "cols", "transcript_match"],
            &params(&["/data/AlphaMissense_hg38.tsv.gz"]),
        ));
        assert_eq!(message, "AlphaMissense: no file specified; use file=<path>");
        let message = init_message(PluginParams::parse_with(
            "AlphaMissense",
            grammar,
            &["file", "cols"],
            &params(&["file=a", "x"]),
        ));
        assert_eq!(
            message,
            "AlphaMissense: parameter 'x' is not key=value while 'file=a' is; \
             give every parameter as key=value"
        );
        // An empty list is the empty hash; the plugin reports its missing file.
        let parsed = PluginParams::parse_with("AlphaMissense", grammar, &["file"], &[]).unwrap();
        assert_eq!(parsed, PluginParams::default());
        let parsed = PluginParams::parse_with(
            "SpliceAI",
            Grammar::KeyValueOnly { file_key: "snv" },
            &["snv", "indel", "cutoff"],
            &params(&["snv=a.synthetic.vcf.gz", "cutoff=0.5"]),
        )
        .unwrap();
        assert_eq!(parsed.get("cutoff"), Some("0.5"));
        assert_eq!(
            init_message(PluginParams::parse_with(
                "SpliceAI",
                Grammar::KeyValueOnly { file_key: "snv" },
                &["snv"],
                &params(&["a.synthetic.vcf.gz"]),
            )),
            "SpliceAI: no file specified; use snv=<path>"
        );
    }

    /// CADD's positional reading takes a trailing `1` as `force_annotate` and
    /// never as a file (CADD.pm:119-120); pLI's `sigfig` is a flag wherever it
    /// stands (pLI.pm:105-106). Under the keyed reading a bare `1` is the
    /// parse failure `params_to_hash` raises.
    #[test]
    fn positional_grammar_accepts_trailing_flag() {
        let cadd = Grammar::PositionalOrKeyValue { flags: &["1"] };
        let parsed = PluginParams::parse_with(
            "CADD",
            cadd,
            &["snv", "indels", "sv", "force_annotate"],
            &params(&["snv_synthetic.tsv.gz", "indels_synthetic.tsv.gz", "1"]),
        )
        .unwrap();
        assert_eq!(
            parsed.positional(),
            params(&["snv_synthetic.tsv.gz", "indels_synthetic.tsv.gz"]).as_slice()
        );
        assert!(parsed.flag("1"));
        assert_eq!(parsed.trailing_flag(), Some("1"));

        let parsed = PluginParams::parse_with(
            "CADD",
            cadd,
            &["snv"],
            &params(&["1", "snv_synthetic.tsv.gz"]),
        )
        .unwrap();
        assert_eq!(
            parsed.positional(),
            params(&["snv_synthetic.tsv.gz"]).as_slice()
        );
        assert!(parsed.flag("1"));
        assert_eq!(parsed.trailing_flag(), None);

        let message = init_message(PluginParams::parse_with(
            "CADD",
            cadd,
            &["snv"],
            &params(&["snv=a_synthetic.tsv.gz", "1"]),
        ));
        assert_eq!(
            message,
            "CADD: parameter '1' is not key=value while 'snv=a_synthetic.tsv.gz' is; \
             give every parameter as key=value, or every one by position"
        );

        let pli = Grammar::PositionalOrKeyValue { flags: &["sigfig"] };
        let parsed = PluginParams::parse_with(
            "pLI",
            pli,
            &[],
            &params(&["pLI_values.txt", "sigfig", "transcript"]),
        )
        .unwrap();
        assert_eq!(
            parsed.positional(),
            params(&["pLI_values.txt", "transcript"]).as_slice()
        );
        assert!(parsed.flag("sigfig"));
        assert_eq!(parsed.trailing_flag(), None);
        assert!(
            !PluginParams::parse_with("pLI", pli, &[], &params(&["pLI_values.txt"]))
                .unwrap()
                .flag("sigfig")
        );
    }

    /// dbNSFP reads positional tokens and `name=value` pairs side by side and
    /// dies on a name it does not know (dbNSFP.pm:192, 586, 601); a plugin that
    /// reads by position alone names that in its refusal of any key.
    #[test]
    fn unknown_named_parameter_is_an_error() {
        let dbnsfp = Grammar::PositionalAndKeyValue { flags: &[] };
        let declared = &[
            "consequence",
            "output_prefix",
            "pep_match",
            "transcript_match",
        ];
        let parsed = PluginParams::parse_with(
            "dbNSFP",
            dbnsfp,
            declared,
            &params(&[
                "dbNSFP_synthetic.gz",
                "consequence=ALL",
                "SIFT_score",
                "pep_match=0",
            ]),
        )
        .unwrap();
        assert_eq!(
            parsed.positional(),
            params(&["dbNSFP_synthetic.gz", "SIFT_score"]).as_slice()
        );
        assert_eq!(
            parsed.keyed(),
            pairs(&[("consequence", "ALL"), ("pep_match", "0")])
        );
        let message = init_message(PluginParams::parse_with(
            "dbNSFP",
            dbnsfp,
            declared,
            &params(&["dbNSFP_synthetic.gz", "transcript=1"]),
        ));
        assert_eq!(
            message,
            "dbNSFP: unknown parameter 'transcript'; dbNSFP accepts consequence, output_prefix, pep_match, transcript_match"
        );
        let message = init_message(PluginParams::parse_with(
            "pLI",
            Grammar::PositionalOrKeyValue { flags: &["sigfig"] },
            &[],
            &params(&["file=pLI_values.txt"]),
        ));
        assert_eq!(
            message,
            "pLI: unknown parameter 'file'; pLI reads its parameters by position"
        );
    }
}
