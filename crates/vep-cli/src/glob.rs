// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The expansion Ensembl VEP applies to its plugins directory, `glob($plugins_dir)`
//! (Runner.pm:908): Perl's `File::Glob::bsd_glob` under the `glob` operator's
//! default flags, `GLOB_BRACE | GLOB_NOMAGIC | GLOB_QUOTE | GLOB_TILDE |
//! GLOB_ALPHASORT`. The steps run in that function's order: a backslash quotes
//! the byte after it, the first brace group splits the pattern into alternatives,
//! and each alternative has its leading `~` replaced, is compiled, is walked
//! against the file system a segment at a time and has its matches sorted.
//!
//! The `glob` operator's splitting of its argument on whitespace
//! (`File::Glob::csh_glob`) is not applied: the value is one pattern, spaces
//! included.

use std::fs;

/// Every path `pattern` expands to, in Ensembl VEP's order: the brace alternatives
/// as written, each alternative's matches sorted case-insensitively with byte
/// order breaking ties (`GLOB_ALPHASORT`), and an alternative holding none of
/// `*`, `?` and `[` standing as itself whether or not it exists (`GLOB_NOMAGIC`).
/// `home` is the directory `~` stands for; `~user` is the user's home directory
/// from the password database, left as written for a user it does not name.
pub(crate) fn expand(pattern: &str, home: &str) -> Vec<String> {
    let mut paths = Vec::new();
    for alternative in brace_alternatives(&quote(pattern.as_bytes())) {
        let alternative = expand_tilde(alternative, home);
        let (tokens, magic) = compile(&alternative);
        if !magic {
            paths.push(literal(&tokens));
            continue;
        }
        let mut matches = Vec::new();
        walk(&mut String::new(), &tokens, &mut matches);
        matches.sort_by(|a, b| alpha_order(a.as_bytes(), b.as_bytes()));
        paths.extend(matches);
    }
    paths
}

/// A pattern byte and whether a backslash quoted it (`GLOB_QUOTE`), which keeps a
/// metacharacter plain; a trailing backslash is itself.
#[derive(Clone, Copy)]
struct Quoted {
    byte: u8,
    quoted: bool,
}

impl Quoted {
    /// Whether this is the unquoted metacharacter `byte`.
    fn plain(self, byte: u8) -> bool {
        !self.quoted && self.byte == byte
    }
}

fn quote(pattern: &[u8]) -> Vec<Quoted> {
    let mut out = Vec::with_capacity(pattern.len());
    let mut bytes = pattern.iter();
    while let Some(&byte) = bytes.next() {
        if byte != b'\\' {
            out.push(Quoted {
                byte,
                quoted: false,
            });
            continue;
        }
        out.push(Quoted {
            byte: bytes.next().map_or(b'\\', |&next| next),
            quoted: true,
        });
    }
    out
}

/// The alternatives of the first brace group (`GLOB_BRACE`): the text before the
/// `{`, one comma-separated member, the text after the matching `}`, each
/// expanded again for the groups it holds, in written order. `{}` alone stands
/// as itself; a `{` with no matching `}` leaves the text before it as the one
/// alternative; past the opening `{`, a `[` is skipped to the first `]` after
/// it, so a brace or comma between the two does not count.
fn brace_alternatives(pattern: &[Quoted]) -> Vec<Vec<Quoted>> {
    if pattern.len() == 2 && pattern[0].plain(b'{') && pattern[1].plain(b'}') {
        return vec![pattern.to_vec()];
    }
    let Some(open) = pattern.iter().position(|q| q.plain(b'{')) else {
        return vec![pattern.to_vec()];
    };
    let Some(close) = matching_brace(pattern, open) else {
        return vec![pattern[..open].to_vec()];
    };
    let mut alternatives = Vec::new();
    let mut depth = 0usize;
    let mut start = open + 1;
    let mut i = open + 1;
    while i <= close {
        let q = pattern[i];
        if q.plain(b'[') {
            i = class_end(pattern, i);
        } else if q.plain(b'{') {
            depth += 1;
        } else if q.plain(b'}') && depth > 0 {
            depth -= 1;
        } else if i == close || (q.plain(b',') && depth == 0) {
            let mut member = pattern[..open].to_vec();
            member.extend_from_slice(&pattern[start..i]);
            member.extend_from_slice(&pattern[close + 1..]);
            alternatives.extend(brace_alternatives(&member));
            start = i + 1;
        }
        i += 1;
    }
    alternatives
}

/// The index of the `}` balancing the `{` at `open`, if one exists.
fn matching_brace(pattern: &[Quoted], open: usize) -> Option<usize> {
    let mut depth = 0usize;
    let mut i = open + 1;
    while i < pattern.len() {
        let q = pattern[i];
        if q.plain(b'[') {
            i = class_end(pattern, i);
        } else if q.plain(b'{') {
            depth += 1;
        } else if q.plain(b'}') {
            if depth == 0 {
                return Some(i);
            }
            depth -= 1;
        }
        i += 1;
    }
    None
}

/// The index of the first `]` after the `[` at `open`, or `open` itself when
/// none follows, so the brace scans step over the `[` alone.
fn class_end(pattern: &[Quoted], open: usize) -> usize {
    pattern[open + 1..]
        .iter()
        .position(|q| q.plain(b']'))
        .map_or(open, |offset| open + 1 + offset)
}

/// A leading `~` or `~user` replaced by the home directory (`GLOB_TILDE`): `home`
/// for `~`, the password database's directory for `~user`, the alternative as it
/// is when the database names no such user. An empty `home` makes `~` the root
/// and `~/x` `/x`, as `expand_tilde_text` reads them.
fn expand_tilde(alternative: Vec<Quoted>, home: &str) -> Vec<Quoted> {
    if !alternative.first().is_some_and(|q| q.plain(b'~')) {
        return alternative;
    }
    let name_end = alternative
        .iter()
        .position(|q| q.plain(b'/'))
        .unwrap_or(alternative.len());
    let user: Vec<u8> = alternative[1..name_end].iter().map(|q| q.byte).collect();
    let dir = if user.is_empty() {
        if home.is_empty() && name_end == alternative.len() {
            "/".to_string()
        } else {
            home.to_string()
        }
    } else {
        match std::str::from_utf8(&user).ok().and_then(home_of_user) {
            Some(dir) => dir,
            None => return alternative,
        }
    };
    let mut out: Vec<Quoted> = dir
        .bytes()
        .map(|byte| Quoted {
            byte,
            quoted: false,
        })
        .collect();
    out.extend_from_slice(&alternative[name_end..]);
    out
}

/// The home directory of `user` in the password database, `None` when it names
/// no such user.
#[cfg(unix)]
fn home_of_user(user: &str) -> Option<String> {
    use std::ffi::{CStr, CString};

    let name = CString::new(user).ok()?;
    let mut buffer = vec![0u8; 16 * 1024];
    loop {
        // SAFETY: `passwd` holds integers and pointers alone, for which all-zero
        // bytes are a valid value.
        let mut entry: libc::passwd = unsafe { std::mem::zeroed() };
        let mut found: *mut libc::passwd = std::ptr::null_mut();
        // SAFETY: `name` is NUL-terminated, `entry` and `found` are live locals
        // the call writes to, and `buffer` is valid for `buffer.len()` bytes;
        // on success the entry's string fields point into `buffer`, which stays
        // alive until `pw_dir` has been copied out below.
        let status = unsafe {
            libc::getpwnam_r(
                name.as_ptr(),
                &mut entry,
                buffer.as_mut_ptr().cast(),
                buffer.len(),
                &mut found,
            )
        };
        if status == libc::ERANGE {
            buffer.resize(buffer.len() * 2, 0);
            continue;
        }
        if status != 0 || found.is_null() || entry.pw_dir.is_null() {
            return None;
        }
        // SAFETY: `pw_dir` is a NUL-terminated string inside `buffer`.
        let dir = unsafe { CStr::from_ptr(entry.pw_dir) };
        return dir.to_str().ok().map(str::to_owned);
    }
}

#[cfg(not(unix))]
fn home_of_user(_user: &str) -> Option<String> {
    None
}

/// One compiled element of an alternative.
#[derive(PartialEq, Eq)]
enum Token {
    /// One byte as itself, a quoted metacharacter included; `/` separates
    /// segments.
    Byte(u8),
    /// `*`: any run of bytes. Adjacent stars compile to one.
    Any,
    /// `?`: one byte.
    One,
    /// `[...]`: one byte inside one of the ranges, or outside all of them when
    /// negated by a leading `!`.
    Set { negate: bool, ranges: Vec<(u8, u8)> },
}

/// The alternative compiled, and whether it holds a metacharacter
/// (`GLOB_MAGCHAR`).
fn compile(alternative: &[Quoted]) -> (Vec<Token>, bool) {
    let mut tokens = Vec::with_capacity(alternative.len());
    let mut magic = false;
    let mut i = 0;
    while i < alternative.len() {
        let q = alternative[i];
        i += 1;
        if q.plain(b'[') {
            if let Some((set, next)) = compile_set(alternative, i) {
                tokens.push(set);
                magic = true;
                i = next;
                continue;
            }
            tokens.push(Token::Byte(b'['));
        } else if q.plain(b'?') {
            tokens.push(Token::One);
            magic = true;
        } else if q.plain(b'*') {
            magic = true;
            if tokens.last() != Some(&Token::Any) {
                tokens.push(Token::Any);
            }
        } else {
            tokens.push(Token::Byte(q.byte));
        }
    }
    (tokens, magic)
}

/// The set whose `[` precedes `start`, with the index after its `]`. The first
/// member, after an optional `!`, is taken as it is, so `[]]` holds `]`; a `-`
/// between two members that is not followed by the `]` makes a range; `None`
/// when no `]` follows the first member, which leaves the `[` a plain byte.
fn compile_set(alternative: &[Quoted], start: usize) -> Option<(Token, usize)> {
    let mut i = start;
    let negate = alternative.get(i).is_some_and(|q| q.plain(b'!'));
    if negate {
        i += 1;
    }
    let mut member = *alternative.get(i)?;
    if !alternative[i + 1..].iter().any(|q| q.plain(b']')) {
        return None;
    }
    i += 1;
    let mut ranges = Vec::new();
    loop {
        let range_end = alternative
            .get(i)
            .filter(|q| q.plain(b'-'))
            .and_then(|_| alternative.get(i + 1))
            .filter(|q| !q.plain(b']'));
        match range_end {
            Some(end) => {
                ranges.push((member.byte, end.byte));
                i += 2;
            }
            None => ranges.push((member.byte, member.byte)),
        }
        let next = *alternative.get(i)?;
        i += 1;
        if next.plain(b']') {
            break;
        }
        member = next;
    }
    Some((Token::Set { negate, ranges }, i))
}

/// The text of an alternative without a metacharacter, every byte as itself.
fn literal(tokens: &[Token]) -> String {
    let bytes: Vec<u8> = tokens
        .iter()
        .map(|token| match token {
            Token::Byte(byte) => *byte,
            Token::Any | Token::One | Token::Set { .. } => {
                unreachable!("a literal has no metacharacter")
            }
        })
        .collect();
    String::from_utf8_lossy(&bytes).into_owned()
}

/// The existing paths matching `tokens` under `path`, unsorted: a segment without
/// a metacharacter is appended as written, with every `/` that follows it; the
/// first segment with one is matched against the entries of the directory built
/// so far; a pattern consumed whole names a path when `lstat` finds it.
fn walk(path: &mut String, mut tokens: &[Token], matches: &mut Vec<String>) {
    loop {
        if tokens.is_empty() {
            if fs::symlink_metadata(path.as_str()).is_ok() {
                matches.push(path.clone());
            }
            return;
        }
        let segment_end = tokens
            .iter()
            .position(|token| *token == Token::Byte(b'/'))
            .unwrap_or(tokens.len());
        let (segment, rest) = tokens.split_at(segment_end);
        if segment.iter().any(|token| !matches!(token, Token::Byte(_))) {
            match_segment(path, segment, rest, matches);
            return;
        }
        path.push_str(&literal(segment));
        let separators = rest
            .iter()
            .take_while(|token| **token == Token::Byte(b'/'))
            .count();
        path.extend(std::iter::repeat_n('/', separators));
        tokens = &rest[separators..];
    }
}

/// Each entry of the directory `path` names (the working directory for an empty
/// `path`) that `segment` matches continues the walk over `rest`. The entries
/// are the directory's own, `.` and `..` among them; one whose name opens with
/// `.` is matched only by a segment opening with a plain `.`. A directory that
/// cannot be read matches nothing.
fn match_segment(path: &mut String, segment: &[Token], rest: &[Token], matches: &mut Vec<String>) {
    let directory = if path.is_empty() { "." } else { path.as_str() };
    let Ok(entries) = fs::read_dir(directory) else {
        return;
    };
    let mut names = vec![".".to_string(), "..".to_string()];
    names.extend(
        entries
            .filter_map(Result::ok)
            .filter_map(|entry| entry.file_name().into_string().ok()),
    );
    let dot_first = segment.first() == Some(&Token::Byte(b'.'));
    let prefix_len = path.len();
    for name in names {
        if name.starts_with('.') && !dot_first {
            continue;
        }
        if !matches_name(name.as_bytes(), segment) {
            continue;
        }
        path.push_str(&name);
        walk(path, rest, matches);
        path.truncate(prefix_len);
    }
}

/// Whether `segment` matches the whole of `name`.
fn matches_name(mut name: &[u8], mut segment: &[Token]) -> bool {
    while let Some((token, rest)) = segment.split_first() {
        match token {
            Token::Any => {
                if rest.is_empty() {
                    return true;
                }
                loop {
                    if matches_name(name, rest) {
                        return true;
                    }
                    let Some((_, shorter)) = name.split_first() else {
                        return false;
                    };
                    name = shorter;
                }
            }
            Token::One => {
                let Some((_, shorter)) = name.split_first() else {
                    return false;
                };
                name = shorter;
            }
            Token::Set { negate, ranges } => {
                let Some((&byte, shorter)) = name.split_first() else {
                    return false;
                };
                let inside = ranges
                    .iter()
                    .any(|&(low, high)| low <= byte && byte <= high);
                if inside == *negate {
                    return false;
                }
                name = shorter;
            }
            Token::Byte(expected) => {
                let Some((&byte, shorter)) = name.split_first() else {
                    return false;
                };
                if byte != *expected {
                    return false;
                }
                name = shorter;
            }
        }
        segment = rest;
    }
    name.is_empty()
}

/// `GLOB_ALPHASORT`'s order: the first byte pair that differs once ASCII case is
/// folded decides, a path's end counting as the byte 0, and paths equal under the
/// fold in byte order. The deciding pair is compared as the target's C `char`,
/// the type the reference subtracts, so a byte outside ASCII sorts before every
/// ASCII byte where `char` is signed and after them where it is unsigned.
fn alpha_order(a: &[u8], b: &[u8]) -> std::cmp::Ordering {
    let differ = a
        .iter()
        .zip(b)
        .position(|(x, y)| !x.eq_ignore_ascii_case(y))
        .unwrap_or(a.len().min(b.len()));
    let fold = |bytes: &[u8]| {
        bytes
            .get(differ)
            .map_or(0, |byte| as_c_char(byte.to_ascii_lowercase()))
    };
    fold(a).cmp(&fold(b)).then_with(|| a.cmp(b))
}

/// `byte` read as the target's C `char`: negative from `0x80` up where that type
/// is signed, the byte's own value where it is unsigned.
fn as_c_char(byte: u8) -> i32 {
    let signed = i32::from(std::ffi::c_char::MIN) < 0;
    if signed && byte >= 0x80 {
        i32::from(byte) - 256
    } else {
        i32::from(byte)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A temporary tree holding the directories `p1`, `p2`, `a`, `B`, `A1`, `a1`,
    /// `Ab`, `aB`, `.hidden` and `star*dir`, and the file `f.txt`; `expand` on it
    /// returns the matches with the tree's path replaced by `T`.
    struct Tree {
        dir: tempfile::TempDir,
    }

    impl Tree {
        fn new() -> Self {
            let dir = tempfile::tempdir().unwrap();
            for name in [
                "p1", "p2", "a", "B", ".hidden", "star*dir", "A1", "a1", "Ab", "aB",
            ] {
                fs::create_dir(dir.path().join(name)).unwrap();
            }
            fs::write(dir.path().join("f.txt"), "").unwrap();
            Self { dir }
        }

        fn root(&self) -> &str {
            self.dir.path().to_str().unwrap()
        }

        fn expand(&self, pattern: &str) -> Vec<String> {
            let root = self.root();
            expand(&format!("{root}/{pattern}"), "/home/someone")
                .into_iter()
                .map(|path| path.replacen(root, "T", 1))
                .collect()
        }
    }

    fn strings(items: &[&str]) -> Vec<String> {
        items.iter().map(|item| item.to_string()).collect()
    }

    /// `*` matches every entry but the dot-prefixed ones, sorted with ASCII case
    /// folded and byte order breaking the ties (`a` before `A1`, `A1` before `a1`,
    /// `Ab` before `aB`); `?` is one byte; `**` is `*`; `*.txt` reaches the file.
    #[test]
    fn star_and_question_match_as_bsd_glob_does() {
        let tree = Tree::new();
        assert_eq!(
            tree.expand("*"),
            strings(&[
                "T/a",
                "T/A1",
                "T/a1",
                "T/Ab",
                "T/aB",
                "T/B",
                "T/f.txt",
                "T/p1",
                "T/p2",
                "T/star*dir"
            ])
        );
        assert_eq!(tree.expand("**"), tree.expand("*"));
        assert_eq!(tree.expand("p?"), strings(&["T/p1", "T/p2"]));
        assert_eq!(tree.expand("?1"), strings(&["T/A1", "T/a1", "T/p1"]));
        assert_eq!(tree.expand("*.txt"), strings(&["T/f.txt"]));
    }

    /// A name whose first byte is outside ASCII sorts as the target's C `char`
    /// reads that byte: after every ASCII name where the type is unsigned, before
    /// them where it is signed (`é` opens with the byte `0xC3`).
    #[test]
    fn a_byte_outside_ascii_sorts_as_the_platform_char() {
        let tree = Tree::new();
        fs::create_dir(tree.dir.path().join("é1")).unwrap();
        let expected = if i32::from(std::ffi::c_char::MIN) < 0 {
            ["T/é1", "T/A1", "T/a1", "T/p1"]
        } else {
            ["T/A1", "T/a1", "T/p1", "T/é1"]
        };
        assert_eq!(tree.expand("*1"), strings(&expected));
        assert_eq!(tree.expand("{p,é}1"), strings(&["T/p1", "T/é1"]));
    }

    /// A pattern with no metacharacter is itself whether or not it exists; one
    /// with a metacharacter and no match is nothing; a `[` no `]` closes is a
    /// plain byte, and so is the pattern after an unmatched `}`.
    #[test]
    fn a_pattern_without_a_metacharacter_stands_as_itself() {
        let tree = Tree::new();
        assert_eq!(tree.expand("p1"), strings(&["T/p1"]));
        assert_eq!(tree.expand("nope"), strings(&["T/nope"]));
        assert_eq!(tree.expand("f.txt/"), strings(&["T/f.txt/"]));
        assert_eq!(tree.expand("zz*"), Vec::<String>::new());
        assert_eq!(tree.expand("p*/nope"), Vec::<String>::new());
        assert_eq!(tree.expand("[p"), strings(&["T/[p"]));
        assert_eq!(tree.expand("p1}"), strings(&["T/p1}"]));
    }

    /// Brace members expand in written order and are never sorted across each
    /// other; a member with no metacharacter stands as itself; `{}` alone is
    /// itself; an unmatched `{` leaves the text before it; a member may be empty
    /// or hold a `/`; a `{` inside a class still opens a group, here unmatched.
    #[test]
    fn braces_expand_in_written_order() {
        let tree = Tree::new();
        assert_eq!(
            tree.expand("{p2,p1,nope}"),
            strings(&["T/p2", "T/p1", "T/nope"])
        );
        assert_eq!(tree.expand("{p,B}*"), strings(&["T/p1", "T/p2", "T/B"]));
        assert_eq!(tree.expand("{x,y}"), strings(&["T/x", "T/y"]));
        assert_eq!(
            tree.expand("{p*,nope}"),
            strings(&["T/p1", "T/p2", "T/nope"])
        );
        assert_eq!(
            tree.expand("{A1,Ab,aB,a1}"),
            strings(&["T/A1", "T/Ab", "T/aB", "T/a1"])
        );
        assert_eq!(tree.expand("{p{1,2},a}"), strings(&["T/p1", "T/p2", "T/a"]));
        assert_eq!(tree.expand("{,p1}"), strings(&["T/", "T/p1"]));
        assert_eq!(tree.expand("{p1/,a}"), strings(&["T/p1/", "T/a"]));
        assert_eq!(tree.expand("p1{}"), strings(&["T/p1"]));
        assert_eq!(expand("{}", "/home/someone"), strings(&["{}"]));
        assert_eq!(tree.expand("p{1"), strings(&["T/p"]));
        assert_eq!(tree.expand("[{]p1"), strings(&["T/["]));
    }

    /// A class matches one byte of its members and ranges, `!` negates it, the
    /// first member may be `]`, and a `-` before the `]` is a member.
    #[test]
    fn classes_match_one_byte() {
        let tree = Tree::new();
        assert_eq!(tree.expand("[ab]*"), strings(&["T/a", "T/a1", "T/aB"]));
        assert_eq!(tree.expand("[a-b]1"), strings(&["T/a1"]));
        assert_eq!(tree.expand("[!pfaAB]*"), strings(&["T/star*dir"]));
        assert_eq!(tree.expand("[!]]1"), strings(&["T/A1", "T/a1", "T/p1"]));
        assert_eq!(tree.expand("[a-]*"), strings(&["T/a", "T/a1", "T/aB"]));
        fs::create_dir(tree.dir.path().join("]b")).unwrap();
        assert_eq!(tree.expand("[]]b"), strings(&["T/]b"]));
    }

    /// A backslash quotes the byte after it: `\*` and `[*]` both match a literal
    /// star, `\{` opens no group, `\~` is no home directory, and a trailing
    /// backslash is itself.
    #[test]
    fn a_backslash_quotes_the_next_byte() {
        let tree = Tree::new();
        assert_eq!(tree.expand("star\\*dir"), strings(&["T/star*dir"]));
        assert_eq!(tree.expand("star[*]dir"), strings(&["T/star*dir"]));
        assert_eq!(tree.expand("\\{p1,p2}"), strings(&["T/{p1,p2}"]));
        assert_eq!(tree.expand("p1\\"), strings(&["T/p1\\"]));
        assert_eq!(expand("\\~root", "/home/someone"), strings(&["~root"]));
    }

    /// An entry opening with `.` is matched only by a segment opening with `.`,
    /// where `.` and `..` are entries too; a class opening the segment does not
    /// count.
    #[test]
    fn dot_entries_need_a_leading_dot() {
        let tree = Tree::new();
        assert_eq!(tree.expand(".*"), strings(&["T/.", "T/..", "T/.hidden"]));
        assert_eq!(tree.expand("\\.h*"), strings(&["T/.hidden"]));
        assert_eq!(tree.expand("[.]*"), Vec::<String>::new());
        assert!(!tree.expand("*").iter().any(|path| path.contains("/.")));
    }

    /// A trailing `/` keeps only the directories and stays in the match; repeated
    /// separators stay as written; a relative pattern matches under the working
    /// directory; a symbolic link matches whatever it points to.
    #[test]
    #[cfg(unix)]
    fn separators_and_links_follow_lstat() {
        let tree = Tree::new();
        std::os::unix::fs::symlink(tree.dir.path().join("p1"), tree.dir.path().join("link"))
            .unwrap();
        std::os::unix::fs::symlink(
            tree.dir.path().join("nowhere"),
            tree.dir.path().join("dangling"),
        )
        .unwrap();
        assert_eq!(tree.expand("p*/"), strings(&["T/p1/", "T/p2/"]));
        assert_eq!(tree.expand("f*/"), Vec::<String>::new());
        assert_eq!(tree.expand("/p*"), strings(&["T//p1", "T//p2"]));
        assert_eq!(tree.expand("[dl]*"), strings(&["T/dangling", "T/link"]));
        let root = tree.root();
        let relative: Vec<String> = expand("p*", "/home/someone");
        assert!(
            !relative.iter().any(|path| path.starts_with(root)),
            "{relative:?}"
        );
        assert_eq!(
            expand(&format!("{root}/*/sub"), "/home/someone"),
            Vec::<String>::new()
        );
        fs::create_dir(tree.dir.path().join("p1").join("sub")).unwrap();
        assert_eq!(tree.expand("*/sub"), strings(&["T/link/sub", "T/p1/sub"]));
    }

    /// `~` is `home`, `~/x` is under it, an empty `home` makes `~` the root and
    /// `~/x` `/x`; `~user` is the user's directory in the password database
    /// (`root`'s on every supported system) and stays as written for a user it
    /// does not name; a member of a brace group has its own `~`.
    #[test]
    fn tilde_names_a_home_directory() {
        assert_eq!(expand("~", "/home/someone"), strings(&["/home/someone"]));
        assert_eq!(
            expand("~/x", "/home/someone"),
            strings(&["/home/someone/x"])
        );
        assert_eq!(expand("~", ""), strings(&["/"]));
        assert_eq!(expand("~/x", ""), strings(&["/x"]));
        assert_eq!(
            expand("~no_such_user_here/x", "/home/someone"),
            strings(&["~no_such_user_here/x"])
        );
        let root_home = if cfg!(target_os = "macos") {
            "/var/root"
        } else {
            "/root"
        };
        assert_eq!(
            expand("~root/x", "/home/someone"),
            strings(&[&format!("{root_home}/x")])
        );
        assert_eq!(expand("~root", "/home/someone"), strings(&[root_home]));
        assert_eq!(
            expand("{~/a,~root/b}", "/home/someone"),
            strings(&["/home/someone/a", &format!("{root_home}/b")])
        );
    }
}
