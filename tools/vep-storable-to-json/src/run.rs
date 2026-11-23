// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The run: the conversion's steps in their order, with the shards converted on the
//! host's cores and every file and line written in that order.

use std::collections::{BTreeMap, BTreeSet};
use std::ffi::OsStr;
use std::io::Write;
use std::os::unix::ffi::OsStrExt;
use std::path::Path;

use crate::info::{cache_version, contig_order, source_type, Info};
use crate::json::{print_bytes, write_pretty, write_string, PerlStr, Pretty};
use crate::perl::Die;
use crate::pool::run_ordered;
use crate::regions::{self, Noted, SequenceRegion};
use crate::regulatory::{self, CellTypes, RegulatoryShard};
use crate::shard::load_shard;
use crate::transcripts::{self, Walk};

/// What the command line gave.
pub struct Request {
    pub cache_dir: Vec<u8>,
    pub output_dir: Vec<u8>,
    /// The chromosomes of the third argument; empty for every chromosome.
    pub chromosomes: Vec<Vec<u8>>,
    pub threads: usize,
}

/// The run stopped: the exit code, after its line went to standard error.
struct Stop(i32);

fn path(bytes: &[u8]) -> &Path {
    Path::new(OsStr::from_bytes(bytes))
}

/// `"$left/$right"`.
fn join(left: &[u8], right: &[u8]) -> Vec<u8> {
    [left, b"/", right].concat()
}

/// The operating system's text for an error, as Perl's `$!` gives it.
fn os_text(error: &std::io::Error) -> String {
    let text = error.to_string();
    match text.rfind(" (os error ") {
        Some(at) if text.ends_with(')') => text[..at].to_string(),
        _ => text,
    }
}

/// The exit code of a Perl `die` outside `eval`: the error number when there is one.
fn die_code(error: &std::io::Error) -> i32 {
    match error.raw_os_error() {
        Some(code) if code & 0xff != 0 => code & 0xff,
        _ => 255,
    }
}

/// A name of `<digits>-<digits><suffix>`, with Perl's `$` allowing one final newline.
fn shard_stem<'a>(name: &'a [u8], suffix: &[u8]) -> Option<&'a [u8]> {
    let body = name.strip_suffix(b"\n").unwrap_or(name);
    let stem = body.strip_suffix(suffix)?;
    let dash = stem.iter().position(|byte| *byte == b'-')?;
    let digits = |part: &[u8]| !part.is_empty() && part.iter().all(u8::is_ascii_digit);
    (digits(&stem[..dash]) && digits(&stem[dash + 1..])).then_some(stem)
}

/// The names of a directory that match `<digits>-<digits><suffix>`, sorted.
fn shard_names(entries: &[Vec<u8>], suffix: &[u8]) -> Vec<Vec<u8>> {
    let mut names: Vec<Vec<u8>> = entries
        .iter()
        .filter(|name| shard_stem(name, suffix).is_some())
        .cloned()
        .collect();
    names.sort();
    names
}

/// `s/<suffix>$/.json/` on a name that matched.
fn json_name(name: &[u8], suffix: &[u8]) -> Vec<u8> {
    let (body, newline): (&[u8], &[u8]) = match name.strip_suffix(b"\n") {
        Some(body) => (body, b"\n"),
        None => (name, b""),
    };
    let stem = body.strip_suffix(suffix).unwrap_or(body);
    [stem, b".json", newline].concat()
}

fn read_dir_names(directory: &[u8]) -> std::io::Result<Vec<Vec<u8>>> {
    let mut names = Vec::new();
    for entry in std::fs::read_dir(path(directory))? {
        names.push(entry?.file_name().as_bytes().to_vec());
    }
    Ok(names)
}

/// `File::Copy::copy`: the bytes, into a file made with the default permissions.
fn copy_file(from: &[u8], to: &[u8]) -> std::io::Result<()> {
    let mut source = std::fs::File::open(path(from))?;
    let mut target = std::fs::File::create(path(to))?;
    std::io::copy(&mut source, &mut target)?;
    Ok(())
}

enum Task {
    Transcripts(Vec<u8>),
    Regulatory(Vec<u8>),
}

enum Outcome {
    /// The walk over a shard that was read, or why it could not be read.
    Transcripts(Result<Walk, Die>),
    Regulatory(RegulatoryShard),
}

struct Run<'a, O: Write, E: Write> {
    request: &'a Request,
    out: &'a mut O,
    err: &'a mut E,
}

impl<O: Write, E: Write> Run<'_, O, E> {
    fn say(&mut self, parts: &[&[u8]]) {
        for part in parts {
            self.out.write_all(part).ok();
        }
    }

    fn warn(&mut self, parts: &[&[u8]]) {
        // Standard error is not buffered: a line is written when it is made.
        self.err.write_all(&parts.concat()).ok();
        self.err.flush().ok();
    }

    /// A fatal error: one `ERROR:` line and exit code 1. For a regulatory shard that cannot be
    /// converted `parts` end in a line end, which leaves an empty line after it.
    fn fail(&mut self, parts: &[&[u8]]) -> Stop {
        self.warn(&[b"ERROR: [storable_to_json] ", &parts.concat(), b"\n"]);
        Stop(1)
    }

    /// A message of a Perl `die`, as the bytes Perl prints for it.
    fn die_text(error: &Die) -> Vec<u8> {
        print_bytes(&error.0)
    }

    /// `File::Path::make_path`, which dies on failure.
    fn make_path(&mut self, directory: &[u8]) -> Result<(), Stop> {
        std::fs::create_dir_all(path(directory)).map_err(|error| {
            self.warn(&[
                b"mkdir ",
                directory,
                b": ",
                os_text(&error).as_bytes(),
                b"\n",
            ]);
            Stop(die_code(&error))
        })
    }

    fn run(&mut self) -> Result<(), Stop> {
        let request = self.request;
        let cache_dir = &request.cache_dir[..];
        let output_dir = &request.output_dir[..];

        let info_txt = join(cache_dir, b"info.txt");
        let mut info = match std::fs::read(path(&info_txt)) {
            Ok(text) => Info::from_lines(&text),
            Err(_) => {
                self.warn(&[
                    b"WARN: no ",
                    &info_txt,
                    b"; info.json carries no source versions\n",
                ]);
                Info::default()
            }
        };
        if let Some(version) = cache_version(cache_dir) {
            info.top.insert(b"cache_version".to_vec(), version);
        }
        if let Some(set) = source_type(cache_dir) {
            info.top
                .insert(b"source_type".to_vec(), Pretty::String(set.to_vec()));
        }
        if info.string(b"serialiser_type") == Some(b"sereal") {
            return Err(self.fail(&[
                cache_dir,
                b" uses the sereal serialiser; only a Storable cache converts",
            ]));
        }

        // Every directory is a chromosome.
        let entries = read_dir_names(cache_dir).map_err(|error| {
            self.fail(&[
                b"cannot open ",
                cache_dir,
                b": ",
                os_text(&error).as_bytes(),
            ])
        })?;
        let mut chr_dirs: Vec<Vec<u8>> = entries
            .into_iter()
            .filter(|name| !name.starts_with(b".") && path(&join(cache_dir, name)).is_dir())
            .collect();
        chr_dirs.sort();
        if !request.chromosomes.is_empty() {
            let present: BTreeSet<&[u8]> = chr_dirs.iter().map(Vec::as_slice).collect();
            let absent: Vec<&[u8]> = request
                .chromosomes
                .iter()
                .map(Vec::as_slice)
                .filter(|name| !present.contains(name))
                .collect();
            if !absent.is_empty() {
                return Err(self.fail(&[
                    b"no chromosome directory named ",
                    &absent.join(&b" "[..]),
                    b" under ",
                    cache_dir,
                ]));
            }
            let wanted: BTreeSet<&[u8]> = request.chromosomes.iter().map(Vec::as_slice).collect();
            chr_dirs.retain(|name| wanted.contains(name.as_slice()));
        }
        self.say(&[
            b"Processing ",
            chr_dirs.len().to_string().as_bytes(),
            b" chromosomes\n",
        ]);
        // A run that stops early leaves none of these files, whatever an earlier run left.
        for name in [
            &b"info.json"[..],
            b"chr_synonyms.txt",
            b"sequence_regions.tsv",
            b"regulatory/binding_matrices.json",
        ] {
            std::fs::remove_file(path(&join(output_dir, name))).ok();
        }

        // Every chromosome directory is listed before a shard is converted, so that the
        // workers have the shards of every chromosome to take. A directory that cannot
        // be listed is reported where the run reaches it.
        let listings: Vec<std::io::Result<Vec<Vec<u8>>>> = chr_dirs
            .iter()
            .map(|chr| read_dir_names(&join(cache_dir, chr)))
            .collect();
        let cell_types = self.collect_cell_types(&chr_dirs, &listings)?;

        let mut tasks = Vec::new();
        for (chr, listing) in chr_dirs.iter().zip(&listings) {
            let Ok(entries) = listing else {
                continue;
            };
            let chr_path = join(cache_dir, chr);
            for name in shard_names(entries, b".gz") {
                tasks.push(Task::Transcripts(join(&chr_path, &name)));
            }
            for name in shard_names(entries, b"_reg.gz") {
                tasks.push(Task::Regulatory(join(&chr_path, &name)));
            }
        }
        let cell_types_ref = &cell_types;
        let work = move |task: &Task| match task {
            Task::Transcripts(input) => Outcome::Transcripts(convert_transcripts(input)),
            Task::Regulatory(input) => {
                Outcome::Regulatory(convert_regulatory(input, cell_types_ref))
            }
        };

        let mut totals = Totals::default();
        let mut matrices = Matrices::default();
        run_ordered(&tasks, request.threads, work, |pool| {
            let mut next_task = 0;
            for (chr, listing) in chr_dirs.iter().zip(&listings) {
                let mut take = || {
                    let outcome = pool.take(next_task);
                    next_task += 1;
                    outcome
                };
                self.chromosome(chr, listing, &mut take, &mut totals, &mut matrices)?;
            }
            Ok(())
        })?;

        if totals.shards == 0 {
            return Err(self.fail(&[b"no transcript shards (<start>-<end>.gz) under ", cache_dir]));
        }
        let failed = totals.shards - totals.converted;
        // A failed shard stops the run before chr_synonyms.txt and info.json are written.
        if failed > 0 {
            return Err(self.fail(&[
                failed.to_string().as_bytes(),
                b" of ",
                totals.shards.to_string().as_bytes(),
                b" shards failed to convert (see the WARN lines above)",
            ]));
        }

        let synonyms_in = join(cache_dir, b"chr_synonyms.txt");
        let synonyms_out = join(output_dir, b"chr_synonyms.txt");
        if path(&synonyms_in).exists() {
            copy_file(&synonyms_in, &synonyms_out).map_err(|error| {
                self.fail(&[
                    b"copy of chr_synonyms.txt failed: ",
                    os_text(&error).as_bytes(),
                ])
            })?;
            self.say(&[b"Wrote ", &synonyms_out, b"\n"]);
        } else {
            self.warn(&[
                b"WARN: no ",
                &synonyms_in,
                b"; the cache carries no chromosome synonyms\n",
            ]);
        }

        let characters = |bytes: &[u8]| PerlStr::bytes(bytes).characters().into_owned();
        let order = contig_order(&chr_dirs);
        info.top.insert(
            b"contig_order".to_vec(),
            Pretty::Array(
                order
                    .iter()
                    .map(|contig| Pretty::String(characters(contig)))
                    .collect(),
            ),
        );
        // A prediction marker stays only with a matrix written and a version line.
        for (tool, written) in [(&b"sift"[..], totals.sift), (b"polyphen", totals.polyphen)] {
            if !(written && info.has_source_version(tool)) {
                info.top.remove(tool);
            }
        }
        // The regulatory keys describe the written tree alone.
        for key in [&b"regulatory"[..], b"cell_types", b"regulatory_schema"] {
            info.top.remove(key);
        }
        if totals.regulatory_files > 0 {
            self.write_binding_matrices(&matrices)?;
            info.top.insert(b"regulatory".to_vec(), Pretty::True);
            info.top.insert(
                b"cell_types".to_vec(),
                Pretty::Array(
                    cell_types
                        .names
                        .iter()
                        .map(|name| Pretty::String(name.clone()))
                        .collect(),
                ),
            );
            info.top
                .insert(b"regulatory_schema".to_vec(), Pretty::Integer(1));
        }
        self.write_sequence_regions(&order, &totals.sequence_regions)?;
        self.write_info_json(&info)?;

        self.say(&[
            b"\nDone. Converted: ",
            totals.converted.to_string().as_bytes(),
            b" / ",
            totals.shards.to_string().as_bytes(),
            b" files, ",
            totals.transcripts.to_string().as_bytes(),
            b" total transcripts, ",
            totals.variation_files.to_string().as_bytes(),
            b" variation files copied, ",
            totals.regulatory_files.to_string().as_bytes(),
            b" regulatory files\n",
        ]);
        Ok(())
    }

    /// The cell-type list: every regulatory shard of the chromosomes is read once
    /// for the keys of its features' activity maps, before any shard is written.
    fn collect_cell_types(
        &mut self,
        chr_dirs: &[Vec<u8>],
        listings: &[std::io::Result<Vec<Vec<u8>>>],
    ) -> Result<CellTypes, Stop> {
        let cache_dir = &self.request.cache_dir[..];
        // The shards in conversion order, each chromosome's listing before them.
        enum Step {
            Unlisted(Vec<u8>, String),
            Shard(Vec<u8>),
        }
        let mut steps = Vec::new();
        for (chr, listing) in chr_dirs.iter().zip(listings) {
            let chr_path = join(cache_dir, chr);
            match listing {
                Err(error) => steps.push(Step::Unlisted(chr_path, os_text(error))),
                Ok(entries) => {
                    for name in shard_names(entries, b"_reg.gz") {
                        steps.push(Step::Shard(join(&chr_path, &name)));
                    }
                }
            }
        }
        let work = |step: &Step| match step {
            Step::Unlisted(..) => Ok(BTreeSet::new()),
            Step::Shard(input) => {
                let image = load_regulatory(input)?;
                Ok(regulatory::cell_type_names(&image))
            }
        };
        let mut names = BTreeSet::new();
        run_ordered(&steps, self.request.threads, work, |pool| {
            for (index, step) in steps.iter().enumerate() {
                if let Step::Unlisted(chr_path, reason) = step {
                    return Err(self.fail(&[b"cannot open ", chr_path, b": ", reason.as_bytes()]));
                }
                let Step::Shard(input) = step else {
                    continue;
                };
                match pool.take(index) {
                    Ok(found) => names.extend(found),
                    Err(error) => {
                        let error: Die = error;
                        return Err(self.fail(&[input, b": ", &Self::die_text(&error), b"\n"]));
                    }
                }
            }
            Ok(())
        })?;
        Ok(CellTypes::new(names))
    }

    /// One chromosome directory: its transcript shards, its variation files, its
    /// regulatory shards, and its line.
    fn chromosome(
        &mut self,
        chr: &[u8],
        listing: &std::io::Result<Vec<Vec<u8>>>,
        take: &mut dyn FnMut() -> Outcome,
        totals: &mut Totals,
        matrices: &mut Matrices,
    ) -> Result<(), Stop> {
        let request = self.request;
        let chr_path = join(&request.cache_dir, chr);
        let out_chr = join(&join(&request.output_dir, b"transcripts"), chr);
        let entries = match listing {
            Ok(entries) => entries,
            Err(error) => {
                return Err(self.fail(&[
                    b"cannot open ",
                    &chr_path,
                    b": ",
                    os_text(error).as_bytes(),
                ]))
            }
        };
        let files = shard_names(entries, b".gz");
        if files.is_empty()
            && entries
                .iter()
                .any(|name| shard_stem(name, b".sereal").is_some())
        {
            return Err(self.fail(&[
                &chr_path,
                b" holds sereal shards; only a Storable cache converts",
            ]));
        }
        self.make_path(&out_chr)?;
        totals.shards += files.len() as u64;
        let mut chr_transcripts = 0u64;
        for file in &files {
            let input = join(&chr_path, file);
            let output = join(&out_chr, &json_name(file, b".gz"));
            let Outcome::Transcripts(outcome) = take() else {
                unreachable!("a transcript shard's outcome");
            };
            let written = outcome.and_then(|walk| {
                // A chromosome's sequence is read from its shards in order until one
                // states a length: the slices of a later shard are not read, so a
                // match that dies there ends no walk.
                if !totals.sequence_regions.contains_key(chr) {
                    match walk.sequence {
                        Noted::Dies(die) => return Err(die),
                        Noted::Region(region) => {
                            totals.sequence_regions.insert(chr.to_vec(), region);
                        }
                        Noted::Nothing => {}
                    }
                }
                let shard = walk.shard?;
                // A shard's transcripts are counted before they are encoded.
                chr_transcripts += shard.transcripts;
                totals.sift |= shard.sift;
                totals.polyphen |= shard.polyphen;
                match shard.encode_error {
                    Some(error) => Err(error),
                    None => std::fs::write(path(&output), &shard.json)
                        .map_err(|error| Die::new(format!("write failed: {}", os_text(&error)))),
                }
            });
            match written {
                Ok(()) => totals.converted += 1,
                Err(error) => {
                    std::fs::remove_file(path(&output)).ok();
                    self.warn(&[b"WARN: ", &input, b": ", &Self::die_text(&error), b"\n\n"]);
                }
            }
        }
        totals.transcripts += chr_transcripts;

        let variation_files = self.copy_variation_files(
            &chr_path,
            &join(&join(&request.output_dir, b"variations"), chr),
            entries,
        )?;
        totals.variation_files += variation_files;

        let regulatory_files = shard_names(entries, b"_reg.gz");
        if !regulatory_files.is_empty() {
            let out_chr = join(&join(&request.output_dir, b"regulatory"), chr);
            self.make_path(&out_chr)?;
            for file in &regulatory_files {
                let input = join(&chr_path, file);
                let output = join(&out_chr, &json_name(file, b"_reg.gz"));
                let Outcome::Regulatory(shard) = take() else {
                    unreachable!("a regulatory shard's outcome");
                };
                let written = matrices
                    .record(shard.matrices, &input)
                    .and(match shard.error {
                        Some(error) => Err(error),
                        None => Ok(()),
                    })
                    .and_then(|()| {
                        std::fs::write(path(&output), &shard.json)
                            .map_err(|error| Die::new(format!("write failed: {}", os_text(&error))))
                    });
                if let Err(error) = written {
                    std::fs::remove_file(path(&output)).ok();
                    return Err(self.fail(&[&input, b": ", &Self::die_text(&error), b"\n"]));
                }
            }
        }
        totals.regulatory_files += regulatory_files.len() as u64;
        self.say(&[
            b"  chr",
            chr,
            b": ",
            files.len().to_string().as_bytes(),
            b" files, ",
            chr_transcripts.to_string().as_bytes(),
            b" transcripts, ",
            variation_files.to_string().as_bytes(),
            b" variation files, ",
            regulatory_files.len().to_string().as_bytes(),
            b" regulatory files\n",
        ]);
        Ok(())
    }

    /// The variation files: `all_vars.gz` with its `.csi` or `.tbi` index, else
    /// every `<start>-<end>_var.gz`, byte for byte.
    fn copy_variation_files(
        &mut self,
        chr_path: &[u8],
        out_chr: &[u8],
        entries: &[Vec<u8>],
    ) -> Result<u64, Stop> {
        let present = |name: &[u8]| entries.iter().any(|entry| entry == name);
        let names: Vec<Vec<u8>> = if present(b"all_vars.gz") {
            let index = [&b"all_vars.gz.csi"[..], b"all_vars.gz.tbi"]
                .into_iter()
                .find(|name| present(name));
            let Some(index) = index else {
                return Err(self.fail(&[
                    chr_path,
                    b"/all_vars.gz has no all_vars.gz.csi or all_vars.gz.tbi beside it",
                ]));
            };
            vec![b"all_vars.gz".to_vec(), index.to_vec()]
        } else {
            shard_names(entries, b"_var.gz")
        };
        if names.is_empty() {
            return Ok(0);
        }
        self.make_path(out_chr)?;
        for name in &names {
            let from = join(chr_path, name);
            copy_file(&from, &join(out_chr, name)).map_err(|error| {
                self.fail(&[b"copy of ", &from, b" failed: ", os_text(&error).as_bytes()])
            })?;
        }
        Ok(names.len() as u64)
    }

    /// `regulatory/binding_matrices.json`: one object, the matrices by stable id in
    /// sorted order.
    fn write_binding_matrices(&mut self, matrices: &Matrices) -> Result<(), Stop> {
        let out_path = join(
            &self.request.output_dir,
            b"regulatory/binding_matrices.json",
        );
        let mut table = vec![b'{'];
        for (index, (id, (json, _))) in matrices.by_id.iter().enumerate() {
            if index > 0 {
                table.push(b',');
            }
            let id = PerlStr {
                bytes: id,
                utf8: true,
            };
            // An id reached a shard as a JSON string before it reaches this table.
            write_string(id, &mut table).ok();
            table.push(b':');
            table.extend_from_slice(json);
        }
        table.extend_from_slice(b"}\n");
        self.write_metadata(&out_path, &table)?;
        self.say(&[
            b"Wrote ",
            &out_path,
            b" (",
            matrices.by_id.len().to_string().as_bytes(),
            b" binding matrices)\n",
        ]);
        Ok(())
    }

    /// `sequence_regions.tsv`: the header line, then the row of each chromosome of the
    /// contig order whose transcripts state a length.
    fn write_sequence_regions(
        &mut self,
        order: &[Vec<u8>],
        stated: &BTreeMap<Vec<u8>, SequenceRegion>,
    ) -> Result<(), Stop> {
        let out_path = join(&self.request.output_dir, b"sequence_regions.tsv");
        let mut table = regions::HEADER.to_vec();
        let (mut rows, mut wide_rows) = (0u64, 0u64);
        for name in order {
            if let Some(region) = stated.get(name) {
                wide_rows += u64::from(region.write_row(name, &mut table));
                rows += 1;
            }
        }
        let file = self.create_metadata(&out_path)?;
        // Perl warns as it prints each row, after the file is opened.
        for _ in 0..wide_rows {
            self.warn(&[b"Wide character in print\n"]);
        }
        self.fill_metadata(file, &out_path, &table)?;
        self.say(&[
            b"Wrote ",
            &out_path,
            b" (",
            rows.to_string().as_bytes(),
            b" sequences)\n",
        ]);
        Ok(())
    }

    fn write_info_json(&mut self, info: &Info) -> Result<(), Stop> {
        let out_path = join(&self.request.output_dir, b"info.json");
        let text = match write_pretty(&info.document()) {
            Ok(text) => {
                let bytes = print_bytes(&text);
                if bytes == text && !text.is_ascii() {
                    // A character past 255: Perl warns and writes the text as UTF-8.
                    self.warn(&[b"Wide character in print\n"]);
                }
                bytes
            }
            Err(error) => {
                // JSON::XS dies here outside any `eval`.
                self.warn(&[&error.0, b"\n"]);
                return Err(Stop(255));
            }
        };
        self.write_metadata(&out_path, &text)?;
        self.say(&[b"Wrote ", &out_path, b"\n"]);
        Ok(())
    }

    /// A file of the output's root, written whole or not left behind.
    fn write_metadata(&mut self, out_path: &[u8], bytes: &[u8]) -> Result<(), Stop> {
        let file = self.create_metadata(out_path)?;
        self.fill_metadata(file, out_path, bytes)
    }

    fn create_metadata(&mut self, out_path: &[u8]) -> Result<std::fs::File, Stop> {
        std::fs::File::create(path(out_path)).map_err(|error| {
            self.fail(&[
                b"cannot write ",
                out_path,
                b": ",
                os_text(&error).as_bytes(),
            ])
        })
    }

    fn fill_metadata(
        &mut self,
        mut file: std::fs::File,
        out_path: &[u8],
        bytes: &[u8],
    ) -> Result<(), Stop> {
        file.write_all(bytes).map_err(|error| {
            std::fs::remove_file(path(out_path)).ok();
            self.fail(&[
                b"write of ",
                out_path,
                b" failed: ",
                os_text(&error).as_bytes(),
            ])
        })
    }
}

#[derive(Default)]
struct Totals {
    shards: u64,
    converted: u64,
    transcripts: u64,
    variation_files: u64,
    regulatory_files: u64,
    sift: bool,
    polyphen: bool,
    /// The sequence each chromosome directory's transcripts state, by directory name.
    sequence_regions: BTreeMap<Vec<u8>, SequenceRegion>,
}

/// The binding matrices the shards name: each once under its stable id, as its
/// canonical JSON, with the shard that first carried it.
#[derive(Default)]
struct Matrices {
    by_id: BTreeMap<Vec<u8>, (Vec<u8>, Vec<u8>)>,
}

impl Matrices {
    /// Records every matrix of a shard, in feature order: a repeat that differs from
    /// the first copy is corrupt.
    fn record(&mut self, records: Vec<regulatory::MatrixRecord>, file: &[u8]) -> Result<(), Die> {
        for record in records {
            match self.by_id.get(&record.id) {
                Some((json, source)) => {
                    if *json != record.json {
                        return Err(Die([
                            b"binding matrix ",
                            &record.id[..],
                            b" differs from its copy in ",
                            &PerlStr::bytes(source).characters(),
                        ]
                        .concat()));
                    }
                }
                None => {
                    self.by_id.insert(record.id, (record.json, file.to_vec()));
                }
            }
        }
        Ok(())
    }
}

fn convert_transcripts(input: &[u8]) -> Result<Walk, Die> {
    match load_shard(path(input))? {
        Some(image) => Ok(transcripts::walk(&image)),
        None => Err(Die::new(
            "fd_retrieve returned no data (truncated Storable image)",
        )),
    }
}

fn load_regulatory(input: &[u8]) -> Result<crate::storable::Image, Die> {
    match load_shard(path(input))? {
        Some(image) => {
            regulatory::check_shape(&image)?;
            Ok(image)
        }
        None => Err(Die::new("invalid regulatory Storable image")),
    }
}

fn convert_regulatory(input: &[u8], cell_types: &CellTypes) -> RegulatoryShard {
    match load_regulatory(input) {
        Ok(image) => regulatory::convert(&image, cell_types),
        Err(error) => RegulatoryShard {
            json: Vec::new(),
            matrices: Vec::new(),
            error: Some(error),
        },
    }
}

/// Runs the conversion and returns the exit code. Standard output is written through
/// `out` and is the caller's to flush.
pub fn run(request: &Request, out: &mut impl Write, err: &mut impl Write) -> i32 {
    let mut run = Run { request, out, err };
    match run.run() {
        Ok(()) => 0,
        Err(Stop(code)) => code,
    }
}
