# Cache checksum manifest

`vep-cache-builder` and `vep-cache-converter` end a build by writing `CHECKSUMS.sha256` in the cache directory, beside `info.json`. The file lists every file of the directory except itself, one `sha256sum` line each: the digest, two spaces, the path relative to the cache directory. Run there, `sha256sum --check CHECKSUMS.sha256` verifies every file. [Cache setup](cache-setup.md) describes the directory the manifest covers.

## What a run checks

A run reads the manifest once as it resolves the cache and first checks that every listed file is present. It then hashes each cache file it reads whole as it reads it and compares the digest with the entry. Those files are the transcript shards, `info.json`, `chr_synonyms.txt`, each variation index (`.csi` or `.tbi`), each region file (`<start>-<end>_var.gz`), and the `.fai` and `.gzi` of a FASTA found in the cache directory. A file whose digest differs from its entry, and a file the manifest does not list, end the run with exit 1 and one line:

```
ERROR: [cache] checksum mismatch for transcripts/21/1-1000000.json.gz in cache
```

A listed file that is missing ends the run before anything is read, as `ERROR: [cache] <path> is listed in CHECKSUMS.sha256 but missing from <dir>`. A missing `info.json` or `transcripts/` is refused earlier, with the cache resolution's own lines.

## What it does not protect

`all_vars.gz` and a FASTA are read by byte range through their index, so a run does not hash them. Their index is hashed, the file must be listed, and a copy that ends before the last block or base its index names is refused. A change inside such a file that keeps its index and length is not detected. Under `--format id` the variation files are read start to end and hashed whole. The manifest is not a signature: whoever can change a file can rewrite it.

## Caches without a manifest

A cache without the file loads and annotates unchanged. The run warns once, in the warnings file and on standard error unless `--quiet`:

```
WARNING: cache /path/to/cache carries no checksum manifest
```

Caches built or converted by releases before 0.5.0 carry none, and deleting the manifest returns a cache to that state. Build or convert the cache again, or write the manifest in the cache directory:

```
find -L . -type f ! -name CHECKSUMS.sha256 | sort | xargs sha256sum \
  > CHECKSUMS.sha256
```

A FASTA copied into a cache directory that carries a manifest is refused until the manifest lists its files; the command above lists them. A `--fasta` or `--synonyms` file is the user's own, wherever it sits, and is never checked.
