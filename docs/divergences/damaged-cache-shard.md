# A transcript shard the cache cannot read

Two differences in which records are written, on cache states no measured run has: a transcript file of the cache cut short after its header, and a cache holding no transcript file at all; nothing is set aside.

Ensembl VEP's transcript cache holds one Storable file per 1 Mb region, gzip-compressed. A region's file is opened as a region is first asked for, through `PerlIO::gzip`, the `gzip` binary or `Compress::Zlib` (`AnnotationSource/Cache/BaseSerialized.pm` lines 134, 140 and 146), and a file that cannot be opened ends the run with the open's error. The bytes are then handed to Storable's `fd_retrieve` (lines 135, 141 and 155), which dies on a stream whose header is not Storable's, so a file of foreign bytes ends the run too.

A file cut after its header is read differently. `fd_retrieve` wraps the C retrieval in an `eval` and dies only on an exception (`Storable.pm` lines 435 to 449 of the Storable 3.23 in the release 116.2 image). A stream that ends inside the data returns `undef` and raises none, through each of the three decompressors. The caller stores whatever came back as the region's transcripts (`BaseSerialized.pm` lines 261 to 265), and the transcript loader reads `undef` as an empty table (`AnnotationSource/Cache/Transcript.pm` line 236). The region is recorded as holding no transcript, and the run exits 0. When the `gzip` binary is the decompressor its `unexpected end of file` appears on standard error; the other two decompressors report nothing.

A variant whose overlapping transcripts the cut file alone held is written as `intergenic_variant`. A variant within 5,000 bases of the region's edge can still meet a transcript that spans the boundary. The regions read for a variant extend by that distance on each side (`AnnotationSource.pm` line 190, `AnnotationType/Transcript.pm` lines 190 to 196). The cache stores a transcript in every region file it overlaps, which is why the loader drops a transcript already seen (`AnnotationType/Transcript.pm` lines 282 to 287).

vep-rs loads every shard of a chromosome when the first record on it arrives. A shard whose bytes cannot be read, inflated or parsed ends the run before any record of that chromosome is annotated, with the chromosome and the loader's own words for the shard, exit 1. The rows of batches on other chromosomes already written stay in the output:

```
ERROR: [cache] Cannot read the transcript cache for chr21: Failed to inflate
<cache>/transcripts/21/42000001-43000000.json.gz: incomplete deflate stream
```

The line is wrapped here and is one line on standard error. A `c.`, `n.` or `p.` notation whose lookup reaches the shard first gives the same refusal after `Cannot read the transcript cache for HGVS input: `. A cache that carries its [checksum manifest](../cache-checksums.md) refuses the cut shard earlier, as a checksum mismatch, before any byte of it is inflated.

To reproduce in Ensembl VEP 116.2: copy the GRCh38 cache and cut one region file of chromosome 21 short with `head -c` to half its bytes. Run [the command of the closing section](reproducing.md#reproducing-these-examples) over a record more than 5,000 bases inside that region; the record is written as `intergenic_variant`. In vep-rs: copy a JSON cache without its `CHECKSUMS.sha256`, cut one shard of chromosome 21 the same way, and run any record on chromosome 21 over the copy.

## A cache holding no transcript shard

Ensembl VEP's valid chromosomes are the cache's directories whose names do not start with a dot, whatever they hold (`CacheDir.pm` line 431), and a region whose file does not exist yields no feature (`AnnotationSource/Cache/BaseSerialized.pm` line 263, `if -e $file`). So a cache none of whose chromosome directories holds a transcript file runs to the end: every record on a directory's name is written as `intergenic_variant`, exit 0, and nothing reports the cache empty.

vep-rs lists the shards of every chromosome directory as the cache opens, before the input is read, and refuses a cache that holds none, a dot-prefixed directory's files not counted:

```
ERROR: [cache] no transcript shards under <cache>/transcripts (serialiser_type
storable)
```

The line is one line on standard error, exit 1, and names the serialiser `info.json` declares (`storable` when it declares none). One directory holding no shard beside others that do is a chromosome with no transcripts in both engines, and a record on it is `intergenic_variant` in both.

To reproduce in Ensembl VEP 116.2: a cache directory holding `info.txt` and an empty `21/`; [the command of the closing section](reproducing.md#reproducing-these-examples) writes any record on chromosome 21 as `intergenic_variant`. In vep-rs: a cache directory holding `info.json` and an empty `transcripts/21/`.
