# A BGZF block whose trailer overstates its size

A difference in which records are written, on an input shape no measured dataset has; nothing is set aside.

A BGZF file is a sequence of gzip members of at most 64 KiB of inflated data each, and every member's trailer declares its inflated size (`ISIZE`). Ensembl VEP reads a compressed input file through `gzip -dc`, or through `IO::Uncompress::Gunzip` where no `gzip` is on the path (`Utils.pm` lines 440-446, the multi-member branch the parser asks for), decompressors that inflate a member's data before they read its trailer. A block whose `ISIZE` declares more than the data inflates to is read through: `gzip -dc` writes the inflated records, reports `invalid compressed data--length error` and ends the stream there; `IO::Uncompress::Gunzip` writes them and reports nothing. Either way the records reach the parser and are annotated.

vep-rs reads the file block by block and checks each block's header and trailer before the block is inflated. A block whose `ISIZE` exceeds 65536, or whose declared size leaves no room for its trailer, ends the run before any record of it is written:

```
ERROR: [input] malformed BGZF block at offset <n> in <file>
```

`<n>` is the block's byte offset in the file. A refused first block writes no record. A block refused later in the file ends the run at that block: the rows the pipeline had already written stay in the output, and the records read but still unannotated are dropped with the run. The exit code is 1 under one decompression thread or several.

A block the file ends inside its data is reported in the block reader's own words under `ERROR: [run]`. A file that ends inside a block's header is refused as a malformed block at that block's offset, which the block reader alone would read as the end of the file. Both of Ensembl VEP's decompressors write the records of the whole blocks before the cut and then report it (`unexpected end of file`; `Header Error: Truncated in FEXTRA Body Section`), so Ensembl VEP annotates those records where vep-rs ends the run.

A member that is not BGZF at all, a plain gzip member after the BGZF blocks, is inflated whole by both engines, as `gzip -dc` inflates every member of a stream, and its records are annotated.
