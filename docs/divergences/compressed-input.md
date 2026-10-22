# Compressed bytes on standard input, and a compressed stream cut short

A difference in which inputs are read, on two shapes no measured dataset has: gzip or BGZF bytes arriving on standard input or through a pipe named by path, and a compressed stream cut inside a member; nothing is set aside.

Ensembl VEP decompresses an input only when it is a regular file whose bytes Perl's `-B` test calls binary (`Parser.pm` line 242). Such a file is opened through `get_compressed_filehandle`, which reads it with `gzip -dc` (`Utils.pm` lines 429-453, the pipe at line 441). vep-rs reads the first bytes of whatever it is given and decompresses a gzip or BGZF stream wherever it comes from.

## Compressed bytes on standard input

`-i STDIN` hands Ensembl VEP's parser the raw `*STDIN` handle (`Parser.pm` lines 235-237), and a pipe or FIFO named by path fails the regular-file test at line 242, so neither is decompressed. The parser reads the compressed bytes as text lines (`TextParser.pm` line 125) and finds no record in them. The run writes a header-only output and exits 0. Perl's `isn't numeric` and `Use of uninitialized value` warnings appear on standard error, and the warnings file holds skipped-line warnings for the bytes it read as records.

vep-rs decompresses gzip or BGZF bytes on standard input, and on a pipe named by path, as it decompresses a file, so the piped stream annotates exactly as the file does. `--format` is still required on standard input, as in Ensembl VEP. The 39 records of the `tests/golden/116/GRCh38-input-compressed` corpus piped as BGZF bytes give Ensembl VEP 116.2 no row, where the file gives it 977; vep-rs writes the 981 rows it writes for the file.

## A stream cut inside a member

Ensembl VEP reads the `gzip -dc` pipe to its end and never reads the program's exit status: nothing in its parser or runner closes the input handle or tests `$?` after a close (the one such check is the output compressor's, `Runner.pm` lines 212-224). A BGZF or plain gzip file cut inside a member therefore yields every record gzip inflated before the cut, and the run exits 0. gzip's own `unexpected end of file` appears on standard error twice, once for each parser object that opens the file (`Parser.pm` lines 133 and 159). The output looks complete and is not.

vep-rs ends the run at the cut: one line `ERROR: [run] Failed to read input ...` on standard error, exit 1, and nothing past the header written (a cut inside the block that holds the header leaves no output file). It does so on the BGZF path under every decompression worker count and on the plain gzip path alike. The same corpus written one BGZF block per record and cut inside the block of record 19 gives Ensembl VEP 116.2 the 849 rows of the 18 complete records. As a single-member gzip stream cut 60 bytes short it gives the 879 rows of the 36 complete records. vep-rs writes no row for either file.

A stream cut exactly at a block boundary, a BGZF file without its end-of-file block, and bytes after the last member of a plain gzip stream are read alike by both engines.

To reproduce: Ensembl VEP 116.2 with the GRCh38 cache and [the command of the closing section](reproducing.md#reproducing-these-examples), adding `-i STDIN --format vcf --max_sv_size -1` with a `bgzip` copy of the input on standard input. For the second difference, name a compressed copy of the input cut short by `head -c`; the row counts above hold under `--max_sv_size -1`, the corpus's own flag. vep-rs takes the same flags and inputs over the same cache.
