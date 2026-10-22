# Empty lines in a VCF input

A difference in which records are written, on an input shape no measured dataset has; nothing is set aside.

Ensembl VEP's parsers for the HGVS, ID, SPDI, region and CAID formats skip a run of empty lines and warn `Skipped N empty line(s) from line M` (`Parser.pm` lines 266-293, `skip_empty_lines`).

Its parser for the Ensembl default format calls the same routine but holds each line as a list of columns. The test for an empty line (`Parser.pm` line 273) reads that list as a word, so an empty line is never skipped. The record is read with no allele string, `get_SO_term` is called for it (`Parser/VEP_input.pm` line 137) and dies on a method the parser lacks (`Parser.pm` line 665), and the run ends there with an exception.

Its VCF parser does not call that routine: an empty line is a record with no column, for which it yields no variant (`Parser/VCF.pm` line 226). The input buffer reads that as the end of its fill (`InputBuffer.pm` line 212). A fill that begins on an empty line holds nothing, and the run ends with the output written so far (`Runner.pm` line 424). So a VCF whose first data line is empty, or that has two empty lines in a row, is annotated up to that point and no further, with no message.

A line of blanks alone is not empty to that parser: it is a record of one column, refused with `Invalid start '' or end '-1' coordinate` (`Invalid start '  '` when a tab leads the blanks). It counts as a record too, so the `line N` of every later warning is one higher for each such line.

vep-rs skips a run of empty or whitespace-only lines in a VCF, annotates every record around it, and writes the warning in `skip_empty_lines`'s words. Its `M` is the number of the record that follows the run. In the Ensembl default format vep-rs skips the run with the same warning and annotates the records around it, where Ensembl VEP's run ends; in the HGVS, ID, SPDI and region formats both engines skip and warn alike.
