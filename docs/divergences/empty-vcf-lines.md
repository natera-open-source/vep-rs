# Empty lines in a VCF input

A difference in which records are written, on an input shape no measured dataset has; nothing is set aside.

Ensembl VEP's parsers for the Ensembl default, HGVS, ID, SPDI and region formats skip a run of empty lines and warn `Skipped N empty line(s) from line M` (`Parser.pm` lines 266-293, `skip_empty_lines`).

Its VCF parser does not call that routine: an empty line is a record with no column, for which it yields no variant (`Parser/VCF.pm` line 226). The input buffer reads that as the end of its fill (`InputBuffer.pm` line 212). A fill that begins on an empty line holds nothing, and the run ends with the output written so far (`Runner.pm` line 424). So a VCF whose first data line is empty, or that has two empty lines in a row, is annotated up to that point and no further, with no message.

vep-rs skips a run of empty or whitespace-only lines in a VCF, annotates every record around it, and writes the warning in `skip_empty_lines`'s words. Its `M` is the number of the record that follows the run.
