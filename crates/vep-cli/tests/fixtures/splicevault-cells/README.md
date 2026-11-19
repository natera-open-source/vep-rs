# SpliceVault cells fixture

A made-up file of the twelve columns the Ensembl VEP plugin splits a row into
(`SpliceVault.pm:201`), its header
line a comment, bgzip-compressed and indexed with `tabix -s 1 -b 2 -e 2`. Its
eight rows are fabricated, at positions of the GRCh38 SpliceVault fixture corpus
and naming its transcripts: a row at `chr21:13618457` for `T>C` whose tenth
cell (`top_events`) is empty; a row at `chr21:13629903` for `A>AT` whose second
event carries a `,` inside its percent part (`Top2:CA;+112;1,5%;inFrame`); a row
at `chr21:13630940` for `G>T` whose twelfth cell (`max_depth`) is empty and is
followed by a thirteenth cell; a row at `chr21:13630943` for `T>C` of eleven
cells, with no `max_depth` cell; a row at `chr21:13631108` for `G>A` with
`max_depth` `93`; a row at `chr21:18294269` for `C>A` of twelve cells, the
twelfth empty (the line ends in a tab); a row at `chr21:18297830` for `C>G`
whose eleventh cell (`sample_count`) opens with a blank (` 33837`) and whose
twelfth ends with one (`26 `); and a row at `chr21:18315257` for `C>A` whose
first event ends in its separator (`Top1:CA;-480;7%;`), an empty fifth part.
`variants.vcf` holds the eight records.

`SpliceVault.pm` line 201 splits a row into twelve variables, so Perl keeps a
trailing empty cell as the empty string and leaves a missing cell undefined,
and lines 206 to 220 write both into the plugin's result. The events cell is
split into a list on `|` (line 216), empty when the cell is, and the text
formats join that list with `,` (`Utils.pm` lines 185 to 197) or `&` (`VCF.pm`
line 391).

Line 184 splits an event into five variables on `:` alone, so a trailing empty
part is the empty string and a `,` inside a part is text. In the JSON format
`numberify` (`Utils.pm` line 231) converts every value `looks_like_number`
admits, blanks at either end included, and writes the number without them.

Ensembl VEP 116.2 wrote, on this fixture: for the two rows whose twelfth cell
is empty, `SpliceVault_site_max_depth=` in the default format, an empty cell in
the tab format, an empty CSQ field in the VCF format and
`"SpliceVault_site_max_depth": ""` inside the `splicevault` object in JSON; for
the eleven-cell row, `SpliceVault_site_max_depth=` in the default format, `-` in
the tab format, an empty CSQ field and `null` inside that object (`VEP_output.pm`
lines 117 to 119; `Tab.pm` line 128; the JSON format encodes the plugin's
nested hash as it stands, `JSON.pm` line 327 dropping top-level keys alone);
for the scored row, `93` in every format; for the row whose events cell is
empty, `SpliceVault_top_events=` in the default format, an empty cell in the
tab format (the empty list is defined), an empty CSQ field and `{}` in JSON;
for the padded row, `SpliceVault_site_max_depth=26 ` and
`SpliceVault_site_sample_count= 33837` in the default format, `26 ` and ` 33837`
in the tab format, `26_` and `_33837` in the VCF format (`VCF.pm` line 403
rewrites whitespace) and the numbers `26` and `33837` in JSON; for the row whose
first event ends in its separator, `Top1:CA:-480:7%:,Top2:CA:-76:39%:Frameshift`
in the default and tab formats, the same with `&` in the VCF format, and in JSON
`"Top1": {"frame": "", ...}`, the empty part written as the empty string; and
for the row whose second event carries a `,`,
`Top1:CD:-330:73_314%:inFrame,Top2:CA:+112:1,5%:inFrame` in the default and tab
formats, `Top1:CD:-330:73_314%:inFrame&Top2:CA:+112:1&5%:inFrame` in the VCF
format (`VCF.pm` line 401 rewrites every `,`), and in JSON two events, `Top2`
with `"percent_of_supporting_samples": "1,5%"`.
