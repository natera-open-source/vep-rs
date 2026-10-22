# Control bytes in output values

A difference in the text of a field, not in any tuple: the comparators key rows by Location, Allele and Feature, which are ASCII, and nothing is set aside.

A value read from an input file can carry an ASCII control byte: an identifier, a chromosome name or a sample column with a carriage return, a form feed, an escape or a NUL in it, a plugin data field with one. Ensembl VEP writes such a value as it read it.

Its default output joins the columns with tabs and rewrites a `;` in an Extra value as `%3B` (`OutputFactory/VEP_output.pm` lines 115-135); the tab output joins the columns and rewrites nothing (`OutputFactory/Tab.pm` line 128). A VCF `CSQ` value has `,` and `|` rewritten as `&`, `;` as `%3B` and a run of whitespace as one `_` (`OutputFactory/VCF.pm` lines 401-404). Every other byte reaches the file as read, so a line end in a default or tab value ends the row and a NUL sits in the text.

vep-rs writes every byte from `0x00` to `0x1F`, and `0x7F`, in a default or tab value as its percent-encoding: `%09` for a tab, `%0A` for a line end, `%0D` for a carriage return, `%00` for a NUL. The rule is the one the Extra column already applies to `;`. Inside a VCF `CSQ` value the whitespace rule comes first, so a tab, a line end, a carriage return, a form feed or a vertical tab is `_` in both engines, and every other control byte is percent-encoded. The JSON output of both engines escapes control characters as JSON requires.

The columns a VCF line copies from the input, the ID among them, stay the input's bytes in both engines, as every format's header lines do. The VCF line vep-rs builds for a record of another input format has each column encoded like any value; Ensembl VEP writes that line's columns as read.

A value with no control byte is written the same by both engines. A row whose field carries one differs in that field's text alone and is charged to vep-rs in both figures.
