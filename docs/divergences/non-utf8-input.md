# Input bytes that are not UTF-8

A difference in the text of a field, not in any tuple: the comparators key rows by Location, Allele and Feature, which are ASCII, and nothing is set aside.

Ensembl VEP opens its input with no decoding layer (`Parser.pm` lines 235-242 for standard input and a compressed file, `ensembl-io`'s `TextParser.pm` line 70 for a plain file) and writes the bytes it read. A Latin-1 `é` (the byte `0xE9`) in a header line, an ID, an INFO value or a sample name comes out as that byte in every format. In the JSON output that byte leaves the document invalid as UTF-8.

vep-rs reads the input as bytes and writes the VCF output as bytes. The header block carries the input's `##` lines and sample names unchanged, and every record line, annotated or a non-variant record kept by `--allow_non_variant`, is the line it was read from with the annotation appended to INFO. So a VCF written from a VCF carries every input byte.

A record is otherwise read as UTF-8 with each invalid sequence replaced by U+FFFD, sample columns included. An ID, an INFO value or a genotype field that is not UTF-8 prints the replacement character in `Uploaded_variation`, in the JSON `input` string (`serde_json` writes UTF-8 alone) and in a warning that quotes the line. Of those columns the default and tab outputs print the ID alone. A record of another input format is read as text throughout, so its name reaches the ID column of the VCF line built for it with the replacement character too.

A file that is UTF-8 throughout, non-ASCII characters included, is written by both engines as read.
