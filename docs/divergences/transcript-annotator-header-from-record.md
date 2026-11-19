# A TranscriptAnnotator table whose record at position 1 of contig `1` becomes its header

**The defect in one sentence.** Ensembl VEP's TranscriptAnnotator plugin reads the table's column names from the last line a `tabix -h` query for the one-base region at position 1 of the contig named `1` prints, and that command prints the header lines and then every record in the region. A table with a record at that position therefore takes the record's cells as its column names. vep-rs reads the column names from the table's last header line, whatever records the table holds.

**Why this is a defect.** The region is used as a query expected to match no record, so that the output is the header alone. A record at that position is a legal row of the table, and the plugin never checks whether the last line printed is a header line. The columns then carry the record's values as names (`<prefix>_0.1`, `<prefix>_Bogus`), and each field of every annotated row is written under a name no header stated. Five warnings that the column names should start with `chr`, `pos`, `ref`, `alt` and `tr` are the only sign.

**What vep-rs does instead.** vep-rs reads the column names from the table's last header line, the last line beginning with `#`, whatever records the table holds and wherever they lie.

**Where it is in Ensembl VEP's source** (`VEP_plugins` at release/116):

- `TranscriptAnnotator.pm`, lines 99-105, `new`: `tabix -f -h` is opened on the file with the region of contig `1` from position 1 to position 1, `$last = $_ while <IN>` keeps its last line, and the leading `#` is stripped from whatever that line is. Line 109 splits it into the column names, and lines 113-118 warn when the first five do not start with `chr`, `pos`, `ref`, `alt` and `tr`.

**What correcting it would change.** The column names would come from the header whatever the table's first contig holds at position 1, as vep-rs reads them.

**The class it produces.** `transcript_annotator_header_from_record_at_1_1`: every `<prefix>_<column>` field of the table (the JSON keys alike) on every format, since the field names themselves differ, on a run whose table has a record at position 1 of contig `1`. The class is documented only: the committed tables hold no such record, and no ordinary table does, since position 1 of a chromosome is telomeric sequence, so no manifest names the class and no comparator applies it.

## 1. The corpus table with one record at position 1 of contig `1`

**In short.** The `GRCh38-plugin-transcript-annotator` corpus's table, with the record `1 1 N A ENST00000000001 0.1 Bogus RecordAtOneOne` added before its chromosome 21 rows, is run on the corpus records `rs463608`, `rs2074312880` and `rs2074361249`. Ensembl VEP's header lines are `## oneone_0.1 : column from <file>`, `## oneone_Bogus : ...` and `## oneone_RecordAtOneOne : ...`, and its `ENST00000618007` row for `rs463608` reads `oneone_0.1=0.61;oneone_Bogus=Tolerated;oneone_RecordAtOneOne=Intronic`. vep-rs writes `oneone_Comment=Intronic;oneone_SIFT_pred=Tolerated;oneone_SIFT_score=0.61` under the header's names. The values agree on every row; only the names differ.

**Try it yourself.** Prepend the record to `tests/golden/116/GRCh38-plugin-transcript-annotator/data/transcript_annotator.tsv.gz`, `bgzip` and `tabix -s 1 -b 2 -e 2` the copy, and run both engines with `--plugin TranscriptAnnotator,file=<copy>` on the corpus's `variants.vcf`; the commands of [Reproducing these examples](reproducing.md#reproducing-these-examples) give the engines and their caches.
