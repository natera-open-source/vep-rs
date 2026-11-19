# The process id a forked run writes before a warning

**The defect in one sentence.** Under `--fork`, Ensembl VEP's parent process relays a child's warnings into the warnings file as one block headed by the child's process id, so the block's first line reads `WARNING: <pid> : WARNING: <text>` where the same run without the flag writes `WARNING: <text>`. vep-rs writes `WARNING: <text>`, with no process id, whatever the number of workers.

**Why this is a defect.** A warning names a condition of the input or the data, and the same input gives the same warning whatever the number of processes a run was given. The process id is a number the operating system hands out, different from run to run and from host to host, and it carries no fact about the record or the cache. Written into the warnings file it makes two runs of one input on one engine differ by nothing but the day they ran, and it puts a second `WARNING:` into the line a script parses. The lines after a block's first are written bare, so one block mixes the two shapes.

**What vep-rs does instead.** vep-rs writes every warning as the condition alone, `WARNING: <text>`, at every worker count, the shape Ensembl VEP writes without `--fork`. Its annotation rows, their order and its statistics are the same at every count, as Ensembl VEP's are.

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.2):

- `Bio/EnsEMBL/VEP/Runner.pm`, lines 387-392, `next_output_line`: a true `--fork` value takes the forked path; no flag or `--fork 0` is the serial one.
- `Bio/EnsEMBL/VEP/Runner.pm`, lines 659-662, `_forked_process`: the child points the warning handle at its standard error, then closes and reopens standard error into a string, so every warning the child raises while annotating its chunk accumulates there.
- `Bio/EnsEMBL/VEP/Runner.pm`, lines 711-719: the child sends the string back to the parent beside its output.
- `Bio/EnsEMBL/VEP/Runner.pm`, line 565, `_forked_buffer_to_output`: the parent passes the child's process id, ` : ` and the whole captured string to `warning_msg` as one message.
- `Bio/EnsEMBL/VEP/BaseVEP.pm`, lines 802-822, `warning_msg`: a message that does not begin with `warn` gets `WARNING: ` before it and is written to the warnings file once, so the relayed block opens `WARNING: <pid> : ` and the child's own lines follow inside it.

**What correcting it would change.** The parent would write the child's lines as the child wrote them, so a run's warnings file would read the same at every fork count, and a warning would carry one `WARNING:`.

**The class it produces.** `fork_warning_pid_prefix`: a warnings-file line that differs between the engines by the `<pid> : ` Ensembl VEP's parent wrote before a forked child's text, the text after it equal to vep-rs's line. The class covers the first line of each relayed block; the block's later lines are bare in both engines. It covers a child's interpreter warning as well, which the relay writes as `WARNING: <pid> : <text>` with no second `WARNING:`. The format comparator's warnings mode sets the prefix aside under this class and counts the lines, and the golden datasets record each reference warning as the condition alone, so their comparison reads the serial shape.

**How often the shape occurs.** Tier 2: every committed reference run was made under `--fork 1`, so every flag set takes the forked path, and 12 of the 392 committed flag sets raise a warning inside a child's annotation. Six, in 2 datasets, have a RefSeq transcript whose sequence disagrees with the assembly; one has a plugin that returns the wrong kind of value; five raise an interpreter warning while annotating, one under `--coding_only` and four, in 2 datasets, in the Downstream plugin. The reference warnings of those 12 flag sets, in 6 datasets, carry the prefix.

## 1. A transcript-assembly mismatch under `--fork 1`

**In short.** The record `rs922296917` of the GRCh38 RefSeq dataset (`tests/golden/116/GRCh38-refseq/variants.vcf`), a deletion whose RefSeq transcript carries a base the assembly does not, raises `Transcript-assembly mismatch in rs922296917`. Ensembl VEP 116.2 writes `WARNING: Transcript-assembly mismatch in rs922296917` to the warnings file under `--fork 0` and `WARNING: 11 : WARNING: Transcript-assembly mismatch in rs922296917` under `--fork 1` and `--fork 4`. vep-rs writes the first form at every worker count. The six annotation rows are the same in every run of either engine.

### 1.1 What a reader needs to know first

A RefSeq cache carries transcripts aligned to the genome whose sequence may disagree with the assembly at some bases. When a variant's transcript does, Ensembl VEP writes one `Transcript-assembly mismatch in <record>` warning per record and annotates from the transcript's own sequence, as [RefSeq and merged caches](../output/refseq-and-merged.md) describes. `--fork N` runs the annotation of a buffer in N child processes; `--fork 0` and no flag run it in the parent. vep-rs's `--fork` is a thread count, parallel by default, and the warnings file is the same file either engine writes beside the output. Every committed reference output was generated under `--fork 1`.

### 1.2 The record

```
21	45505371	rs922296917	CGGCCCTCCGGGCCCCCCTG	C	.	.	.
```

Both engines run `--refseq` against the dataset's RefSeq cache and FASTA, Ensembl VEP once at each of `--fork 0`, `--fork 1` and `--fork 4`, vep-rs at `--fork 1` and `--fork 4`.

### 1.3 What each engine writes

Ensembl VEP 116.2's warnings file under `--fork 0`:

```
WARNING: Transcript-assembly mismatch in rs922296917
```

and under `--fork 1` and `--fork 4`:

```
WARNING: 11 : WARNING: Transcript-assembly mismatch in rs922296917
```

vep-rs's warnings file under `--fork 1` and `--fork 4`:

```
WARNING: Transcript-assembly mismatch in rs922296917
```

The six rows of the default output are identical across the five runs. The committed reference of the dataset, `expected/default.warnings.log`, holds the forked form on its first line and the dataset's second mismatch bare on its second: the two warnings were raised in one child and relayed as one block. Without `--quiet`, the child's own copy of the line on standard error is captured as well, so the block carries the warning twice.

### 1.4 Why Ensembl VEP writes that

Under `--fork` the child annotates its chunk with standard error redirected into a string, so the warning the transcript raises lands there instead of in the file. The child sends the string to the parent with its rows. The parent passes `11 : ` followed by the string to `warning_msg`, which prefixes `WARNING: ` because the message begins with a digit, and writes the whole block as one warning. The `11` is the process id the child was given inside the container; another host gives another number. Under `--fork 0` the warning is raised in the parent, whose handle is the file, and is written as the child would have written it.

### 1.5 Why vep-rs's row is the correct one

The warning describes the transcript and the record, and both are the same whether one process or four annotated the buffer. The process id describes neither; it is a run-time accident the file preserves. vep-rs's line is the line Ensembl VEP itself writes for this record when nothing is relayed, and it is the same line at every worker count, so two runs of one input agree on their warnings as they agree on their rows.

### 1.6 Try it yourself

```
vep -i record.vcf -o out.txt --offline --cache --dir_cache <corpus>/cache_root \
  --cache_version 116 --assembly GRCh38 --refseq \
  --fasta <corpus>/reference.fa.gz --fork 4
cat out.txt_warnings.txt
```

vep-rs prints `WARNING: Transcript-assembly mismatch in rs922296917`. Ensembl VEP 116.2 prints the same line after `WARNING: <pid> : `, and the bare line when `--fork 4` is dropped.

### 1.7 How the comparator treats this class

The reproduction `scripts/divergences/reproductions/fork_warning_pid_prefix.toml` holds this record and run with both engines' warning lines, the process id masked. In the format comparator's warnings mode the class, named in `--masks`, removes a leading `<pid> : ` from each of the reference's lines before the two files are compared and counts the lines it changed under `masks_applied`. The golden datasets record each reference warning in their manifests as the condition alone, so the harness compares vep-rs's lines against the serial shape; the committed `expected/*.warnings.log` files keep Ensembl VEP's text as written. The comparison oracle's three warnings cases raise no warning inside a child, so none reaches the class.
