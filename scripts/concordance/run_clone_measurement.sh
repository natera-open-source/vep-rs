#!/usr/bin/env bash
#
# run_clone_measurement.sh -- the measurement entrypoint for a single engine.
#
# Runs the full 10-dataset measurement (6 SNP/indel + 4 structural-variant) for ONE
# engine (vep-rs | fastvep | perl) against an already-populated data tree,
# emitting a concordance report per cell (the tuple comparator for the SNP/indel
# cells, tools/vep-sv-compare for the structural-variant ones) and a
# machine-readable wall_times.csv.
#
# Scope: measurement only. It does not fetch inputs, build binaries, or manage
# hosts. Point it at a data tree you have already populated, and give it the
# engine binary path or the Ensembl VEP Docker image reference.
#
# The exact command lines this script embeds:
#
#   vep-rs   SNP/indel: vep -i <canonical.vcf.gz> -o <out> --offline \
#                  --json_cache <cache dir> --fasta <fa> --fork 16 \
#                  --buffer_size 5000 --force --no_stats --quiet
#   vep-rs   SV: same, per-VCF loop, same --fork as SNP/indel,
#                  plus --species <species> --assembly <ASM>
#   fastvep  SNP/indel: RAYON_NUM_THREADS=16 fastvep annotate --gff3 <g> \
#                  --fasta <fa> --transcript-cache <cache.bin> \
#                  -i <vcf> -o <raw> --output-format tab
#                then vep-fastvep-normalize <raw> <output.txt>
#                then the tuple comparator ... --glob 'output.txt'
#   fastvep  SV: same annotate form, per-VCF loop; concordance is best-effort
#                  and is skipped when the loop produces no non-empty output.
#   perl   both: docker run ... vep -i <in> -o <out> --offline --cache \
#                  --dir_cache /work/cache --cache_version <ver> \
#                  --species <species> --assembly <asm> --format vcf \
#                  --buffer_size 5000 --fork 16 --fasta /work/ref.fa \
#                  --no_headers --no_stats --quiet --force_overwrite
#                Ensembl VEP is the ground truth, so it has no concordance branch and
#                its scratch output is deleted. It contributes wall time only.
#
# Expected --data-dir layout. Every path below is overridable by the matching
# VEP_* environment variable; these are only the defaults relative to --data-dir.
#
#   vcf/<asm>/all_variants_canonical/<base>.canonical.vcf.gz   SNP/indel inputs
#   ground_truth/perl/snp_indel/<dataset>/output.txt             SNP/indel truth
#   ground_truth/perl/sv_per_vcf/canonical_inputs/<asm>/       SV inputs
#   ground_truth/perl/sv_per_vcf/<asm>/                        SV truth
#   caches/vep-rs/<asm>/                                       vep-rs JSON cache
#   caches/fastvep/<asm>.bin                                   fastVEP cache
#   caches/perl/vep-cache/<species>/<version>_<asm>/           Perl Storable cache
#   reference/<asm>/<genome>.fa[.fai]                           reference FASTA
#   reference/<asm>/<annotation>.gff3                           GFF3 annotation
#
# Every dataset row names its species beside its assembly. A genome's FASTA, vep-rs JSON cache
# and Ensembl cache version are its row of harness_genomes.tsv beside this script
# (VEP_HARNESS_GENOMES names another table; harness_genomes.sh reads it), the FASTA and the
# cache as paths under --data-dir; the two human rows read the defaults above and 116. The
# fastVEP cache and its GFF3 are named for the two human assemblies alone.
#
# "canonical" means Ensembl-style contig names. Raw all_variants/ inputs score
# F1=0 on the GRCh38 chr21 cells, whose contigs are `chr21` while the ground
# truth and the caches expect `21`; canonical inputs fix that symmetrically for
# every engine. Generate them with `tools/vep-prepare-benchmark-vcfs`; the
# vep-rs and fastVEP engines read them as given, and the Ensembl VEP engine derives its
# own canonical, sites-only copy from all_variants/ with bcftools at run time.
#
# The two ground-truth trees are WRITTEN by generate_reference_output.sh in this
# directory, which runs the Ensembl VEP invocation above over the same inputs (SV at the
# ground truth's own --fork 4). This script only reads them.
#
# Usage:
#   scripts/concordance/run_clone_measurement.sh \
#     --engine vep-rs --vep-binary /path/to/vep \
#     --data-dir <data> --output-dir <data>/work \
#     [--datasets all|snp-indel|sv|"s01 s04" | --cells s01,s02,s07,s08,s09,s10] [--walltime-only] \
#     [--skip-fields] [--fork 16] \
#     [--reference-release 115.2|116.2|116.3] [--run-id <id>] [--host-id <id>] \
#     [--vep-extra-args "<flags>"]
#
#   scripts/concordance/run_clone_measurement.sh --engine fastvep \
#     --fastvep-binary /path/to/fastvep --data-dir <data> --output-dir <data>/work \
#     [--fastvep-extra-args "<flags>"]
#
#   scripts/concordance/run_clone_measurement.sh --engine perl \
#     --perl-image ensemblorg/ensembl-vep:release_116.3 [--perl-image-digest <sha256:hex>] \
#     --data-dir <data> --output-dir <data>/work
set -euo pipefail

# Locate the repo so the builds under tools/ resolve relative to THIS file
# (this script lives at scripts/concordance/, so the repo root is ../..).
# The harness scripts ship beside this one in the same checkout, so nothing else
# has to be located.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
# shellcheck source=scripts/concordance/harness_genomes.sh
source "$SCRIPT_DIR/harness_genomes.sh"

# Defaults
ENGINE=""
VEP_BINARY=""
FASTVEP_BINARY=""
PERL_IMAGE="ensemblorg/ensembl-vep:release_116.3"
# The tag's manifest-list digest on Docker Hub. The Ensembl VEP engine refuses to time an image
# whose RepoDigests lack it, so a moved tag cannot be measured under the pinned engine_version.
PERL_IMAGE_DIGEST="${VEP_PERL_IMAGE_DIGEST:-sha256:2a46b0ba4ac544a87136aad92af92050ccbbf09098c1ce306d7e88ef7ed19aa3}"
DATA_DIR=""
OUTPUT_DIR=""
DATASETS_FILTER="all"
WALLTIME_ONLY="false"
SKIP_FIELDS="false"
# Default parallelism. 16 matches the --fork value the reported measurements
# used; override with --fork to match your own core count.
FORK="${VEP_FORK:-16}"
BUFFER_SIZE="${VEP_BUFFER_SIZE:-5000}"
RUN_ID="manual"
HOST_ID="local"
REPLICATE=""
CONTAINER_PROBE=0
VEP_EXTRA_ARGS=""
FASTVEP_EXTRA_ARGS=""
# The Ensembl VEP release of the ground truth, handed to both comparators as
# --reference-release. Empty means derived per cell from the ground-truth directory's
# provenance (provenance.json or <output>.provenance.json, reference_release.py); a
# cell whose ground truth names no release fails rather than scoring under a guess.
REFERENCE_RELEASE="${VEP_REFERENCE_RELEASE:-}"

usage() {
	cat <<'EOF'
Run the full 10-dataset measurement for ONE engine on a populated data tree.

Required:
  --engine {vep-rs|fastvep|perl}
  --data-dir <path>      Populated data tree (vcf/, caches/, reference/, ground_truth/)
  --output-dir <path>    Where to write per-dataset results + wall_times.csv

Engine binary (one, matching --engine):
  --vep-binary <path>      vep-rs binary             (engine vep-rs)
  --fastvep-binary <path>  fastVEP binary            (engine fastvep). A name without
                           a slash is found through PATH, as the cells start it,
                           never in the working directory. Without
                           --fastvep-extra-args a path holding '=' is an error:
                           env, which starts the binary, reads such a path as an
                           assignment.
  --perl-image <ref>       Ensembl VEP Docker image   (engine perl; default
                           ensemblorg/ensembl-vep:release_116.3). The cache
                           release the Perl cells read is each genome's row of
                           the genome table (116 for the two human rows) or,
                           for every genome, VEP_PERL_CACHE_VERSION; it must
                           match the image's release: a release_115.2 image
                           needs VEP_PERL_CACHE_VERSION=115.
                           The Perl cells count and rewrite whole inputs and
                           outputs with this tree's release build
                           tools/vep-file-census/target/release/vep-file-census
                           (cargo build --release --locked --manifest-path
                           tools/vep-file-census/Cargo.toml); a Perl run without
                           it stops before the first engine run.
  --perl-image-digest <sha256:hex>  Digest the pulled Perl image must carry (engine
                           perl; default the pinned tag's manifest-list digest, or
                           VEP_PERL_IMAGE_DIGEST)

Optional:
  --datasets <filter>    all (default) | snp-indel | sv | space-separated "s01 s04"
  --cells <ids>          A comma-separated cell list (s01,s02,s07,s08,s09,s10), the
                         comma form of --datasets; an unknown or repeated id is an
                         error before any cell runs
  --walltime-only        Skip concordance compare / GT-sync / normalize; time
                         the annotation only (annotation output still written so
                         the real I/O happens). Default: full concordance + timing.
  --skip-fields          Skip the per-column field comparator (report/fields.json)
                         on the vep-rs SNP/indel cells; F1 is computed either way.
                         For runs that need only F1 and wall time. Default: run it.
                         The field comparator is this tree's release build,
                         tools/vep-fields-compare/target/release/vep-fields-compare,
                         and the SNP/indel tuple comparator of every cell that
                         scores F1 is
                         tools/vep-outputs-compare/target/release/vep-outputs-compare
                         (each: cargo build --release --locked --manifest-path
                         tools/<tool>/Cargo.toml). A comparator the run would
                         start that is not an executable stops the run before
                         the first engine run; nothing is built here.
  --fork <N>             SNP/indel parallelism (default 16; fastVEP uses
                         RAYON_NUM_THREADS=N)
  --buffer-size <N>      Engine buffer size (default 5000)
  --reference-release <115.2|116.2|116.3>
                         The Ensembl VEP release the ground truth was produced with,
                         passed to the tuple comparator and vep-sv-compare
                         as --reference-release. Default: read per cell from the
                         ground-truth directory's provenance (the image tag, else the
                         cache version); a cell with neither fails. Also VEP_REFERENCE_RELEASE.
  --run-id <id>          Provenance: identifier for this measurement campaign,
                         written to wall_times.csv
  --host-id <id>         Provenance: identifier for the machine running this
                         invocation, written to wall_times.csv
  --replicate <n>        Provenance: 1..N index of this run within the campaign
                         (defaults to the host id on a single run)
  --container-probe <N>  Perl engine only: after the datasets, time N bare starts of
                         the Perl image (`docker run --rm <image> true`, no
                         annotation) under /usr/bin/time -v and write them to
                         <output-dir>/container_start.csv (default 0: none). Every
                         Perl cell pays one container start per invocation; the
                         probe measures that cost on its own.
  --vep-extra-args "<flags>"
                         vep-rs engine only: arguments appended, after the harness's
                         own, to every vep-rs invocation of the run (the SNP/indel
                         cells and the SV per-VCF loop alike), split on whitespace.
                         They reach the binary through <output-dir>/vep-wrapper.sh,
                         so the invocations documented in this file's header stay
                         as written; the wrapper costs one bash start per invocation
                         inside the timed window, which a run without this flag does
                         not pay. The wrapper logs every argument vector it runs to
                         <output-dir>/vep-wrapper.argv; after the last cell the log is
                         folded into <output-dir>/provenance.json beside the binary,
                         its sha256 and the extra arguments (that record is written
                         for every vep-rs run, extra arguments or not). An empty value
                         appends nothing and writes no wrapper.
  --fastvep-extra-args "<flags>"
                         fastVEP engine only: arguments appended, after the harness's
                         own, to every fastVEP `annotate` run (the warmup and the timed
                         run of each SNP/indel cell, and each input of the SV per-VCF
                         loop), split on whitespace. Each word is an option (-x, --name,
                         --name=value) or a value of [A-Za-z0-9_=.,:/-] that opens with
                         no dash; any other word stops the run before anything is
                         written. They reach the binary through
                         <output-dir>/fastvep-wrapper.sh, at the cost of one bash start
                         per invocation inside the timed window, and the wrapper logs
                         every argument vector it runs to
                         <output-dir>/fastvep-wrapper.argv. <output-dir>/provenance.json
                         then names the binary, the words (fastvep_extra_args, joined by
                         single spaces, and fastvep_extra_argv), the wrapper
                         (fastvep_wrapper) and, after the last cell, the executed
                         argument vectors (fastvep_executed_argv). An empty value
                         appends nothing and writes no wrapper and no record.
                         With words, an <output-dir> whose path holds '=' is an
                         error before anything is written: env, which starts the
                         wrapper, reads such a path as an assignment.
                         --sa-dir names a directory of one assembly's annotation
                         databases, so a run whose words carry it takes cells of one
                         assembly: cells of two are an error before any cell runs.
EOF
	exit 1
}

while [[ $# -gt 0 ]]; do
	case "$1" in
	--engine)
		ENGINE="$2"
		shift 2
		;;
	--vep-binary)
		VEP_BINARY="$2"
		shift 2
		;;
	--fastvep-binary)
		FASTVEP_BINARY="$2"
		shift 2
		;;
	--perl-image)
		PERL_IMAGE="$2"
		shift 2
		;;
	--perl-image-digest)
		PERL_IMAGE_DIGEST="$2"
		shift 2
		;;
	--data-dir)
		DATA_DIR="$2"
		shift 2
		;;
	--output-dir)
		OUTPUT_DIR="$2"
		shift 2
		;;
	--datasets)
		DATASETS_FILTER="$2"
		shift 2
		;;
	--cells)
		DATASETS_FILTER="${2//,/ }"
		shift 2
		;;
	--walltime-only)
		WALLTIME_ONLY="true"
		shift
		;;
	--skip-fields)
		SKIP_FIELDS="true"
		shift
		;;
	--fork)
		FORK="$2"
		shift 2
		;;
	--buffer-size)
		BUFFER_SIZE="$2"
		shift 2
		;;
	--reference-release)
		REFERENCE_RELEASE="$2"
		shift 2
		;;
	--run-id)
		RUN_ID="$2"
		shift 2
		;;
	--host-id)
		HOST_ID="$2"
		shift 2
		;;
	--replicate)
		REPLICATE="$2"
		shift 2
		;;
	--container-probe)
		CONTAINER_PROBE="$2"
		shift 2
		;;
	--vep-extra-args)
		VEP_EXTRA_ARGS="$2"
		shift 2
		;;
	--fastvep-extra-args)
		FASTVEP_EXTRA_ARGS="$2"
		shift 2
		;;
	-h | --help) usage ;;
	*)
		echo "ERROR: unknown argument: $1" >&2
		usage
		;;
	esac
done

# should_run matches a cell id between single spaces while the validation loop below splits on
# any whitespace, so a tab- or newline-separated list would validate and then skip every cell:
# one space per gap, no padding.
DATASETS_FILTER=$(printf '%s' "$DATASETS_FILTER" | tr -s '[:space:]' ' ')
DATASETS_FILTER="${DATASETS_FILTER# }"
DATASETS_FILTER="${DATASETS_FILTER% }"

# Validate args
case "$REFERENCE_RELEASE" in
"" | 115.2 | 116.2 | 116.3) ;;
*)
	echo "ERROR: --reference-release must be 115.2, 116.2 or 116.3 (got '$REFERENCE_RELEASE')" >&2
	usage
	;;
esac
[[ "$CONTAINER_PROBE" =~ ^[0-9]+$ ]] || {
	echo "ERROR: --container-probe takes a non-negative integer (got '$CONTAINER_PROBE')" >&2
	usage
}
[[ "$CONTAINER_PROBE" -eq 0 || "$ENGINE" == "perl" ]] || {
	echo "ERROR: --container-probe applies to --engine perl only (got engine '$ENGINE')" >&2
	usage
}
[[ -z "${VEP_EXTRA_ARGS//[[:space:]]/}" || "$ENGINE" == "vep-rs" ]] || {
	echo "ERROR: --vep-extra-args applies to --engine vep-rs only (got engine '$ENGINE')" >&2
	usage
}
# A word of --fastvep-extra-args is an option (one or two dashes, then a letter) or a value that
# opens with no dash; the class holds no character the shell or JSON reads, so a word needs no
# escaping wherever it is rendered. `read -d ''` splits the whole value on whitespace, a newline
# included, without globbing a word. The refusals below come before the output directory is made.
FASTVEP_EXTRA_OPTION='^-{1,2}[A-Za-z][A-Za-z0-9_=.,:/-]*$'
FASTVEP_EXTRA_VALUE='^[A-Za-z0-9_=.,:/][A-Za-z0-9_=.,:/-]*$'
FASTVEP_EXTRA_ARGV=()
read -r -d '' -a FASTVEP_EXTRA_ARGV <<<"$FASTVEP_EXTRA_ARGS" || true
[[ ${#FASTVEP_EXTRA_ARGV[@]} -eq 0 || "$ENGINE" == "fastvep" ]] || {
	echo "ERROR: --fastvep-extra-args applies to --engine fastvep only (got engine '$ENGINE')" >&2
	usage
}
for word in ${FASTVEP_EXTRA_ARGV[@]+"${FASTVEP_EXTRA_ARGV[@]}"}; do
	# Matched in the C locale: under a collating locale `[A-Za-z]` and `[0-9]` also match
	# letters and digits outside ASCII.
	(LC_ALL=C && [[ "$word" =~ $FASTVEP_EXTRA_OPTION || "$word" =~ $FASTVEP_EXTRA_VALUE ]]) || {
		echo "ERROR: [run_clone_measurement] --fastvep-extra-args word '$word' is neither an option (-x, --name, --name=value) nor a value of [A-Za-z0-9_=.,:/-] that opens with no dash" >&2
		exit 1
	}
done
# The SNP/indel cells run fastVEP as the command operand of `env`, which reads an operand holding
# `=` as an assignment and runs the next word instead. The operand is
# <output-dir>/fastvep-wrapper.sh with words and the binary's own path without them; the wrapper
# execs the binary, so with words the binary's path never reaches `env`.
[[ ${#FASTVEP_EXTRA_ARGV[@]} -eq 0 || "$OUTPUT_DIR" != *=* ]] || {
	echo "ERROR: [run_clone_measurement] --output-dir '$OUTPUT_DIR' holds '=': with --fastvep-extra-args the SNP/indel cells run <output-dir>/fastvep-wrapper.sh under env, which reads a path holding '=' as an assignment; name an output directory without one" >&2
	exit 1
}
[[ "$ENGINE" != "fastvep" || ${#FASTVEP_EXTRA_ARGV[@]} -gt 0 || "$FASTVEP_BINARY" != *=* ]] || {
	echo "ERROR: [run_clone_measurement] --fastvep-binary '$FASTVEP_BINARY' holds '=': the SNP/indel cells run it under env, which reads a path holding '=' as an assignment; name a path without one" >&2
	exit 1
}
case "$ENGINE" in
vep-rs | fastvep | perl) ;;
"")
	echo "ERROR: --engine required" >&2
	usage
	;;
*)
	echo "ERROR: invalid --engine: $ENGINE (vep-rs|fastvep|perl)" >&2
	usage
	;;
esac
[[ -z "$DATA_DIR" ]] && {
	echo "ERROR: --data-dir required" >&2
	usage
}
[[ -z "$OUTPUT_DIR" ]] && {
	echo "ERROR: --output-dir required" >&2
	usage
}
[[ -d "$DATA_DIR" ]] || {
	echo "ERROR: data-dir not a directory: $DATA_DIR" >&2
	exit 1
}
case "$ENGINE" in
vep-rs)
	[[ -n "$VEP_BINARY" ]] || {
		echo "ERROR: --vep-binary required for engine vep-rs" >&2
		exit 1
	}
	# A readable regular file, not only `-x`: a directory and an execute-only file both pass
	# `-x`, and provenance.json reads the file to hash it.
	[[ -f "$VEP_BINARY" && -r "$VEP_BINARY" && -x "$VEP_BINARY" ]] || {
		echo "ERROR: vep-rs binary not a readable executable file: $VEP_BINARY" >&2
		exit 1
	}
	# The command lines run the operand as given, so a name without a slash would resolve
	# through PATH while the test above and the hash read the file in the working directory;
	# the absolute path makes the three agree.
	[[ "$VEP_BINARY" == /* ]] || VEP_BINARY="$PWD/$VEP_BINARY"
	command -v python3 >/dev/null 2>&1 || {
		echo "ERROR: [run_clone_measurement] python3 is not on PATH: every vep-rs run writes <output-dir>/provenance.json through it; install it" >&2
		exit 1
	}
	;;
fastvep)
	[[ -n "$FASTVEP_BINARY" ]] || {
		echo "ERROR: --fastvep-binary required for engine fastvep" >&2
		exit 1
	}
	# The cells start the binary by the word given, so a name without a slash runs the file
	# PATH resolves, never a like-named file of the working directory, and the executable test
	# reads that file. `type -P` returns a non-executable PATH match when no executable one
	# exists, hence `-x` on its result.
	if [[ "$FASTVEP_BINARY" == */* ]]; then
		[[ -x "$FASTVEP_BINARY" ]] || {
			echo "ERROR: fastVEP binary not executable: $FASTVEP_BINARY" >&2
			exit 1
		}
	else
		fastvep_on_path=$(type -P -- "$FASTVEP_BINARY") || fastvep_on_path=""
		[[ -x "$fastvep_on_path" ]] || {
			echo "ERROR: [run_clone_measurement] --fastvep-binary '$FASTVEP_BINARY' names no executable on PATH: a name without a slash is started through PATH, never from the working directory; give a path (./$FASTVEP_BINARY for a file there)" >&2
			exit 1
		}
	fi
	[[ ${#FASTVEP_EXTRA_ARGV[@]} -eq 0 ]] || command -v python3 >/dev/null 2>&1 || {
		echo "ERROR: [run_clone_measurement] python3 is not on PATH: a fastVEP run with --fastvep-extra-args writes <output-dir>/provenance.json through it; install it" >&2
		exit 1
	}
	;;
perl)
	command -v docker >/dev/null 2>&1 || {
		echo "ERROR: docker required for engine perl" >&2
		exit 1
	}
	[[ "$PERL_IMAGE_DIGEST" =~ ^sha256:[0-9a-f]{64}$ ]] || {
		echo "ERROR: --perl-image-digest must be sha256:<64 hex> (got '$PERL_IMAGE_DIGEST')" >&2
		exit 1
	}
	docker image inspect "$PERL_IMAGE" >/dev/null 2>&1 || docker pull "$PERL_IMAGE" >/dev/null
	perl_digests=$(docker image inspect --format '{{join .RepoDigests ","}}' "$PERL_IMAGE" 2>/dev/null || echo "")
	[[ "$perl_digests" == *"@$PERL_IMAGE_DIGEST"* ]] || {
		echo "ERROR: [perl-image] digest mismatch: ${perl_digests:-none} (expected $PERL_IMAGE_DIGEST)" >&2
		exit 1
	}
	;;
esac

echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] run_clone_measurement engine=$ENGINE datasets=$DATASETS_FILTER walltime_only=$WALLTIME_ONLY skip_fields=$SKIP_FIELDS fork=$FORK"
echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] data_dir=$DATA_DIR output_dir=$OUTPUT_DIR repo=$REPO_DIR"

# Resource locations, all relative to DATA_DIR by default and all overridable by
# the matching environment variable. Override any one of them when your data
# tree is laid out differently; see the layout note in the file header.
#
# The ground-truth trees must hold Ensembl VEP output produced with the SAME cache
# release and reference FASTA the engines are given here, or the concordance
# figures compare two different annotations.
SNP_INDEL_GT_DIR="${VEP_SNP_INDEL_GT_DIR:-$DATA_DIR/ground_truth/perl/snp_indel}"
SV_GT_DIR="${VEP_SV_GT_DIR:-$DATA_DIR/ground_truth/perl/sv_per_vcf}"

# The genome table: one row a genome, read through harness_genome_field. A genome's FASTA
# and vep-rs JSON cache are the row's `fasta` and `json_cache`, two paths under DATA_DIR,
# and its Ensembl cache version the row's `cache_version`.
GENOMES_TABLE="${VEP_HARNESS_GENOMES:-$SCRIPT_DIR/harness_genomes.tsv}"
# genome_table_lines <column> <prefix>: one "<species>:<assembly>	<prefix><field>" line for
# every row of the genome table whose <column> field is not empty, in the table's order.
# Each field is read through harness_genome_field, so a table that reader refuses is
# refused here with the reader's own line; a table that cannot be read, or that names no
# genome under a species and an assembly column, returns 1 as well.
genome_table_lines() {
	local column="$1" prefix="$2" keys key field
	[[ -f "$GENOMES_TABLE" && -r "$GENOMES_TABLE" ]] || {
		echo "ERROR: [run_clone_measurement] the genome table cannot be read: $GENOMES_TABLE" >&2
		return 1
	}
	keys=$(awk -F'\t' '
		NR == 1 {
			for (i = 1; i <= NF; i++) {
				if ($i == "species") species_at = i
				if ($i == "assembly") assembly_at = i
			}
			next
		}
		species_at && assembly_at { print $species_at ":" $assembly_at }
	' <"$GENOMES_TABLE")
	[[ -n "$keys" ]] || {
		echo "ERROR: [run_clone_measurement] the genome table names no genome under a species and an assembly column: $GENOMES_TABLE" >&2
		return 1
	}
	while IFS= read -r key; do
		field=$(harness_genome_field "$GENOMES_TABLE" "${key%%:*}" "${key#*:}" "$column") || return 1
		[[ -z "$field" ]] || printf '%s\t%s%s\n' "$key" "$prefix" "$field"
	done <<<"$keys"
}
# resource_for <table> <species> <assembly> <resource> <remedy>: the table's path for the
# genome; a table is "<species>:<assembly>	<path>" lines, the last line for a genome winning.
# A genome with no line is named on stderr with the resource it lacks and the remedy and
# returns 1, so the cell fails before any engine runs with an empty path.
resource_for() {
	local path
	path=$(printf '%s\n' "$1" | awk -F'\t' -v k="$2:$3" '$1 == k {v = $2} END {print v}')
	[[ -n "$path" ]] || {
		echo "ERROR: [run_clone_measurement] no $4 is named for the genome $2:$3; $5" >&2
		return 1
	}
	printf '%s' "$path"
}
# Per-genome reference FASTA (shared across all engines), indexed (.fa + .fai): the genome
# table's rows, then VEP_FASTA_GRCH37 and VEP_FASTA_GRCH38, which replace the two human rows.
FASTA_TABLE=$(genome_table_lines fasta "$DATA_DIR/") || exit 1
[[ -z "${VEP_FASTA_GRCH37:-}" ]] || FASTA_TABLE="$FASTA_TABLE
homo_sapiens:GRCh37	$VEP_FASTA_GRCH37"
[[ -z "${VEP_FASTA_GRCH38:-}" ]] || FASTA_TABLE="$FASTA_TABLE
homo_sapiens:GRCh38	$VEP_FASTA_GRCH38"
fasta_for() { resource_for "$FASTA_TABLE" "$1" "$2" "FASTA" "add its row to $GENOMES_TABLE"; }

# Engine-specific resources.
# vep-rs JSON caches. Build them with `vep-cache-builder`, or convert an existing
# Perl Storable cache with `tools/vep-storable-to-json`. The cache must carry populated
# 3' UTR sequences; a cache built without them scores measurably lower. The genome
# table's rows name them, then VEP_RS_CACHE_GRCH37 and VEP_RS_CACHE_GRCH38, which replace
# the two human rows.
VEPRS_CACHE_TABLE=$(genome_table_lines json_cache "$DATA_DIR/") || exit 1
[[ -z "${VEP_RS_CACHE_GRCH37:-}" ]] || VEPRS_CACHE_TABLE="$VEPRS_CACHE_TABLE
homo_sapiens:GRCh37	$VEP_RS_CACHE_GRCH37"
[[ -z "${VEP_RS_CACHE_GRCH38:-}" ]] || VEPRS_CACHE_TABLE="$VEPRS_CACHE_TABLE
homo_sapiens:GRCh38	$VEP_RS_CACHE_GRCH38"
# fastVEP transcript cache (.bin) and the GFF3 passed beside it. fastVEP loads the cache alone
# when its file exists and does not parse the GFF3, so the cache is the gene model source of
# every fastVEP cell, whatever gene models the GFF3 holds.
FASTVEP_CACHE_GRCH37="${VEP_FASTVEP_CACHE_GRCH37:-$DATA_DIR/caches/fastvep/grch37.bin}"
FASTVEP_CACHE_GRCH38="${VEP_FASTVEP_CACHE_GRCH38:-$DATA_DIR/caches/fastvep/grch38.bin}"
GFF3_GRCH37="${VEP_GFF3_GRCH37:-$DATA_DIR/reference/grch37/annotation.gff3}"
GFF3_GRCH38="${VEP_GFF3_GRCH38:-$DATA_DIR/reference/grch38/annotation.gff3}"
# The comparators of the SNP/indel cells, both release builds of this tree: the tuple
# comparator of the cells that score F1, tools/vep-outputs-compare, and the per-column field
# comparator of the vep-rs cells, tools/vep-fields-compare. Nothing is built here and no other
# program runs in a comparator's place: one the run would start that is not an executable
# file stops the run here, before any engine runs. A run that scores no F1 (--walltime-only,
# the Ensembl VEP engine) asks for neither, and the field comparator is asked for only where it
# runs: the vep-rs engine without --skip-fields.
#
# The fastVEP cells that score F1 run a third release build of this tree first,
# tools/vep-fastvep-normalize, found and refused the same way. It corrects the
# `Feature_type=Transcript` that fastVEP through v0.4.0 writes on intergenic rows, where
# Ensembl VEP and vep-rs both emit `-`; feature type is part of the comparison key, so
# without it every intergenic row of such a run scores discordant. From v0.5.0 fastVEP
# writes `-` there itself, and the rewrite changes no row of its output.
OUTPUTS_COMPARE=()
FIELDS_COMPARE=()
FASTVEP_NORMALIZE=()
if [[ "$WALLTIME_ONLY" != "true" && "$ENGINE" != "perl" ]]; then
	outputs_bin="$REPO_DIR/tools/vep-outputs-compare/target/release/vep-outputs-compare"
	if [[ -f "$outputs_bin" && -x "$outputs_bin" ]]; then
		OUTPUTS_COMPARE=("$outputs_bin")
	else
		echo "ERROR: [run_clone_measurement] the tuple comparator $outputs_bin is not an executable: build the tool with 'cargo build --release --locked --manifest-path $REPO_DIR/tools/vep-outputs-compare/Cargo.toml'" >&2
		exit 1
	fi
	if [[ "$ENGINE" == "vep-rs" && "$SKIP_FIELDS" != "true" ]]; then
		fields_bin="$REPO_DIR/tools/vep-fields-compare/target/release/vep-fields-compare"
		if [[ -f "$fields_bin" && -x "$fields_bin" ]]; then
			FIELDS_COMPARE=("$fields_bin")
		else
			echo "ERROR: [run_clone_measurement] the field comparator $fields_bin is not an executable: build the tool with 'cargo build --release --locked --manifest-path $REPO_DIR/tools/vep-fields-compare/Cargo.toml'" >&2
			exit 1
		fi
	fi
	if [[ "$ENGINE" == "fastvep" ]]; then
		normalize_bin="$REPO_DIR/tools/vep-fastvep-normalize/target/release/vep-fastvep-normalize"
		if [[ -f "$normalize_bin" && -x "$normalize_bin" ]]; then
			FASTVEP_NORMALIZE=("$normalize_bin")
		else
			echo "ERROR: [run_clone_measurement] the fastVEP output normaliser $normalize_bin is not an executable: build the tool with 'cargo build --release --locked --manifest-path $REPO_DIR/tools/vep-fastvep-normalize/Cargo.toml'" >&2
			exit 1
		fi
	fi
fi
# The Ensembl VEP cells make passes of their own over whole inputs and outputs: the contig-prefix
# pass of `canonicalize` and the row and record counts of `assert_perl_annotated`. They run
# a fourth release build of this tree, tools/vep-file-census, across the host's cores; it
# is found and refused the same way, and only the Ensembl VEP engine asks for it.
FILE_CENSUS=()
if [[ "$ENGINE" == "perl" ]]; then
	census_bin="$REPO_DIR/tools/vep-file-census/target/release/vep-file-census"
	if [[ -f "$census_bin" && -x "$census_bin" ]]; then
		FILE_CENSUS=("$census_bin")
	else
		echo "ERROR: [run_clone_measurement] the file census $census_bin is not an executable: build the tool with 'cargo build --release --locked --manifest-path $REPO_DIR/tools/vep-file-census/Cargo.toml'" >&2
		exit 1
	fi
fi
# Perl Storable cache (the --dir_cache mount; contains <species>/...). VEP resolves
# <species>/<cache_version>_<assembly> underneath it from --species, --cache_version and
# --assembly.
#
# The flag is `--dir_cache`. `--full_cache_dir` is not a VEP flag, and VEP does not
# error on it: it ignores it, falls back to
# a default cache directory that holds nothing, finds no annotation source, and emits one
# warning per record instead of annotations, so the run produces zero data lines and exits 0.
#
# Nothing downstream would notice: the timed Ensembl VEP run's output is written to a scratch
# directory and deleted, and concordance is scored against a pinned ground truth rather than
# this run's output, so an Ensembl VEP wall time can measure a VEP that annotated nothing. The
# `assert_perl_annotated` gate below exists to catch that, and `assert_perl_cache_complete`
# fails any cell whose cache is a chromosome subset.
PERL_CACHE_DIR="${VEP_PERL_CACHE_DIR:-$DATA_DIR/caches/perl/vep-cache}"
# Ensembl VEP cache release. Must match the release the ground truth was built with.
# VEP_PERL_CACHE_VERSION names one release for every genome; without it each genome reads
# the cache_version of its row of the genome table, since an Ensembl Genomes cache carries
# its own release number, and a genome with no row reads 116.
PERL_CACHE_VERSION="${VEP_PERL_CACHE_VERSION:-116}"
CACHE_VERSION_TABLE=""
if [[ -z "${VEP_PERL_CACHE_VERSION:-}" ]]; then
	CACHE_VERSION_TABLE=$(genome_table_lines cache_version "") || exit 1
	# A table that names no cache_version leaves one empty line here, which is no genome.
	# The table is awk's input as a here-string: a pipe's writer that is still writing
	# when its reader leaves at its first match is ended by the closed pipe, which
	# pipefail makes fatal.
	bad_version=$(awk -F'\t' 'NF && $2 !~ /^[1-9][0-9]*$/ { print $1 " " $2; exit }' <<<"$CACHE_VERSION_TABLE")
	[[ -z "$bad_version" ]] || {
		echo "ERROR: [run_clone_measurement] the genome table gives ${bad_version%% *} the cache_version '${bad_version#* }', which is no positive integer: $GENOMES_TABLE" >&2
		exit 1
	}
fi
# cache_version_for <species> <assembly>: the Ensembl cache release the genome's Ensembl VEP cells
# read: its line of CACHE_VERSION_TABLE, or PERL_CACHE_VERSION for a genome with no line.
cache_version_for() {
	local version
	version=$(printf '%s\n' "$CACHE_VERSION_TABLE" | awk -F'\t' -v k="$1:$2" '$1 == k {v = $2} END {print v}')
	printf '%s' "${version:-$PERL_CACHE_VERSION}"
}

# fastVEP is measured on the two human assemblies alone: its transcript cache and the GFF3
# passed beside it are named for those two genomes, and a fastVEP cell of any other
# genome fails as genome_resources_missing.
FASTVEP_CACHE_TABLE="homo_sapiens:GRCh37	$FASTVEP_CACHE_GRCH37
homo_sapiens:GRCh38	$FASTVEP_CACHE_GRCH38"
GFF3_TABLE="homo_sapiens:GRCh37	$GFF3_GRCH37
homo_sapiens:GRCh38	$GFF3_GRCH38"
# Each takes <species> <assembly>.
veprs_cache_for() { resource_for "$VEPRS_CACHE_TABLE" "$1" "$2" "vep-rs JSON cache" "add its row to $GENOMES_TABLE"; }
fastvep_cache_for() { resource_for "$FASTVEP_CACHE_TABLE" "$1" "$2" "fastVEP cache" "fastVEP is measured on the two human assemblies alone"; }
gff3_for() { resource_for "$GFF3_TABLE" "$1" "$2" "GFF3" "fastVEP is measured on the two human assemblies alone"; }

# Dataset table: id:name:assembly:type:vcf_base:species
declare -a DATASETS=(
	"s01:ClinVar_GRCh37:GRCh37:snp_indel:clinvar_grch37:homo_sapiens"
	"s02:ClinVar_GRCh38:GRCh38:snp_indel:clinvar_grch38:homo_sapiens"
	"s03:gnomAD_v2.1.1_chr21:GRCh37:snp_indel:gnomad_genomes_v2.1.1_chr21:homo_sapiens"
	"s04:gnomAD_v4.1_chr21:GRCh38:snp_indel:gnomad_genomes_v4.1_chr21:homo_sapiens"
	# s05 + s06 are the ONLY datasets whose SOURCE inputs carry per-sample genotypes
	# (1000 Genomes; the high-coverage set has 3,202 samples, FORMAT=GT). Genotypes
	# do NOT affect VEP annotation output (VEP keys on position/allele; tab output has
	# no genotype columns), so the apples-to-apples cells (s05/s06) use the SITES-ONLY
	# inputs (bcftools view -G; the 4 other SNP/indel datasets s01-s04 are sites-only at
	# source). The `*m` rows below measure the SAME datasets WITH genotypes for the
	# sample-scaling measurement (vep-rs + fastVEP; this script's Ensembl VEP branch always
	# strips genotypes, and Ensembl VEP's one genotype-carrying measurement is the probe
	# row of manuscript/data/perl_bench.csv).
	# Concordance F1 is identical between the sites-only and multi-sample forms.
	"s05:1KG_integrated_chr21:GRCh37:snp_indel:1kg_integrated_chr21_sitesonly:homo_sapiens"
	"s06:1KG_highcov_chr21:GRCh38:snp_indel:1kg_highcov_chr21:homo_sapiens"
	"s07:SV_GRCh37:GRCh37:sv:-:homo_sapiens"
	"s08:SV_GRCh38:GRCh38:sv:-:homo_sapiens"
	"s09:SV_real_world_GRCh37:GRCh37:sv:-:homo_sapiens"
	"s10:SV_real_world_GRCh38:GRCh38:sv:-:homo_sapiens"
	# Sample-scaling rows (with-genotypes; NOT part of --datasets all / the headline table).
	# Run explicitly, e.g. --datasets "s05m s06m". vcf_base names the genotype-carrying input.
	"s05m:1KG_integrated_chr21_multisample:GRCh37:snp_indel:1kg_integrated_chr21:homo_sapiens"
	"s06m:1KG_highcov_chr21_multisample:GRCh38:snp_indel:1kg_highcov_chr21_multisample:homo_sapiens"
)

# A cell list names ids from the table above; an unknown id would match nothing, time nothing
# and exit 0, and a repeated id would run its cell once under a list that names it twice, so
# both are refused here, before any cell runs.
case "$DATASETS_FILTER" in
all | snp-indel | sv) ;;
*)
	known_cells=$(printf '%s\n' "${DATASETS[@]}" | cut -d: -f1 | tr '\n' ' ')
	[[ -n "$DATASETS_FILTER" ]] || {
		echo "ERROR: [run_clone_measurement] empty cell list in --datasets/--cells; the cells are: $known_cells" >&2
		exit 1
	}
	seen_cells=""
	# An id shaped like `s0*` is refused as text; with globbing on, a file named `s01` in the
	# working directory would satisfy it.
	set -o noglob
	for cell_id in $DATASETS_FILTER; do
		[[ " $known_cells " == *" $cell_id "* ]] || {
			echo "ERROR: [run_clone_measurement] unknown cell '$cell_id' in --datasets/--cells; the cells are: $known_cells" >&2
			exit 1
		}
		[[ " $seen_cells " != *" $cell_id "* ]] || {
			echo "ERROR: [run_clone_measurement] repeated cell '$cell_id' in --datasets/--cells" >&2
			exit 1
		}
		seen_cells="$seen_cells $cell_id"
	done
	set +o noglob
	;;
esac

# SV per-assembly dedupe: the canonical per-VCF SV GT carries all 16 VCFs per
# assembly in one directory (13 synthetic + 3 real-world). s07/s09 (GRCh37) and
# s08/s10 (GRCh38) measure the SAME per-assembly SV set, so the annotation +
# comparison is done ONCE per assembly under s07 (GRCh37) / s08 (GRCh38);
# s09/s10 record a pointer marker. vep-sv-compare reports overall +
# per-file F1, so the synthetic vs real-world split is recoverable.
# bash 3.2 has no associative arrays, and macOS ships 3.2 while `bash -n`
# and shellcheck both PASS on `declare -A` -- it fails at RUNTIME, mid-run. A
# delimited string is exact here: at most two keys (the two assemblies) and the
# values are dataset ids with no colons.
SV_DONE=""

# should_run <dataset_id> <stype> : honor --datasets filter.
should_run() {
	local sid="$1" stype="$2"
	# Sample-scaling rows (ids ending in `m`, e.g. s05m/s06m) are OPT-IN: they run ONLY
	# when named explicitly (--datasets "s05m s06m"), never under all/snp-indel, so a
	# default run stays the sites-only apples-to-apples set and never times the
	# genotype matrix, whose cost is measured separately.
	case "$DATASETS_FILTER" in
	all) [[ "$sid" != *m ]] ;;
	snp-indel) [[ "$stype" == "snp_indel" && "$sid" != *m ]] ;;
	sv) [[ "$stype" == "sv" ]] ;;
	*) [[ " $DATASETS_FILTER " == *" $sid "* ]] ;;
	esac
}

# fastVEP reads its annotation databases from the one directory --sa-dir names, and a
# directory holds one assembly's, so the cells of a run whose words carry it are of one
# assembly: cells of two would all be annotated from one assembly's databases, and are
# refused here, before any cell runs. `${...[*]-}`, not a bare `"${...[*]}"`: on bash 3.2
# an empty array is an unbound variable under `set -u`.
fastvep_words=" ${FASTVEP_EXTRA_ARGV[*]-} "
if [[ "$fastvep_words" == *" --sa-dir "* || "$fastvep_words" == *" --sa-dir="* ]]; then
	sa_assemblies=""
	sa_count=0
	for entry in "${DATASETS[@]}"; do
		IFS=':' read -r dataset_id _ assembly stype _ <<<"$entry"
		should_run "$dataset_id" "$stype" || continue
		[[ " $sa_assemblies " != *" $assembly "* ]] || continue
		sa_assemblies="$sa_assemblies $assembly"
		sa_count=$((sa_count + 1))
	done
	[[ "$sa_count" -le 1 ]] || {
		echo "ERROR: [run_clone_measurement] --fastvep-extra-args carries --sa-dir, a directory of one assembly's databases, and the selected cells are of $sa_count assemblies:$sa_assemblies; select the cells of one assembly with --datasets or --cells and run once per assembly" >&2
		exit 1
	}
fi

# The comparator of the structural-variant cells that score F1: the release build of
# tools/vep-sv-compare, which nothing here builds, taken from CARGO_TARGET_DIR when that is
# set and holds one (cargo writes a build there, not under the tool's own target directory)
# and from the tool's own target directory otherwise. A run that will score such a cell
# stops here, before any engine runs, when neither holds an executable.
SV_COMPARE="$REPO_DIR/tools/vep-sv-compare/target/release/vep-sv-compare"
if [[ -n "${CARGO_TARGET_DIR:-}" && -f "$CARGO_TARGET_DIR/release/vep-sv-compare" && -x "$CARGO_TARGET_DIR/release/vep-sv-compare" ]]; then
	SV_COMPARE="$CARGO_TARGET_DIR/release/vep-sv-compare"
fi
if [[ "$WALLTIME_ONLY" != "true" && "$ENGINE" != "perl" ]]; then
	for entry in "${DATASETS[@]}"; do
		IFS=':' read -r dataset_id _ _ stype _ <<<"$entry"
		[[ "$stype" == "sv" ]] || continue
		should_run "$dataset_id" "$stype" || continue
		[[ -f "$SV_COMPARE" && -x "$SV_COMPARE" ]] || {
			echo "ERROR: [run_clone_measurement] the structural-variant comparator $SV_COMPARE is not an executable: build the tool with 'cargo build --release --locked --manifest-path $REPO_DIR/tools/vep-sv-compare/Cargo.toml'" >&2
			exit 1
		}
		break
	done
fi

# Output dir + wall_times.csv. The schema is identical across all three engines,
# so runs can be concatenated into one uniform CSV.
#
# The first write under the output directory: every refusal on the arguments stands above
# this line, so a refused run leaves an earlier run's files there as they are.
mkdir -p "$OUTPUT_DIR"
# A re-run into the same output directory starts without the previous run's record, wrappers
# and logs: a vep-rs run writes the record always and its wrapper and log only with extra
# arguments, a fastVEP run writes all three only with extra arguments, so none of them sits
# beside a wall_times.csv it does not describe.
rm -f "$OUTPUT_DIR/provenance.json" "$OUTPUT_DIR/vep-wrapper.sh" "$OUTPUT_DIR/vep-wrapper.argv" \
	"$OUTPUT_DIR/fastvep-wrapper.sh" "$OUTPUT_DIR/fastvep-wrapper.argv"
WALL_CSV="$OUTPUT_DIR/wall_times.csv"
echo "engine,run_id,host_id,replicate,suite_id,suite_name,assembly,run_kind,wall_seconds,exit_code" >"$WALL_CSV"
# replicate is a provenance column: this run's 1..N index across a set of repeated
# runs, passed in by the caller. When unset (a single run) it defaults to the host
# id so the column is never empty.
[[ -n "$REPLICATE" ]] || REPLICATE="$HOST_ID"

# --vep-extra-args reaches the binary through a wrapper that VEP_BINARY is pointed at, so
# the command lines below stay as the header documents them. Without extra arguments no
# wrapper is written, so a timed command line never carries one it was not asked for.
# provenance.json's executed_argv is filled after the dataset loop from the wrapper's log.
VEP_REAL_BINARY="$VEP_BINARY"
VEP_WRAPPER=""
VEP_ARGV_LOG=""
VEP_EXTRA_ARGV=()
if [[ "$ENGINE" == "vep-rs" ]]; then
	# `-d ''` reads to the end of the value, so a newline inside it is one more separator;
	# a plain `read` stops at the first newline and drops what follows. `read` exits 1 at
	# the end of input without a NUL, which is the normal case here.
	read -r -d '' -a VEP_EXTRA_ARGV <<<"$VEP_EXTRA_ARGS" || true
	if [[ ${#VEP_EXTRA_ARGV[@]} -gt 0 ]]; then
		VEP_WRAPPER="$OUTPUT_DIR/vep-wrapper.sh"
		VEP_ARGV_LOG="$OUTPUT_DIR/vep-wrapper.argv"
		: >"$VEP_ARGV_LOG"
		real_binary_q=$(printf '%q' "$VEP_REAL_BINARY")
		extra_args_q=$(printf ' %q' "${VEP_EXTRA_ARGV[@]}")
		argv_log_q=$(printf '%q' "$VEP_ARGV_LOG")
		cat >"$VEP_WRAPPER" <<EOF
#!/usr/bin/env bash
# Written by run_clone_measurement.sh for --vep-extra-args: runs the measured binary
# with the extra arguments after the harness's own, and appends every argument vector
# it runs to the log beside it, one line per invocation, fields separated by \\037.
{ printf '%s\\037' $real_binary_q "\$@"$extra_args_q; echo; } >>$argv_log_q
exec $real_binary_q "\$@"$extra_args_q
EOF
		chmod +x "$VEP_WRAPPER"
		VEP_BINARY="$VEP_WRAPPER"
		echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] vep-rs extra arguments: ${VEP_EXTRA_ARGV[*]} (through $VEP_WRAPPER)"
	fi
	# `${VEP_EXTRA_ARGV[@]+...}`, not a bare `"${VEP_EXTRA_ARGV[@]}"`: on macOS bash 3.2 an
	# empty array is an unbound variable under `set -u`.
	python3 - "$OUTPUT_DIR/provenance.json" "$VEP_REAL_BINARY" "$VEP_EXTRA_ARGS" "$VEP_WRAPPER" ${VEP_EXTRA_ARGV[@]+"${VEP_EXTRA_ARGV[@]}"} <<'PY'
import hashlib, json, sys
out, binary, extra, wrapper, *argv = sys.argv[1:]
record = {"engine": "vep-rs", "vep_binary": binary,
          "vep_binary_sha256": hashlib.sha256(open(binary, "rb").read()).hexdigest(),
          "vep_extra_args": extra, "vep_extra_argv": argv,
          "wrapper": wrapper or None, "executed_argv": []}
json.dump(record, open(out, "w"), indent=1, sort_keys=True)
PY
fi

# --fastvep-extra-args reaches fastVEP the same way, through a wrapper that FASTVEP_BINARY is
# pointed at, so both fastVEP command lines below stay as the header documents them. A run
# without the words writes no wrapper and no record; with them, provenance.json names the
# binary, the words and the wrapper, and its fastvep_executed_argv is filled after the last
# cell from the wrapper's log.
FASTVEP_ARGV_LOG=""
if [[ ${#FASTVEP_EXTRA_ARGV[@]} -gt 0 ]]; then
	FASTVEP_WRAPPER="$OUTPUT_DIR/fastvep-wrapper.sh"
	FASTVEP_ARGV_LOG="$OUTPUT_DIR/fastvep-wrapper.argv"
	: >"$FASTVEP_ARGV_LOG"
	real_q=$(printf '%q' "$FASTVEP_BINARY")
	extra_q=$(printf ' %q' "${FASTVEP_EXTRA_ARGV[@]}")
	log_q=$(printf '%q' "$FASTVEP_ARGV_LOG")
	cat >"$FASTVEP_WRAPPER" <<EOF
#!/usr/bin/env bash
# Written by run_clone_measurement.sh for --fastvep-extra-args: runs the measured binary
# with the extra arguments after the harness's own, and appends every argument vector
# it runs to the log beside it, one line per invocation, fields separated by \\037.
{ printf '%s\\037' $real_q "\$@"$extra_q; echo; } >>$log_q
exec $real_q "\$@"$extra_q
EOF
	chmod +x "$FASTVEP_WRAPPER"
	python3 - "$OUTPUT_DIR/provenance.json" "$FASTVEP_BINARY" "$FASTVEP_WRAPPER" "${FASTVEP_EXTRA_ARGV[@]}" <<'PY'
import json, sys
out, binary, wrapper, *argv = sys.argv[1:]
record = {"engine": "fastvep", "fastvep_binary": binary,
          "fastvep_extra_args": " ".join(argv), "fastvep_extra_argv": argv,
          "fastvep_wrapper": wrapper, "fastvep_executed_argv": []}
json.dump(record, open(out, "w"), indent=1, sort_keys=True)
PY
	FASTVEP_BINARY="$FASTVEP_WRAPPER"
	echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] fastVEP extra arguments: ${FASTVEP_EXTRA_ARGV[*]} (through $FASTVEP_WRAPPER)"
fi

# discard_warmup_output <dataset_dir>
# The discarded warmup's output must be gone before the timed run starts: its pages,
# dirty or clean, are memory the timed run would otherwise reclaim while writing its
# own, and the reclaim lands in the timed wall time while user and system time stay
# unchanged. The warmup's purpose, a warm page cache for the inputs, survives the
# flush and the delete.
discard_warmup_output() {
	sync
	find "$1/rust" "$1/perl_scratch" -maxdepth 1 -type f \( -name '*.txt' -o -name '*.stderr' \) -delete 2>/dev/null || true
}

# time_run <label> <out_time_txt> -- <cmd...>
# Runs the command under GNU /usr/bin/time -v when available, capturing the
# timing block to <out_time_txt>. Echoes "<wall_seconds>|<exit_code>". Never
# aborts the run on a non-zero engine exit (recorded as exit_code in the CSV
# instead).
#
# GNU time is not present on macOS or the BSDs, where /usr/bin/time takes no -v.
# Without it the wall time still comes from the surrounding `date` delta, so the
# measurement holds; only the peak-RSS block in <out_time_txt> is lost. Install
# GNU time (`brew install gnu-time`) and set VEP_TIME_BIN=gtime to recover it.
time_run() {
	local label="$1" out_txt="$2"
	shift 2
	[[ "$1" == "--" ]] && shift
	local start end rc
	start=$(date +%s.%N)
	set +e
	# GNU time propagates the command's exit status, so the fallback to a bare run
	# is decided by whether GNU time is usable, never by the command failing.
	if "${VEP_TIME_BIN:-/usr/bin/time}" -v true >/dev/null 2>&1; then
		"${VEP_TIME_BIN:-/usr/bin/time}" -v "$@" 2>"$out_txt"
	else
		"$@" 2>"$out_txt"
	fi
	rc=$?
	set -e
	end=$(date +%s.%N)
	# Prefer GNU time's "Elapsed (wall clock)"; fall back to the date delta.
	local wall
	wall=$(awk '
		/Elapsed \(wall clock\)/ {
			n=split($NF, a, ":");
			if (n==3) { print a[1]*3600 + a[2]*60 + a[3] }
			else if (n==2) { print a[1]*60 + a[2] }
			else { print $NF }
			found=1; exit
		}
		END { if (!found) print "" }
	' "$out_txt")
	if [[ -z "$wall" ]]; then
		wall=$(awk "BEGIN{printf \"%.3f\", $end - $start}")
	fi
	echo "$wall|$rc"
}

# csv_row <dataset_id> <dataset_name> <assembly> <run_kind> <wall> <rc>
csv_row() {
	echo "$ENGINE,$RUN_ID,$HOST_ID,$REPLICATE,$1,$2,$3,$4,$5,$6" >>"$WALL_CSV"
}

# Ensembl VEP canonical-contig filter (Ensembl VEP engine only). The cached SNP/indel ground
# truth was built from canonical-filtered inputs, so the timed Ensembl VEP run filters the
# same way.
#
# ALSO drops per-sample genotypes (bcftools view -G -> sites-only). VEP annotates
# per (position, allele) and its tab output carries NO genotype columns, so this
# leaves the annotation output (and therefore concordance F1, which keys on
# location/allele/feature/consequence) byte-identical, and it is what makes the
# multi-sample datasets' timed cells comparable across engines: on s06 (1KG_highcov,
# 3,202 samples) the genotype matrix is most of the input bytes, and its cost is
# measured separately by the `*m` datasets. Of the other five SNP/indel datasets of a
# default run, s01 to s04 hold no genotype columns at source and s05 reads its sites-only
# copy, so -G is a harmless no-op for them. The structural-variant datasets do not pass
# through this filter.
#
# FAILS CLOSED at every step. There is no cp fallback: an unfiltered or
# genotype-carrying VCF written to the canonical path would be measured and
# reported as if the canonical filter had been applied, and on the multi-sample
# datasets would time the genotype-carrying input under a sites-only label. A missing bcftools,
# a failed filter, or a failed index each abort the cell instead.
CANONICAL_CONTIGS_LIST="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y,MT,M"
canonicalize() {
	local in_vcf="$1" out_vcf="$2"
	if ! command -v bcftools >/dev/null 2>&1; then
		# FAIL CLOSED. Every canonical-contigs guarantee in the reported results
		# rests on this filter, and the cached Ensembl VEP ground truth was
		# generated from canonical-filtered inputs, so annotating the raw input
		# instead measures a DIFFERENT input than the one the comparison assumes.
		# A bare cp with only a WARN here would silently void the policy claim and
		# (on the multi-sample datasets) time the genotype matrix, whose cost is
		# measured separately; genotypes change no consequence.
		echo "ERROR: [run_clone_measurement] bcftools unavailable; cannot produce a canonical-contig input for $in_vcf" >&2
		return 1
	fi
	[[ -f "${in_vcf}.tbi" ]] || tabix -p vcf "$in_vcf" 2>/dev/null || true
	local first_chrom src
	first_chrom=$(zcat -f "$in_vcf" | grep -v '^#' | head -1 | cut -f1 || echo "")
	src="$in_vcf"
	if [[ "$first_chrom" =~ ^chr ]]; then
		local stripped="${out_vcf%.vcf.gz}.stripped.vcf.gz"
		# The pipe's status is read under pipefail: a pass that stops part-way has written
		# the records before the stop, and bcftools would filter them as a whole input.
		if ! "${FILE_CENSUS[@]}" strip-chr "$in_vcf" | bgzip >"$stripped"; then
			echo "ERROR: [run_clone_measurement] contig-prefix pass failed for $in_vcf" >&2
			rm -f "$stripped"
			return 1
		fi
		tabix -p vcf "$stripped" 2>/dev/null || true
		src="$stripped"
	fi
	# -G drops genotypes (sites-only); -t applies the canonical-contig filter.
	#
	# Both steps FAIL CLOSED. A bare `cp` of the raw input on bcftools failure
	# would leave an unfiltered, genotype-carrying VCF at the canonical path --
	# measured and reported as if the filter had been applied. A canonical
	# filter that can silently not apply is not a filter.
	if ! bcftools view -G -t "$CANONICAL_CONTIGS_LIST" -O z -o "$out_vcf" "$src" 2>/dev/null; then
		echo "ERROR: [run_clone_measurement] canonical-contig + sites-only filter failed for $src" >&2
		rm -f "$out_vcf"
		return 1
	fi
	# The index is required: the record-count assertions and Ensembl VEP both read it.
	if ! tabix -p vcf "$out_vcf" 2>/dev/null; then
		echo "ERROR: [run_clone_measurement] tabix index failed for $out_vcf" >&2
		return 1
	fi
}

# True when this dataset has Ensembl VEP ground truth to compare against. Distinguishes an
# EXPECTED absence from a missing file, because conflating the two makes a run
# report failure while producing perfect data.
#
# The sample-scaling datasets (`s05m`/`s06m`) have no Ensembl VEP ground truth: the ground truth is
# generated from the sites-only inputs, and genotypes change no consequence, so the
# sites-only truth is the truth for both forms (see `should_run` above). Symlinking a nonexistent
# target anyway leaves a DANGLING link, which a subsequent copy of the output tree is
# liable to treat as an error, turning a complete run into a reported failure.
#
# A missing ground truth on any OTHER dataset is a data-tree error and fails the cell loudly.
# Keying on existence alone would let a dataset whose ground truth is absent read as a
# clean skip, which is the same conflation in the opposite direction.
ground_truth_ready() {
	local sid="$1" sname="$2" asm="$3" gt="$4"
	if [[ "$sid" == *m ]]; then
		echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)]   $sid is sample-scaling: wall time only, no Perl ground truth exists"
		return 1
	fi
	if [[ ! -f "$gt/output.txt" ]]; then
		echo "[ERROR] $sid: no Perl ground truth at $gt/output.txt" >&2
		csv_row "$sid" "$sname" "$asm" "ground_truth_missing" "" "1"
		return 1
	fi
	return 0
}

# reference_release_for <ground-truth dir>
#
# Prints the --reference-release for a cell: the --reference-release given to this run,
# else the release named by the ground truth's provenance (reference_release.py reads the
# directory's provenance.json, else the sidecar beside its output, else the set record in
# the parent directory; the image tag first, else the cache version). Fails when nothing
# names one, so a ground truth of unknown release is never scored under a guessed
# registry; the caller records the cell as compare_failed before timing it.
reference_release_for() {
	local gt_dir="$1"
	if [[ -n "$REFERENCE_RELEASE" ]]; then
		printf '%s\n' "$REFERENCE_RELEASE"
		return 0
	fi
	command -v python3 >/dev/null 2>&1 || {
		echo "ERROR: [run_clone_measurement] python3 is not on PATH: the reference release is derived by $SCRIPT_DIR/reference_release.py; install it or pass --reference-release" >&2
		return 1
	}
	python3 "$SCRIPT_DIR/reference_release.py" --provenance-dir "$gt_dir" || {
		echo "ERROR: [run_clone_measurement] no reference release for $gt_dir: its provenance names neither an Ensembl VEP image tag nor a cache version; pass --reference-release" >&2
		return 1
	}
}

# Per-engine SNP/indel measurement
# measure_p1 covers the SNP/indel datasets (dataset type `snp_indel`).
measure_p1() {
	local dataset_id="$1" dataset_name="$2" assembly="$3" vcf_base="$4" species="$5"
	local dataset_dir="$OUTPUT_DIR/$dataset_id"
	mkdir -p "$dataset_dir/rust" "$dataset_dir/report" "$dataset_dir/perf"
	local fasta
	fasta=$(fasta_for "$species" "$assembly") || {
		csv_row "$dataset_id" "$dataset_name" "$assembly" "genome_resources_missing" "" "1"
		return 1
	}
	# Canonical input (Ensembl-style contigs) -- the set the
	# ground_truth/perl/snp_indel/ Ensembl VEP GT was generated from.
	local vcf
	vcf="$DATA_DIR/vcf/$(printf %s "$assembly" | tr "[:upper:]" "[:lower:]")/all_variants_canonical/${vcf_base}.canonical.vcf.gz"
	local gt_dir="$SNP_INDEL_GT_DIR/$dataset_id"
	# The comparator's --reference-release, derived BEFORE the warmup and timed runs so a
	# ground truth that names no release fails the cell before any annotation is timed.
	# Only a cell that will compare derives it: the Ensembl VEP engine has no comparator, a
	# sample-scaling dataset (`*m`) has no ground truth, and a ground truth that is absent
	# altogether is reported by ground_truth_ready after the runs.
	local reference_release=""
	if [[ "$WALLTIME_ONLY" != "true" && "$ENGINE" != "perl" && "$dataset_id" != *m && -f "$gt_dir/output.txt" ]]; then
		reference_release=$(reference_release_for "$gt_dir") || {
			csv_row "$dataset_id" "$dataset_name" "$assembly" "compare_failed" "" "1"
			return 1
		}
	fi

	case "$ENGINE" in
	vep-rs)
		local cache
		cache=$(veprs_cache_for "$species" "$assembly") || {
			csv_row "$dataset_id" "$dataset_name" "$assembly" "genome_resources_missing" "" "1"
			return 1
		}
		local out="$dataset_dir/rust/output.txt"
		local run_kind perf_txt res wall rc
		for run_kind in warmup timed; do
			perf_txt="$dataset_dir/perf/${run_kind}.time.txt"
			res=$(time_run "${dataset_id}_${run_kind}" "$perf_txt" -- \
				"$VEP_BINARY" -i "$vcf" -o "$out" \
				--offline --json_cache "$cache" \
				--fasta "$fasta" \
				--fork "$FORK" --buffer_size "$BUFFER_SIZE" \
				--force --no_stats --quiet)
			wall="${res%%|*}"
			rc="${res##*|}"
			echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)]   $dataset_id $run_kind wall=${wall}s rc=$rc"
			csv_row "$dataset_id" "$dataset_name" "$assembly" "$run_kind" "$wall" "$rc"
			[[ "$run_kind" == warmup ]] && discard_warmup_output "$dataset_dir"
		done
		if [[ "$WALLTIME_ONLY" != "true" ]] && ground_truth_ready "$dataset_id" "$dataset_name" "$assembly" "$gt_dir"; then
			mkdir -p "$dataset_dir/perl"
			ln -sf "$gt_dir/output.txt" "$dataset_dir/perl/output.txt"
			"${OUTPUTS_COMPARE[@]}" \
				--perl-dir "$dataset_dir/perl" \
				--rust-dir "$dataset_dir/rust" \
				--report-dir "$dataset_dir/report" \
				--reference-release "$reference_release" \
				--input-vcf "$vcf" \
				--assembly "$assembly" || {
				echo "ERROR: [run_clone_measurement] $dataset_id SNP/indel compare failed; the cell has no F1 and this run is incomplete" >&2
				csv_row "$dataset_id" "$dataset_name" "$assembly" "compare_failed" "" "1"
				return 1
			}
			# Per-column and Extra-key agreement, beside the F1 report. Not F1-bearing:
			# a failure is logged in the CSV and the cell still completes.
			if [[ "$SKIP_FIELDS" != "true" ]]; then
				"${FIELDS_COMPARE[@]}" \
					--perl "$dataset_dir/perl/output.txt" \
					--rust "$dataset_dir/rust/output.txt" \
					--out "$dataset_dir/report/fields.json" \
					--tmp-dir "$dataset_dir" >/dev/null || {
					echo "WARNING: [run_clone_measurement] $dataset_id field comparison failed; report/fields.json is absent" >&2
					csv_row "$dataset_id" "$dataset_name" "$assembly" "fields_compare_failed" "" "1"
				}
			else
				echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] $dataset_id field comparison skipped (--skip-fields)"
			fi
		fi
		;;
	fastvep)
		local cache gff3
		cache=$(fastvep_cache_for "$species" "$assembly") || {
			csv_row "$dataset_id" "$dataset_name" "$assembly" "genome_resources_missing" "" "1"
			return 1
		}
		# fastVEP, given a --transcript-cache path that names no file, parses the --gff3 file,
		# saves a cache of its own at the path and exits 0: the cell would be timed and scored
		# on the GFF3's gene models with nothing recorded.
		[[ -f "$cache" && -s "$cache" ]] || {
			echo "ERROR: [run_clone_measurement] $dataset_id: no fastVEP transcript cache file with bytes at $cache; where the path names no file fastVEP parses the GFF3 and saves a cache of its own there, so the cell would be timed and scored on the GFF3's gene models" >&2
			csv_row "$dataset_id" "$dataset_name" "$assembly" "fastvep_cache_missing" "" "1"
			return 1
		}
		gff3=$(gff3_for "$species" "$assembly") || {
			csv_row "$dataset_id" "$dataset_name" "$assembly" "genome_resources_missing" "" "1"
			return 1
		}
		local raw="$dataset_dir/rust/raw.txt"
		local run_kind perf_txt res wall rc
		for run_kind in warmup timed; do
			perf_txt="$dataset_dir/perf/${run_kind}.time.txt"
			res=$(time_run "${dataset_id}_${run_kind}" "$perf_txt" -- \
				env RAYON_NUM_THREADS="$FORK" "$FASTVEP_BINARY" annotate \
				--gff3 "$gff3" --fasta "$fasta" --transcript-cache "$cache" \
				-i "$vcf" -o "$raw" --output-format tab)
			wall="${res%%|*}"
			rc="${res##*|}"
			echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)]   $dataset_id $run_kind wall=${wall}s rc=$rc"
			csv_row "$dataset_id" "$dataset_name" "$assembly" "$run_kind" "$wall" "$rc"
			[[ "$run_kind" == warmup ]] && discard_warmup_output "$dataset_dir"
		done
		if [[ "$WALLTIME_ONLY" != "true" ]] && ground_truth_ready "$dataset_id" "$dataset_name" "$assembly" "$gt_dir"; then
			# Normalize the timed raw output (intergenic Feature_type=Transcript -> -).
			# FAIL CLOSED. A `cp` fallback here would substitute un-normalised output
			# and depress fastVEP's F1 for a reason unrelated to annotation quality,
			# with no error surfaced. A missing or broken normaliser must abort the
			# cell instead.
			if [[ -s "$raw" ]]; then
				# `return 1`, never `exit 1`: this runs inside measure_p1, whose call
				# site is guarded (`measure_p1 ... || failed_datasets+=(...)`). Exiting
				# here would abort the whole run, skipping every later dataset and any
				# cleanup the caller performs.
				"${FASTVEP_NORMALIZE[@]}" "$raw" "$dataset_dir/rust/output.txt" || {
					echo "[ERROR] fastvep normalizer failed on $dataset_id; refusing to compare un-normalised output" >&2
					csv_row "$dataset_id" "$dataset_name" "$assembly" "normalizer_failed" "" "1"
					return 1
				}
			fi
			mkdir -p "$dataset_dir/perl"
			ln -sf "$gt_dir/output.txt" "$dataset_dir/perl/output.txt"
			"${OUTPUTS_COMPARE[@]}" \
				--perl-dir "$dataset_dir/perl" \
				--rust-dir "$dataset_dir/rust" \
				--report-dir "$dataset_dir/report" \
				--reference-release "$reference_release" \
				--input-vcf "$vcf" \
				--assembly "$assembly" --glob 'output.txt' || {
				echo "ERROR: [run_clone_measurement] $dataset_id SNP/indel compare failed; the cell has no F1 and this run is incomplete" >&2
				csv_row "$dataset_id" "$dataset_name" "$assembly" "compare_failed" "" "1"
				return 1
			}
		fi
		;;
	perl)
		# Ensembl VEP is the ground-truth generator: it has no concordance branch. The
		# scratch annotation is timed then deleted (the cached GT is the source
		# of truth; this engine contributes wall-time only).
		#
		# Cache completeness is asserted BEFORE any work: a chromosome-subset cache
		# produces a fast, wrong wall time and exits 0.
		if ! assert_perl_cache_complete "$species" "$assembly" "$(cache_version_for "$species" "$assembly")"; then
			csv_row "$dataset_id" "$dataset_name" "$assembly" "perl_cache_incomplete" "" "1"
			return 1
		fi
		local canon_vcf="$dataset_dir/inputs/${vcf_base}.canonical.vcf.gz"
		local raw_vcf
		raw_vcf="$DATA_DIR/vcf/$(printf %s "$assembly" | tr "[:upper:]" "[:lower:]")/all_variants/${vcf_base}.vcf.gz"
		mkdir -p "$dataset_dir/inputs" "$dataset_dir/perl_scratch"
		# canonicalize fails closed (see its definition). A cell whose input is not
		# canonical cannot be compared against the canonical-filtered ground truth,
		# so record the failure and skip rather than time a different input.
		if ! canonicalize "$raw_vcf" "$canon_vcf"; then
			echo "ERROR: [run_clone_measurement] $dataset_id skipped: canonical input unavailable" >&2
			csv_row "$dataset_id" "$dataset_name" "$assembly" "canonicalize_failed" "" "1"
			return 1
		fi
		local run_kind perf_txt res wall rc
		for run_kind in warmup timed; do
			perf_txt="$dataset_dir/perf/${run_kind}.time.txt"
			res=$(perl_annotate_timed "${dataset_id}_${run_kind}" "$perf_txt" \
				"$canon_vcf" "$dataset_dir/perl_scratch/output.txt" "$species" "$assembly")
			wall="${res%%|*}"
			rc="${res##*|}"
			echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)]   $dataset_id $run_kind wall=${wall}s rc=$rc"
			csv_row "$dataset_id" "$dataset_name" "$assembly" "$run_kind" "$wall" "$rc"
			[[ "$run_kind" == warmup ]] && discard_warmup_output "$dataset_dir"
		done
		# Before the scratch directory goes, confirm the timed run actually annotated. This is
		# the only point at which the evidence exists.
		if ! assert_perl_annotated "$dataset_dir/perl_scratch/output.txt" "$canon_vcf" "$dataset_id"; then
			csv_row "$dataset_id" "$dataset_name" "$assembly" "perl_not_annotating" "" "1"
			rm -rf "$dataset_dir/perl_scratch" "$dataset_dir/inputs"
			return 1
		fi
		rm -rf "$dataset_dir/perl_scratch" "$dataset_dir/inputs"
		;;
	esac
}

# perl_annotate_timed <label> <perf_txt> <in.vcf.gz> <out.txt> <species> <assembly>
# Runs Ensembl VEP in Docker under /usr/bin/time. --fork matches the value the other
# engines are given, so the comparison is like for like; the SV ground truth's own
# fork setting belongs to generate_reference_output.sh, not to a timed cell.
perl_annotate_timed() {
	local label="$1" perf_txt="$2" in_vcf="$3" out_txt="$4" species="$5" assembly="$6"
	local fasta out_dir out_base cache_version
	fasta=$(fasta_for "$species" "$assembly")
	cache_version=$(cache_version_for "$species" "$assembly")
	out_dir=$(dirname "$out_txt")
	out_base=$(basename "$out_txt")
	mkdir -p "$out_dir"
	time_run "$label" "$perf_txt" -- \
		docker run --rm \
		--user "$(id -u):$(id -g)" \
		-v "${in_vcf}:/work/input.vcf.gz:ro" \
		-v "${out_dir}:/work/perl_out" \
		-v "${PERL_CACHE_DIR}:/work/cache:ro" \
		-v "${fasta}:/work/ref.fa:ro" \
		-v "${fasta}.fai:/work/ref.fa.fai:ro" \
		"$PERL_IMAGE" \
		vep \
		-i /work/input.vcf.gz \
		-o "/work/perl_out/${out_base}" \
		--offline --cache --dir_cache /work/cache --cache_version "$cache_version" \
		--species "$species" --assembly "$assembly" \
		--format vcf --buffer_size "$BUFFER_SIZE" --fork "$FORK" \
		--fasta /work/ref.fa \
		--no_headers --no_stats --quiet --force_overwrite
}

# assert_perl_annotated <perl_out_file_or_dir> <input_file_or_dir> <dataset_id>
#
# An Ensembl VEP invocation that cannot reach its cache does not fail: it emits one warning per
# record, zero annotation rows, and exit 0. Because the timed run's output is deleted immediately
# after the loop, nothing downstream ever sees it, so without this gate an Ensembl VEP wall time
# can measure a VEP that annotated nothing.
#
# THE FLOOR IS REFERENCE-FREE ON PURPOSE. It needs no ground truth staged, so it cannot be
# skipped on a data tree that carries only inputs. A working VEP emits at least one row per input
# record (one per overlapping transcript, or a single intergenic row) and several per record on
# genome-wide and structural-variant inputs alike. So `rows >= records` fails a broken run, which
# sits at zero, and passes a legitimately sparse one with an order of magnitude spare.
#
# It is not specific to one flag: a wrong cache mount, a wrong --cache_version, a wrong
# --assembly and a missing cache all present identically as "almost no output".
#
# `n=$(...) || n=0`, never `$(... || echo 0)`: the count prints 0 AND exits 1 when no line
# counts, as grep -c does,
# so the latter form makes n the two-line string "0\n0", and every later numeric test on it dies
# with "integer expression expected" and falls through to its else branch -- reporting a pass
# because the comparison ERRORED rather than because it succeeded.
assert_perl_annotated() {
	local perl_out="$1" inputs="$2" dataset_id="$3"
	local out_rows=0 in_records=0 f n

	if [[ -d "$perl_out" ]]; then
		for f in "$perl_out"/*.txt; do
			[[ -e "$f" ]] || continue
			# VEP writes <out>_warnings.txt beside its output; counting it would mask the
			# defect exactly, since the broken run puts one line per record in there.
			case "$f" in *_warnings.txt) continue ;; esac
			n=$("${FILE_CENSUS[@]}" count --not-starting '#' "$f") || n=0
			out_rows=$((out_rows + n))
		done
	elif [[ -e "$perl_out" ]]; then
		out_rows=$("${FILE_CENSUS[@]}" count --not-starting '#' "$perl_out") || out_rows=0
	fi

	for f in "$inputs" "$inputs"/*.vcf "$inputs"/*.vcf.gz; do
		[[ -f "$f" ]] || continue
		case "$f" in
		*.vcf.gz) n=$("${FILE_CENSUS[@]}" count --gunzip --not-starting '#' "$f") || n=0 ;;
		*.vcf) n=$("${FILE_CENSUS[@]}" count --not-starting '#' "$f") || n=0 ;;
		*) continue ;;
		esac
		in_records=$((in_records + n))
	done

	echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)]   $dataset_id perl annotation check: ${out_rows} rows / ${in_records} records"
	# ZERO INPUT RECORDS IS A FAILURE, NOT A PASS. The floor below is `out_rows >=
	# in_records`, so `in_records == 0` would make the whole gate vacuous: an absent or empty
	# input directory yields 0 rows against 0 records, the predicate is false, and the cell
	# is recorded as a clean, very fast measurement. `shopt -s nullglob` on the SV branch
	# makes an empty directory iterate zero times and exit 0, so nothing upstream objects
	# either, and the cache guard passes because the cache is fine. That is a wrong
	# number reached with both gates green.
	if [[ "$in_records" -eq 0 ]]; then
		echo "ERROR: [run_clone_measurement] $dataset_id counted 0 input records" >&2
		echo "  The inputs are missing or empty, so the output-volume floor cannot bound anything." >&2
		echo "  A cell timed in this state records a near-zero wall time as a success." >&2
		return 1
	fi
	if [[ "$out_rows" -lt "$in_records" ]]; then
		echo "ERROR: [run_clone_measurement] $dataset_id Perl emitted $out_rows rows for $in_records input records" >&2
		echo "  An Ensembl VEP run that cannot reach its cache emits warnings rather than annotations and still exits 0." >&2
		echo "  Check --dir_cache, the cache mount, --cache_version and --assembly before trusting this cell." >&2
		return 1
	fi
	return 0
}

# assert_perl_cache_complete <species> <assembly> [<cache version>]
#
# A chromosome-subset cache (a `116_GRCh37/` holding only `21/` against a genome-wide
# `116_GRCh38/`) is a second way for Ensembl VEP to annotate almost nothing and exit 0, independent
# of the flag.
#
# The reference-free `rows >= records` floor of `assert_perl_annotated` catches a chr21-only
# cache only where the input is genome-wide. Over the 16 GRCh37 SV VCFs, whose records are
# chr21 by design, such a cache passes the floor with an order of magnitude spare while
# `09_breakends` loses every mate that points off chromosome 21. An output-volume gate
# cannot see a cache gap the input barely samples, so completeness is asserted on the cache
# itself, before any timing.
#
# Primary contigs only. Ensembl's full GRCh37 cache carries 699 directories, the other 674 being
# unplaced scaffolds and haplotype contigs, and no input here reaches them -- so requiring all
# 699 would fail a legitimately pruned cache. Requiring the 25 primary contigs fails every
# chromosome-subset cache, which is the defect class. The 25 are the human primary set and
# are required under homo_sapiens; a slice of another species is required to hold at least
# one chromosome directory. The slice is read under the cache version given, the caller's
# cache_version_for, and under PERL_CACHE_VERSION without one.
assert_perl_cache_complete() {
	local species="$1" assembly="$2"
	# The default is load-bearing rather than defensive. The CACHE-COMPLETE gate in
	# test_run_clone_measurement.sh extracts this function with awk and evals it, so the gate
	# runs the shipped code and not a copy that can drift from it. Under that eval the harness's
	# own `PERL_CACHE_VERSION` assignment has not run and the test file is `set -u`, so a bare
	# reference would abort the whole test file on an unbound variable, exiting 1 with no
	# failing assertion printed and reading as a harness defect rather than a missing default.
	local ver="${3:-${PERL_CACHE_VERSION:-116}}"
	local slice="$PERL_CACHE_DIR/$species/${ver}_${assembly}"
	local primaries="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT"
	local c missing="" n_chrom=0

	if [[ ! -d "$slice" ]]; then
		echo "ERROR: [run_clone_measurement] Perl cache slice not found: $slice" >&2
		echo "  Expected $species/${ver}_${assembly} under \$PERL_CACHE_DIR ($PERL_CACHE_DIR)." >&2
		return 1
	fi
	if [[ "$species" != "homo_sapiens" ]]; then
		for c in "$slice"/*/; do
			[[ -d "$c" ]] && n_chrom=$((n_chrom + 1))
		done
		if [[ "$n_chrom" -eq 0 ]]; then
			echo "ERROR: [run_clone_measurement] Perl cache slice holds no chromosome directory: $slice" >&2
			return 1
		fi
		echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)]   perl cache check ($species $assembly): $n_chrom chromosome directories; a slice of a species other than homo_sapiens is checked to hold at least one"
		return 0
	fi

	for c in $primaries; do
		[[ -d "$slice/$c" ]] || missing="$missing $c"
	done

	if [[ -n "$missing" ]]; then
		echo "ERROR: [run_clone_measurement] Perl cache for $assembly is missing primary contigs:$missing" >&2
		echo "  Slice: $slice" >&2
		echo "  A chromosome-subset cache does not fail a VEP run. It under-annotates every record" >&2
		echo "  outside the contigs it holds, exits 0, and yields a wall time that is too fast." >&2
		echo "  Do NOT time this cell against it." >&2
		return 1
	fi

	echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)]   perl cache check ($assembly): all 25 primary contigs present"
	return 0
}

# Per-engine SV measurement (per-VCF over the canonical 16-VCF set per assembly)
# assert_sv_inputs_present <inputs_dir> <dataset_id>
# assert_sv_outputs_present <outputs_dir> <n_inputs> <dataset_id>
#
# The per-VCF loop runs under `nullglob`, so an input directory that is absent, empty or
# holds only manifest.json is not an error to bash: the loop body never runs, `time_run`
# reports a 0 s wall time with rc=0, the cell writes no output, and the comparison that
# follows finds nothing to score. A cell that annotated nothing must FAIL, so both ends
# of the loop are asserted: at least one input VCF before it runs, and one non-empty
# output per input after it. Reference-free on purpose: neither needs the ground truth
# staged.
sv_input_count() {
	local dir="$1" n=0 f b
	shopt -s nullglob
	for f in "$dir"/*.vcf "$dir"/*.vcf.gz; do
		b=$(basename "$f"); b=${b%.vcf.gz}; b=${b%.vcf}
		[[ "$b" == "manifest" ]] && continue
		n=$((n + 1))
	done
	shopt -u nullglob
	echo "$n"
}
assert_sv_inputs_present() {
	local dir="$1" sid="$2" n
	n=$(sv_input_count "$dir")
	if [[ "$n" -lt 1 ]]; then
		echo "ERROR: [run_clone_measurement] $sid: no SV input VCF under $dir (SV_GT_DIR=$SV_GT_DIR); the cell would time an empty loop" >&2
		return 1
	fi
	return 0
}
assert_sv_outputs_present() {
	local dir="$1" want="$2" sid="$3" n=0 f
	shopt -s nullglob
	for f in "$dir"/*.txt; do
		if grep -q -v '^#' "$f" 2>/dev/null; then n=$((n + 1)); fi
	done
	shopt -u nullglob
	if [[ "$n" -lt "$want" ]]; then
		echo "ERROR: [run_clone_measurement] $sid: $n non-empty SV output file(s) under $dir for $want input VCF(s); the engine annotated nothing on at least one file" >&2
		return 1
	fi
	return 0
}

measure_sv() {
	local dataset_id="$1" dataset_name="$2" assembly="$3" species="$4"
	local dataset_dir="$OUTPUT_DIR/$dataset_id"
	mkdir -p "$dataset_dir/rust" "$dataset_dir/report" "$dataset_dir/perf"
	# `${var,,}` is bash 4+ and fails at runtime on macOS bash 3.2 while passing
	# `bash -n` and shellcheck.
	local asm_lc
	asm_lc=$(printf %s "$assembly" | tr "[:upper:]" "[:lower:]")
	local sv_inputs="$SV_GT_DIR/canonical_inputs/$asm_lc"
	local sv_gt="$SV_GT_DIR/$asm_lc"
	if ! assert_sv_inputs_present "$sv_inputs" "$dataset_id"; then
		csv_row "$dataset_id" "$dataset_name" "$assembly" "sv_inputs_missing" "" "1"
		return 1
	fi
	local sv_n_inputs
	sv_n_inputs=$(sv_input_count "$sv_inputs")

	# Per-assembly dedupe: annotate once under s07/s08; s09/s10 are pointers.
	if [[ "$SV_DONE" == *":${asm_lc}="* ]]; then
		sv_done_dataset="${SV_DONE#*:${asm_lc}=}"
		sv_done_dataset="${sv_done_dataset%%:*}"
		echo "SV $assembly annotated under ${sv_done_dataset}; see that dataset's report (per-file F1 covers synthetic + real-world)" \
			>"$dataset_dir/POINTER"
		csv_row "$dataset_id" "$dataset_name" "$assembly" "pointer" "0" "0"
		return 0
	fi
	SV_DONE="${SV_DONE}:${asm_lc}=${dataset_id}:"

	# As in measure_p1: the ground truth's release is derived before the timed runs.
	local reference_release=""
	if [[ "$WALLTIME_ONLY" != "true" && "$ENGINE" != "perl" ]]; then
		reference_release=$(reference_release_for "$sv_gt") || {
			csv_row "$dataset_id" "$dataset_name" "$assembly" "compare_failed" "" "1"
			return 1
		}
	fi

	local fasta
	fasta=$(fasta_for "$species" "$assembly") || {
		csv_row "$dataset_id" "$dataset_name" "$assembly" "genome_resources_missing" "" "1"
		return 1
	}

	case "$ENGINE" in
	vep-rs)
		local cache
		cache=$(veprs_cache_for "$species" "$assembly") || {
			csv_row "$dataset_id" "$dataset_name" "$assembly" "genome_resources_missing" "" "1"
			return 1
		}
		local run_kind rust_out perf_txt res wall rc
		for run_kind in warmup timed; do
			rust_out="$dataset_dir/rust"
			perf_txt="$dataset_dir/perf/${run_kind}.time.txt"
			# Wrap the whole per-VCF loop in one /usr/bin/time so wall_seconds is
			# the full SV-cell annotation time.
			#
			# `$FORK` here, NOT a hardcoded 1. Every engine takes the same --fork on
			# every cell: the Ensembl VEP and fastVEP branches below both pass
			# `$FORK` on the SV cells too, so hardcoding 1 here would make vep-rs the
			# only engine of the three annotating serially and turn every SV ratio
			# into a floor rather than a measurement. The GT is a fixed archived
			# artifact, and vep-rs output is fork-invariant by construction:
			# `annotate_batch` dispatches the same `annotate_one` over the same `&mut`
			# batch through either `par_iter_mut` or `iter_mut` and emits in input
			# order (`runner.rs`), and the per-chromosome cache parse runs on the
			# global rayon pool regardless (`json_cache.rs::load_transcripts_for_chr`).
			# What --fork gates is only the per-variant consequence fan-out.
			res=$(time_run "${dataset_id}_${run_kind}" "$perf_txt" -- \
				bash -c '
					set -uo pipefail
					bin="$1"; cache="$2"; fasta="$3"; asm="$4"; buf="$5"; in_dir="$6"; out_dir="$7"; fork="$8"; species="$9"
					shopt -s nullglob
					failed=0
					for vcf in "$in_dir"/*.vcf "$in_dir"/*.vcf.gz; do
						b=$(basename "$vcf"); b=${b%.vcf.gz}; b=${b%.vcf}
						[[ "$b" == "manifest" ]] && continue
						"$bin" -i "$vcf" -o "$out_dir/$b.txt" \
							--offline --json_cache "$cache" --fasta "$fasta" \
							--fork "$fork" --buffer_size "$buf" --force --no_stats \
							--species "$species" --assembly "$asm" --quiet \
							2>"$out_dir/$b.stderr" || {
							failed=$((failed+1))
							echo "ERROR: [run_clone_measurement] per-VCF invocation failed: $b" >&2
						}
					done
					[[ "$failed" -eq 0 ]] || {
						echo "ERROR: [run_clone_measurement] $failed per-VCF invocation(s) failed" >&2
						exit "$failed"
					}
				' _ "$VEP_BINARY" "$cache" "$fasta" "$assembly" "$BUFFER_SIZE" "$sv_inputs" "$rust_out" "$FORK" "$species")
			wall="${res%%|*}"
			rc="${res##*|}"
			echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)]   $dataset_id $run_kind wall=${wall}s rc=$rc"
			# One non-empty output per input, or the cell is a failure however fast it
			# ran: a 0 s SV cell with rc=0 is the signature of an empty input glob.
			if ! assert_sv_outputs_present "$rust_out" "$sv_n_inputs" "$dataset_id"; then
				csv_row "$dataset_id" "$dataset_name" "$assembly" "${run_kind}_no_output" "$wall" "1"
				return 1
			fi
			csv_row "$dataset_id" "$dataset_name" "$assembly" "$run_kind" "$wall" "$rc"
			[[ "$run_kind" == warmup ]] && discard_warmup_output "$dataset_dir"
		done
		if [[ "$WALLTIME_ONLY" != "true" ]]; then
			# `--vep-rs-cache` is what activates the cross-chromosome mask, which removes
			# Ensembl VEP tuples naming a transcript the variant's chromosome does not hold. The
			# cache is the same one vep-rs just annotated against, so it is the right
			# authority. Passing it is NOT optional for a reported run: without it the
			# mask is inactive and the report says so.
			"$SV_COMPARE" \
				--assembly "$asm_lc" \
				--input-dir "$sv_inputs" \
				--perl-dir "$sv_gt" \
				--rust-dir "$dataset_dir/rust" \
				--vep-rs-cache "$cache" \
				--reference-release "$reference_release" \
				--output-dir "$dataset_dir/report" || {
				echo "ERROR: [run_clone_measurement] $dataset_id SV compare failed; the cell has no F1 and this run is incomplete" >&2
				csv_row "$dataset_id" "$dataset_name" "$assembly" "compare_failed" "" "1"
				return 1
			}
		fi
		;;
	fastvep)
		local cache gff3
		cache=$(fastvep_cache_for "$species" "$assembly") || {
			csv_row "$dataset_id" "$dataset_name" "$assembly" "genome_resources_missing" "" "1"
			return 1
		}
		# As in measure_p1: at a cache path that names no file fastVEP would parse the GFF3.
		[[ -f "$cache" && -s "$cache" ]] || {
			echo "ERROR: [run_clone_measurement] $dataset_id: no fastVEP transcript cache file with bytes at $cache; where the path names no file fastVEP parses the GFF3 and saves a cache of its own there, so the cell would be timed and scored on the GFF3's gene models" >&2
			csv_row "$dataset_id" "$dataset_name" "$assembly" "fastvep_cache_missing" "" "1"
			return 1
		}
		gff3=$(gff3_for "$species" "$assembly") || {
			csv_row "$dataset_id" "$dataset_name" "$assembly" "genome_resources_missing" "" "1"
			return 1
		}
		local run_kind rust_out perf_txt res wall rc
		for run_kind in warmup timed; do
			rust_out="$dataset_dir/rust"
			perf_txt="$dataset_dir/perf/${run_kind}.time.txt"
			res=$(time_run "${dataset_id}_${run_kind}" "$perf_txt" -- \
				bash -c '
					set -uo pipefail
					bin="$1"; cache="$2"; fasta="$3"; gff3="$4"; in_dir="$5"; out_dir="$6"; threads="$7"
					shopt -s nullglob
					failed=0
					for vcf in "$in_dir"/*.vcf "$in_dir"/*.vcf.gz; do
						b=$(basename "$vcf"); b=${b%.vcf.gz}; b=${b%.vcf}
						[[ "$b" == "manifest" ]] && continue
						RAYON_NUM_THREADS="$threads" "$bin" annotate \
							--gff3 "$gff3" --fasta "$fasta" --transcript-cache "$cache" \
							-i "$vcf" -o "$out_dir/$b.txt" --output-format tab \
							2>"$out_dir/$b.stderr" || {
							failed=$((failed+1))
							echo "ERROR: [run_clone_measurement] per-VCF invocation failed: $b" >&2
						}
					done
					[[ "$failed" -eq 0 ]] || {
						echo "ERROR: [run_clone_measurement] $failed per-VCF invocation(s) failed" >&2
						exit "$failed"
					}
				' _ "$FASTVEP_BINARY" "$cache" "$fasta" "$gff3" "$sv_inputs" "$rust_out" "$FORK")
			wall="${res%%|*}"
			rc="${res##*|}"
			echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)]   $dataset_id $run_kind wall=${wall}s rc=$rc"
			# Same gate as the vep-rs arm, at a floor of ONE non-empty output rather than
			# one per input: fastVEP's SV CONCORDANCE is best-effort below (it may decline
			# individual synthetic files), but a timing cell whose whole loop produced
			# nothing is the empty-glob signature, not a measurement.
			if ! assert_sv_outputs_present "$rust_out" 1 "$dataset_id"; then
				csv_row "$dataset_id" "$dataset_name" "$assembly" "${run_kind}_no_output" "$wall" "1"
				return 1
			fi
			csv_row "$dataset_id" "$dataset_name" "$assembly" "$run_kind" "$wall" "$rc"
			[[ "$run_kind" == warmup ]] && discard_warmup_output "$dataset_dir"
		done
		# SV concordance is best-effort for fastVEP: a loop that produces no non-empty
		# output is skipped, never hard-failed.
		if [[ "$WALLTIME_ONLY" != "true" ]]; then
			local produced
			produced=$(find "$dataset_dir/rust" -name '*.txt' -size +0c 2>/dev/null | wc -l | tr -d ' ')
			if [[ "$produced" -eq 0 ]]; then
				echo "fastVEP produced no non-empty SV output; SV concordance skipped" \
					>"$dataset_dir/report/SV_CONCORDANCE_SKIPPED.txt"
				echo "[WARN] $dataset_id ($assembly) fastVEP SV concordance skipped: no usable output"
			else
				# `--vep-rs-cache` on the fastVEP arm too, so this arm's adjusted F1 applies
				# the SAME exclusion set as the vep-rs arm's. Without it the cross-chromosome
				# class is inactive here and the two engines' adjusted SV columns sit on
				# different masks.
				# `--scored-engine fastvep` is required with it: the filter's sanity assertion
				# that the scored engine cannot name a transcript the cache lacks for the
				# variant's chromosome holds for vep-rs and not for an engine carrying its own
				# annotation release, and it aborts the whole comparison when it fires.
				#
				# `veprs_cache_for`, NOT the local `$cache`, which in this branch is fastVEP's
				# OWN cache. The flag wants the vep-rs JSON cache: it is read as the authority
				# on which chromosome a transcript sits on, in the layout that loader parses,
				# and it is not the transcript set the scored engine annotated against.
				# The mask authority is the vep-rs JSON cache, which fastVEP itself never
				# reads, so a fastVEP data tree can be complete for annotation and still
				# lack it; assert it by name before the comparator does with a traceback.
				if [[ ! -d "$(veprs_cache_for "$species" "$assembly")" ]]; then
					echo "ERROR: [run_clone_measurement] $dataset_id: vep-rs JSON cache absent at $(veprs_cache_for "$species" "$assembly"); the fastVEP tree must carry it for the cross-chromosome mask (--vep-rs-cache)" >&2
					csv_row "$dataset_id" "$dataset_name" "$assembly" "mask_cache_missing" "" "1"
					return 1
				fi
				"$SV_COMPARE" \
					--assembly "$asm_lc" \
					--input-dir "$sv_inputs" \
					--perl-dir "$sv_gt" \
					--rust-dir "$dataset_dir/rust" \
					--vep-rs-cache "$(veprs_cache_for "$species" "$assembly")" \
					--scored-engine fastvep \
					--reference-release "$reference_release" \
					--output-dir "$dataset_dir/report" || {
					# A compare that RAN and failed is not the documented no-output case
					# handled above: it means a missing input (the vep-rs cache the mask
					# reads, a ground-truth file) or a comparator error, and the cell has
					# no F1.
					echo "ERROR: [run_clone_measurement] $dataset_id fastVEP SV compare failed; the cell has no F1 and this run is incomplete" >&2
					csv_row "$dataset_id" "$dataset_name" "$assembly" "compare_failed" "" "1"
					return 1
				}
			fi
		fi
		;;
	perl)
		# Ensembl VEP SV: time over the canonical per-VCF SV inputs (already canonical at
		# GT-gen time). Scratch output deleted (Ensembl VEP is GT; wall-time only).
		#
		# This is the branch `assert_perl_annotated` cannot protect: these inputs are
		# chr21 by design, so a chr21-only cache clears its rows/records floor while
		# `09_breakends` loses every off-chromosome mate. See the guard's own comment.
		local cache_version
		cache_version=$(cache_version_for "$species" "$assembly")
		if ! assert_perl_cache_complete "$species" "$assembly" "$cache_version"; then
			csv_row "$dataset_id" "$dataset_name" "$assembly" "perl_cache_incomplete" "" "1"
			return 1
		fi
		local run_kind perf_txt res wall rc
		mkdir -p "$dataset_dir/perl_scratch"
		for run_kind in warmup timed; do
			perf_txt="$dataset_dir/perf/${run_kind}.time.txt"
			res=$(time_run "${dataset_id}_${run_kind}" "$perf_txt" -- \
				bash -c '
					set -uo pipefail
					image="$1"; cache="$2"; fasta="$3"; asm="$4"; in_dir="$5"; scratch="$6"; uidgid="$7"; buf="$8"; fork="$9"; cache_version="${10}"; species="${11}"
					mkdir -p "$scratch"
					shopt -s nullglob
					failed=0
					for vcf in "$in_dir"/*.vcf "$in_dir"/*.vcf.gz; do
						b=$(basename "$vcf"); b=${b%.vcf.gz}; b=${b%.vcf}
						[[ "$b" == "manifest" ]] && continue
						if [[ "$vcf" == *.vcf ]]; then
							bgzip -c "$vcf" >"$scratch/$b.vcf.gz" 2>/dev/null || cp "$vcf" "$scratch/$b.vcf.gz"
							invcf="$scratch/$b.vcf.gz"
						else
							invcf="$vcf"
						fi
						docker run --rm --user "$uidgid" \
							-v "$invcf:/work/input.vcf.gz:ro" \
							-v "$scratch:/work/perl_out" \
							-v "$cache:/work/cache:ro" \
							-v "$fasta:/work/ref.fa:ro" \
							-v "$fasta.fai:/work/ref.fa.fai:ro" \
							"$image" \
							vep -i /work/input.vcf.gz -o "/work/perl_out/$b.txt" \
							--offline --cache --dir_cache /work/cache --cache_version "$cache_version" \
							--species "$species" --assembly "$asm" \
							--format vcf --buffer_size "$buf" --fork "$fork" \
							--fasta /work/ref.fa \
							--no_headers --no_stats --quiet --force_overwrite || {
							failed=$((failed+1))
							echo "ERROR: [run_clone_measurement] per-VCF invocation failed: $b" >&2
						}
					done
					[[ "$failed" -eq 0 ]] || {
						echo "ERROR: [run_clone_measurement] $failed per-VCF invocation(s) failed" >&2
						exit "$failed"
					}
				' _ "$PERL_IMAGE" "$PERL_CACHE_DIR" "$fasta" "$assembly" \
				"$sv_inputs" "$dataset_dir/perl_scratch" "$(id -u):$(id -g)" "$BUFFER_SIZE" "$FORK" "$cache_version" "$species")
			wall="${res%%|*}"
			rc="${res##*|}"
			echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)]   $dataset_id $run_kind wall=${wall}s rc=$rc"
			csv_row "$dataset_id" "$dataset_name" "$assembly" "$run_kind" "$wall" "$rc"
			[[ "$run_kind" == warmup ]] && discard_warmup_output "$dataset_dir"
		done
		# Same gate as the SNP/indel path: the per-VCF Ensembl VEP output is about to be discarded.
		if ! assert_perl_annotated "$dataset_dir/perl_scratch" "$sv_inputs" "$dataset_id"; then
			csv_row "$dataset_id" "$dataset_name" "$assembly" "perl_not_annotating" "" "1"
			rm -rf "$dataset_dir/perl_scratch"
			return 1
		fi
		rm -rf "$dataset_dir/perl_scratch"
		;;
	esac
}

# Main dataset loop
# RECORD, DON'T ABORT. A failing cell is recorded (the measure_* helpers already
# write a canonicalize_failed or per-VCF-failure row via csv_row) and the loop
# moves on; a caller rejects the run by reading the CSV rows whose exit_code is
# non-zero.
#
# The `|| failed_datasets+=(...)` guard is load-bearing under `set -euo pipefail`.
# Both functions fail closed and return non-zero, so calling them bare would let
# the FIRST bad cell kill the script -- contradicting the "record the failure and
# skip" contract they document. A caller that also runs under `set -e` with no
# trap inherits that abort, so a single fail-closed cell would skip every
# remaining step of the surrounding run and produce no results at all.
failed_datasets=()
for entry in "${DATASETS[@]}"; do
	IFS=':' read -r dataset_id dataset_name assembly stype vcf_base species <<<"$entry"
	if ! should_run "$dataset_id" "$stype"; then
		echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] dataset $dataset_id skipped (datasets=$DATASETS_FILTER)"
		continue
	fi
	echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] dataset $dataset_id ($dataset_name $species $assembly $stype) starting"
	if [[ "$stype" == "snp_indel" ]]; then
		measure_p1 "$dataset_id" "$dataset_name" "$assembly" "$vcf_base" "$species" ||
			failed_datasets+=("$dataset_id")
	else
		measure_sv "$dataset_id" "$dataset_name" "$assembly" "$species" ||
			failed_datasets+=("$dataset_id")
	fi
done

# The wrapper's argument-vector log, one record per engine invocation, is folded into
# provenance.json after the last cell: under executed_argv for vep-rs and
# fastvep_executed_argv for fastVEP. A run writes at most one of the two logs.
ARGV_LOG="${VEP_ARGV_LOG:-$FASTVEP_ARGV_LOG}"
if [[ -n "$ARGV_LOG" ]]; then
	ARGV_KEY=executed_argv
	[[ "$ENGINE" == vep-rs ]] || ARGV_KEY=fastvep_executed_argv
	python3 - "$OUTPUT_DIR/provenance.json" "$ARGV_LOG" "$ARGV_KEY" <<'PY'
import json, sys
out, log, key = sys.argv[1:]
record = json.load(open(out))
record[key] = [line.split("\037")[:-1] for line in open(log).read().split("\n") if line]
json.dump(record, open(out, "w"), indent=1, sort_keys=True)
PY
	echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] $(grep -c . "$ARGV_LOG" || true) $ENGINE invocation(s) recorded in $OUTPUT_DIR/provenance.json"
fi

# container_probe <n>: n bare starts of the Ensembl VEP image, the same --user flag as the
# timed cells and no mounts, so the wall time is the container start-up alone. Every
# probe is timed: the image is already pulled and warm from the cells above, which is
# the state every timed Ensembl VEP invocation starts from. The rows carry the same
# provenance columns as wall_times.csv, keyed by probe index.
container_probe() {
	local n="$1" csv="$OUTPUT_DIR/container_start.csv" perf_dir="$OUTPUT_DIR/container_probe"
	mkdir -p "$perf_dir"
	echo "engine,run_id,host_id,replicate,image,probe_index,wall_seconds,exit_code" >"$csv"
	local i res wall rc
	for ((i = 1; i <= n; i++)); do
		res=$(time_run "container_probe_$i" "$perf_dir/probe-$i.time.txt" -- 			docker run --rm --user "$(id -u):$(id -g)" "$PERL_IMAGE" true)
		wall=${res%|*}
		rc=${res#*|}
		echo "$ENGINE,$RUN_ID,$HOST_ID,$REPLICATE,$PERL_IMAGE,$i,$wall,$rc" >>"$csv"
		if [[ "$rc" -ne 0 ]]; then
			echo "ERROR: [run_clone_measurement] container probe $i exited $rc; a probe with a non-zero exit must not enter an aggregate" >&2
		fi
	done
	echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] container probe: $n starts of $PERL_IMAGE written to $csv"
}

if [[ "$CONTAINER_PROBE" -gt 0 ]]; then
	container_probe "$CONTAINER_PROBE"
fi

if [[ ${#failed_datasets[@]} -gt 0 ]]; then
	echo "ERROR: [run_clone_measurement] ${#failed_datasets[@]} dataset(s) failed: ${failed_datasets[*]}" >&2
	echo "  Their rows in $WALL_CSV carry a non-zero rc; treat this run as failed." >&2
fi

echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] run_clone_measurement complete (engine=$ENGINE); wall_times.csv at $WALL_CSV"
