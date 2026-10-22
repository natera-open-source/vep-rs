#!/usr/bin/env bash
# storable_to_json.pl keeps a `sift` or `polyphen` marker in info.json only when the
# conversion wrote a matrix of the tool and info.txt names the tool's `source_<tool>`
# version, the one line VEP's availability check reads (Cache/Transcript.pm:155). This
# suite builds a one-transcript Storable cache three ways, converts each and reads the
# info.json back: both markers with a version for PolyPhen alone keeps `polyphen` and
# drops `sift`; a version for both keeps both; no matrices drops both.
#
# The converter needs a perl with Storable and JSON. The suite runs the first perl on
# PATH when it has them, else the perl of the ensemblorg/ensembl-vep:release_116.2
# image when docker holds it, and exits 2 with neither (a skip, not a pass). The host
# needs gzip, for the fixture's matrix, and python3, to read the info.json back.
#
# Run: bash scripts/data/test_storable_to_json.sh
set -uo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
CONVERTER="$HERE/storable_to_json.pl"
IMAGE=ensemblorg/ensembl-vep:release_116.2

PASS=0; FAIL=0
ok() { PASS=$((PASS + 1)); echo "  ok   $1"; }
bad() { FAIL=$((FAIL + 1)); echo "  FAIL $1"; }

for tool in gzip python3; do
    command -v "$tool" >/dev/null 2>&1 || { echo "ERROR: [test_storable_to_json] $tool is not on PATH" >&2; exit 2; }
done

WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT
cp "$CONVERTER" "$WORK/storable_to_json.pl"

# run_perl <script under $WORK> <args...>: the script and its paths are relative to $WORK.
if perl -MStorable -MJSON -e 1 2>/dev/null; then
    run_perl() { (cd "$WORK" && perl "$@"); }
    echo "perl: $(command -v perl)"
elif command -v docker >/dev/null 2>&1 && docker image inspect "$IMAGE" >/dev/null 2>&1; then
    run_perl() {
        docker run --rm --pull never -u "$(id -u):$(id -g)" -v "$WORK":/work -w /work "$IMAGE" perl "$@"
    }
    echo "perl: $IMAGE"
else
    echo "ERROR: [test_storable_to_json] no perl with Storable and JSON on PATH and no $IMAGE in docker" >&2
    exit 2
fi

# A gzip-compressed blob stands in for a prediction matrix: the converter re-encodes a
# matrix as it finds it when the bytes already carry the gzip magic.
printf 'VEP' | gzip -c > "$WORK/matrix.gz"

# make_cache.pl <cache dir> <info.txt lines> <with matrices: 1 or 0>: one chromosome
# directory with one Storable shard holding one transcript.
cat > "$WORK/make_cache.pl" <<'PERL'
use strict;
use warnings;
use Storable qw(nstore);
my ($cache_dir, $info_lines, $with_matrices) = @ARGV;
mkdir $cache_dir or die "mkdir $cache_dir: $!";
mkdir "$cache_dir/21" or die "mkdir $cache_dir/21: $!";
open my $info, '>', "$cache_dir/info.txt" or die $!;
print $info $info_lines;
close $info;
open my $mfh, '<', 'matrix.gz' or die $!;
binmode $mfh;
my $blob = do { local $/; <$mfh> };
close $mfh;
my $matrix = sub { { analysis => $_[0], peptide_length => 3, translation_md5 => 'md5', matrix => $blob } };
my %transcript = (
    stable_id => 'ENST00000000001', version => 1, dbID => 1, start => 1000, end => 2000,
    strand => 1, biotype => 'protein_coding', source => 'ensembl',
    _gene_stable_id => 'ENSG00000000001', _gene_symbol => 'GENE1',
    _gene_symbol_source => 'HGNC', _gene_hgnc_id => 'HGNC:1',
    _variation_effect_feature_cache => {
        codon_table => 1,
        ($with_matrices ? (protein_function_predictions => {
            sift => $matrix->('sift'),
            polyphen_humvar => $matrix->('polyphen'),
        }) : ()),
    },
);
my $shard = "$cache_dir/21/1-1000000";
nstore({ '21' => [\%transcript] }, $shard) or die "nstore: $!";
system('gzip', '-f', $shard) == 0 or die "gzip: $?";
PERL

# marker_state <out dir>: prints `sift=<value or ->` and `polyphen=<value or ->` and the
# source_versions keys of info.json.
marker_state() {
    python3 -c '
import json, sys
info = json.load(open(sys.argv[1]))
print("sift=%s polyphen=%s versions=%s" % (
    info.get("sift", "-"), info.get("polyphen", "-"),
    ",".join(sorted(info.get("source_versions", {})))))
' "$1/info.json"
}

convert() { # convert <case> <info.txt lines> <with matrices>
    local case=$1 lines=$2 matrices=$3
    run_perl make_cache.pl "$case" "$lines" "$matrices" || { bad "$case: fixture cache built"; return 1; }
    mkdir -p "$WORK/$case.out"
    if run_perl storable_to_json.pl "$case" "$case.out" > "$WORK/$case.log" 2>&1; then
        ok "$case: converter exit 0"
    else
        bad "$case: converter failed ($(tail -1 "$WORK/$case.log"))"
        return 1
    fi
    if [[ -f "$WORK/$case.out/info.json" ]]; then
        ok "$case: info.json written"
    else
        bad "$case: info.json written"
        return 1
    fi
}

common=$'species\thomo_sapiens\nassembly\tGRCh38\nsift\tb\npolyphen\tb\n'
polyphen_version=$'source_polyphen\t2.2.3\n'
sift_version=$'source_sift\t6.2.1\n'

# 1. Both markers, a version line for PolyPhen alone: `sift` is dropped, `polyphen` kept.
if convert polyphen-only "${common}${polyphen_version}" 1; then
    state=$(marker_state "$WORK/polyphen-only.out")
    if [[ $state == "sift=- polyphen=b versions=polyphen" ]]; then
        ok "marker without a source_<tool> line is dropped ($state)"
    else
        bad "marker without a source_<tool> line: $state"
    fi
fi

# 2. A version line for both tools: both markers stay.
if convert both-versions "${common}${polyphen_version}${sift_version}" 1; then
    state=$(marker_state "$WORK/both-versions.out")
    if [[ $state == "sift=b polyphen=b versions=polyphen,sift" ]]; then
        ok "markers with their source_<tool> lines stay ($state)"
    else
        bad "markers with their source_<tool> lines: $state"
    fi
fi

# 3. Version lines for both tools but no matrix written: both markers are dropped.
if convert no-matrices "${common}${polyphen_version}${sift_version}" 0; then
    state=$(marker_state "$WORK/no-matrices.out")
    if [[ $state == "sift=- polyphen=- versions=polyphen,sift" ]]; then
        ok "markers without a written matrix are dropped ($state)"
    else
        bad "markers without a written matrix: $state"
    fi
fi

echo "$PASS ok, $FAIL failed"
[[ $FAIL -eq 0 ]]
