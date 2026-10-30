#!/usr/bin/env bash
# storable_to_json.pl keeps a `sift` or `polyphen` marker in info.json only when the
# conversion wrote a matrix of the tool and info.txt names the tool's `source_<tool>`
# version, the one line VEP's availability check reads (Cache/Transcript.pm:155). This
# script builds small Storable caches, converts them and checks prediction markers
# and regulatory feature shards. Regulatory metadata must describe the output tree,
# never just repeat the input marker.
#
# The converter needs a perl with Storable and JSON. The script runs the first perl on
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

# make_cache.pl <cache dir> <info.txt lines> <with matrices: 1 or 0> [regulatory mode]:
# one chromosome with a transcript and optional blessed regulatory fixture objects.
cat > "$WORK/make_cache.pl" <<'PERL'
use strict;
use warnings;
use Storable qw(nstore);
my ($cache_dir, $info_lines, $with_matrices, $regulatory) = @ARGV;
$regulatory //= 'none';
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
if ($regulatory ne 'none') {
    my $data = {};
    if ($regulatory eq 'full') {
        my $bless = sub { bless $_[1], "Bio::EnsEMBL::Funcgen::$_[0]" };
        $data = { '21' => {
            RegulatoryFeature => [
                $bless->('RegulatoryFeature', {
                    stable_id => 'ENSR21_BBX5K', start => '45000787', end => '45001450',
                    strand => '0', feature_type => 'enhancer', cell_types => { A549 => 'ACTIVE' },
                }),
                $bless->('RegulatoryFeature', {
                    stable_id => 'ENSR_FIXTURE', start => 45002000, end => 45002100,
                    strand => 0, feature_type => $bless->('FeatureType', { so_name => 'promoter' }),
                }),
            ],
            MotifFeature => [
                $bless->('MotifFeature', {
                    stable_id => 'ENSM00000587149', start => '45009809', end => '45009823',
                    strand => '-1', score => 5.46194289642971,
                    binding_matrix => $bless->('BindingMatrix', {
                        stable_id => 'ENSPFM0363', name => 'MAFF_fixture',
                        elements => {
                            1 => { A => 137, C => 73, G => 62, T => 250 },
                            2 => { A => 22, C => 67, G => 1, T => 465 },
                        },
                        associated_transcription_factor_complexes => [
                            $bless->('TranscriptionFactorComplex', {
                                display_name => 'MAFK',
                                components => [$bless->('TranscriptionFactor', { name => 'component_not_output' })],
                            }),
                            $bless->('TranscriptionFactorComplex', { display_name => 'JUN::FOS' }),
                        ],
                    }),
                }),
                $bless->('MotifFeature', {
                    stable_id => 'ENSM_ARRAY', start => 45009830, end => 45009831,
                    strand => 1, binding_matrix => $bless->('BindingMatrix', {
                        stable_id => 'ARRAY_MATRIX', elements => [[1, 2, 3, 4], [4, 3, 2, 1]],
                    }),
                }),
                $bless->('MotifFeature', {
                    stable_id => 'ENSM_GFF', start => 45009900, end => 45009905,
                    strand => 1, binding_matrix => $bless->('BindingMatrix', { stable_id => 'GFF_MATRIX' }),
                }),
            ],
        }};
    }
    if ($regulatory eq 'nomatrix') {
        $data = { '21' => { MotifFeature => [
            bless({ stable_id => 'ENSM_NO_MATRIX', start => 45009900, end => 45009905, strand => 1 },
                'Bio::EnsEMBL::Funcgen::MotifFeature'),
        ] } };
    }
    my $reg_shard = "$cache_dir/21/45000001-46000000_reg";
    if ($regulatory eq 'corrupt') {
        open my $bad, '>', $reg_shard or die $!;
        print $bad 'not a Storable image';
        close $bad;
    } else {
        nstore($data, $reg_shard) or die "nstore regulatory: $!";
    }
    system('gzip', '-f', $reg_shard) == 0 or die "gzip regulatory: $?";
}
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

convert() { # convert <case> <info.txt lines> <with matrices> [regulatory mode]
    local case=$1 lines=$2 matrices=$3 regulatory=${4:-none}
    run_perl make_cache.pl "$case" "$lines" "$matrices" "$regulatory" || { bad "$case: fixture cache built"; return 1; }
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

# 4. Real VEP object shapes: metadata is enabled by written shards even if input
# says 0; original matrix counts/positions and complex names survive unchanged.
if convert regulatory "${common}"$'regulatory\t0\ncell_types\tA549\n' 0 full; then
    if python3 - "$WORK/regulatory.out" <<'PY'
import json
import pathlib
import sys
root = pathlib.Path(sys.argv[1])
info = json.loads((root / "info.json").read_text())
assert info["regulatory"] is True
assert "cell_types" not in info
features = json.loads((root / "regulatory/21/45000001-46000000.json").read_text())
assert len(features) == 5
assert features[:2] == [
    dict(stable_id="ENSR21_BBX5K", start=45000787, end=45001450, strand=0,
         feature_type="RegulatoryFeature", biotype="enhancer"),
    dict(stable_id="ENSR_FIXTURE", start=45002000, end=45002100, strand=0,
         feature_type="RegulatoryFeature", biotype="promoter"),
]
assert features[2] == dict(
    stable_id="ENSM00000587149", start=45009809, end=45009823, strand=-1,
    feature_type="MotifFeature", binding_matrix=dict(
        stable_id="ENSPFM0363", name="MAFF_fixture",
        transcription_factors=["MAFK", "JUN::FOS"],
        frequencies={"1": dict(A=137, C=73, G=62, T=250),
                     "2": dict(A=22, C=67, G=1, T=465)}))
assert features[3]["binding_matrix"] == dict(
    stable_id="ARRAY_MATRIX", transcription_factors=[],
    frequencies=[[1, 2, 3, 4], [4, 3, 2, 1]])
assert features[4]["binding_matrix"] == dict(stable_id="GFF_MATRIX", transcription_factors=[])
assert not (root / "info.txt").exists()
PY
    then
        ok "regulatory objects, raw matrices, TF complex names and metadata preserved"
    else
        bad "regulatory object conversion"
    fi
fi

# 5. Even an empty regulatory shard is a written tree; an input marker without
# any source shards must not enable runtime regulatory lookup.
if convert empty-regulatory "$common" 0 empty; then
    if python3 - "$WORK/empty-regulatory.out" <<'PY'
import json, pathlib, sys
root = pathlib.Path(sys.argv[1])
assert json.loads((root / "info.json").read_text())["regulatory"] is True
assert json.loads((root / "regulatory/21/45000001-46000000.json").read_text()) == []
PY
    then ok "empty regulatory shard enables the written tree"
    else bad "empty regulatory shard conversion"; fi
fi
if convert missing-regulatory "${common}"$'regulatory\t1\n' 0; then
    if python3 - "$WORK/missing-regulatory.out" <<'PY'
import json, pathlib, sys
root = pathlib.Path(sys.argv[1])
assert "regulatory" not in json.loads((root / "info.json").read_text())
assert not (root / "regulatory").exists()
PY
    then ok "input regulatory marker without shards is dropped"
    else bad "input-only regulatory marker"; fi
fi

# 6. A corrupt regulatory shard must fail conversion and leave no metadata,
# including when an earlier conversion left a positive marker in the output.
if run_perl make_cache.pl corrupt-regulatory "${common}"$'regulatory\t1\n' 0 corrupt; then
    mkdir -p "$WORK/corrupt-regulatory.out"
    printf '{"regulatory":true}\n' > "$WORK/corrupt-regulatory.out/info.json"
    if run_perl storable_to_json.pl corrupt-regulatory corrupt-regulatory.out > "$WORK/corrupt-regulatory.log" 2>&1; then
        bad "corrupt regulatory shard rejected"
    elif [[ -e "$WORK/corrupt-regulatory.out/info.json" || -e "$WORK/corrupt-regulatory.out/regulatory/21/45000001-46000000.json" ]]; then
        bad "failed regulatory conversion leaves no metadata or partial shard"
    else
        ok "corrupt regulatory shard rejected without metadata or partial JSON"
    fi
else
    bad "corrupt regulatory fixture built"
fi

# 7. A motif feature without a binding matrix is corrupt (the funcgen schema requires
# the matrix and the GFF route blesses a stub), so the conversion fails and names it.
if run_perl make_cache.pl no-matrix "${common}" 0 nomatrix; then
    mkdir -p "$WORK/no-matrix.out"
    if run_perl storable_to_json.pl no-matrix no-matrix.out > "$WORK/no-matrix.log" 2>&1; then
        bad "motif without a binding matrix rejected"
    elif ! grep -q 'MotifFeature ENSM_NO_MATRIX lacks binding_matrix' "$WORK/no-matrix.log"; then
        bad "motif without a binding matrix named ($(tail -1 "$WORK/no-matrix.log"))"
    elif [[ -e "$WORK/no-matrix.out/info.json" ]]; then
        bad "failed conversion leaves no metadata"
    else
        ok "motif without a binding matrix rejected by name"
    fi
else
    bad "matrix-less motif fixture built"
fi

echo "$PASS ok, $FAIL failed"
[[ $FAIL -eq 0 ]]
