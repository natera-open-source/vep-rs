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

# make_cache.pl <cache dir> <info.txt lines> <with matrices: 1 or 0> [regulatory mode] [mapper mode]:
# one chromosome with a transcript and optional blessed regulatory fixture objects. The
# mapper mode `folded` gives the transcript the mapper object Ensembl's TranscriptMapper
# builds with a two-base deletion edit folded into its pairs, beside the exons the
# unedited pairs would come from; `empty` gives it a mapper object without pairs.
cat > "$WORK/make_cache.pl" <<'PERL'
use strict;
use warnings;
use File::Path qw(make_path);
use Storable qw(nstore);
my ($cache_dir, $info_lines, $with_matrices, $regulatory, $mapper_mode) = @ARGV;
$regulatory //= 'none';
$mapper_mode //= 'none';
make_path("$cache_dir/21") or die "make_path $cache_dir/21: $!";
open my $info, '>', "$cache_dir/info.txt" or die $!;
print $info $info_lines;
close $info;
open my $mfh, '<', 'matrix.gz' or die $!;
binmode $mfh;
my $blob = do { local $/; <$mfh> };
close $mfh;
my $matrix = sub { { analysis => $_[0], peptide_length => 3, translation_md5 => 'md5', matrix => $blob } };
my $unit = sub { bless { id => $_[0], start => $_[1], end => $_[2] }, 'Bio::EnsEMBL::Mapper::Unit' };
my $pair = sub {
    my ($cdna_start, $cdna_end, $gen_start, $gen_end) = @_;
    bless {
        from => $unit->('cdna', $cdna_start, $cdna_end),
        to   => $unit->('genome', $gen_start, $gen_end),
        ori  => 1,
    }, 'Bio::EnsEMBL::Mapper::Pair';
};
my @exons = (
    bless({ start => 1000, end => 1099, strand => 1, phase => -1, end_phase => -1, rank => 1 }, 'Bio::EnsEMBL::Exon'),
    bless({ start => 1200, end => 1299, strand => 1, phase => -1, end_phase => -1, rank => 2 }, 'Bio::EnsEMBL::Exon'),
);
# TranscriptMapper::_load_mapper over these exons with the deletion of cDNA 51-52
# folded in: the first exon splits at the edit and the cDNA shortens by two.
my @folded_pairs = ($pair->(1, 50, 1000, 1049), $pair->(51, 98, 1052, 1099), $pair->(99, 198, 1200, 1299));
my $mapper = $mapper_mode eq 'none' ? () : bless({
    start_phase => -1,
    exon_coord_mapper => bless({
        from => 'cdna', to => 'genomic', from_cs => undef, to_cs => undef, _is_sorted => 1,
        pair_count => ($mapper_mode eq 'folded' ? scalar @folded_pairs : 0),
        _pair_cdna    => { CDNA   => ($mapper_mode eq 'folded' ? [@folded_pairs] : []) },
        _pair_genomic => { GENOME => ($mapper_mode eq 'folded' ? [@folded_pairs] : []) },
    }, 'Bio::EnsEMBL::Mapper'),
}, 'Bio::EnsEMBL::TranscriptMapper');
my %transcript = (
    stable_id => 'ENST00000000001', version => 1, dbID => 1, start => 1000, end => 2000,
    strand => 1, biotype => 'protein_coding', source => 'ensembl',
    _gene_stable_id => 'ENSG00000000001', _gene_symbol => 'GENE1',
    _gene_symbol_source => 'HGNC', _gene_hgnc_id => 'HGNC:1',
    ($mapper_mode eq 'none' ? () : (
        _trans_exon_array => [@exons],
        attributes => [bless({ code => '_rna_edit', name => 'RNA Edit', value => '51 52 ' }, 'Bio::EnsEMBL::Attribute')],
    )),
    _variation_effect_feature_cache => {
        codon_table => 1,
        ($mapper_mode eq 'none' ? () : (sorted_exons => [@exons], mapper => $mapper)),
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
    my $bless = sub { bless $_[1], "Bio::EnsEMBL::Funcgen::$_[0]" };
    # The matrix as the dumper stores it on every motif that names it: elements by
    # position then nucleotide, the complexes with their components.
    my $maff = sub { $bless->('BindingMatrix', {
        stable_id => 'ENSPFM0363', name => 'MAFF_fixture', source => 'SELEX', threshold => 4.1,
        length => 2, unit => 'Frequencies',
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
    }) };
    my $motif = sub {
        my ($id, $dbid, $start, $end, $strand, $score, $seq, $matrix) = @_;
        $bless->('MotifFeature', {
            stable_id => $id, dbID => $dbid, start => $start, end => $end, strand => $strand,
            score => $score, _variation_effect_feature_cache => { seq => $seq }, binding_matrix => $matrix,
        });
    };
    if ($regulatory eq 'full') {
        $data = { '21' => {
            RegulatoryFeature => [
                $bless->('RegulatoryFeature', {
                    stable_id => 'ENSR21_BBX5K', dbID => '613467', start => '45000787', end => '45001450',
                    strand => '0', feature_type => 'enhancer', epigenome_count => 0,
                    cell_types => {
                        A549 => 'ACTIVE', HepG2 => 'INACTIVE', HUVEC => 'POISED', K562 => 'REPRESSED',
                        'adrenal_gland_(96-97_dpf)' => 'NA',
                    },
                }),
                $bless->('RegulatoryFeature', {
                    stable_id => 'ENSR_FIXTURE', dbID => 2, start => 45002000, end => 45002100,
                    strand => 0, feature_type => $bless->('FeatureType', { so_name => 'promoter' }),
                    cell_types => { A549 => 'INACTIVE' },
                }),
                $bless->('RegulatoryFeature', {
                    stable_id => 'ENSR_CTCF', dbID => 3, start => 45003000, end => 45003100,
                    strand => 0, feature_type => 'CTCF_binding_site', cell_types => {},
                }),
            ],
            MotifFeature => [
                $motif->('ENSM00000587149', '587149', '45009809', '45009823', '-1', 5.46194289642971, 'TTTGCTGAGTCAGCA', $maff->()),
                $motif->('ENSM_SHARED', 4, 45009830, 45009844, 1, 2.5, 'ACGTACGTACGTACG', $maff->()),
                $motif->('ENSM_LENGTH_ONLY', 5, 45009900, 45009905, 1, 1.5, 'ACGTAC',
                    $bless->('BindingMatrix', { stable_id => 'LENGTH_ONLY_MATRIX', length => 6 })),
            ],
        }};
    }
    if ($regulatory eq 'nomatrix') {
        $data = { '21' => { MotifFeature => [
            bless({ stable_id => 'ENSM_NO_MATRIX', dbID => 6, start => 45009900, end => 45009905, strand => 1,
                score => 1, _variation_effect_feature_cache => { seq => 'ACGTAC' } },
                'Bio::EnsEMBL::Funcgen::MotifFeature'),
        ] } };
    }
    if ($regulatory eq 'badactivity') {
        $data = { '21' => { RegulatoryFeature => [
            $bless->('RegulatoryFeature', {
                stable_id => 'ENSR_SIXTH', dbID => 7, start => 45002000, end => 45002100, strand => 0,
                feature_type => 'enhancer', cell_types => { A549 => 'ACTIVE', HepG2 => 'UNKNOWN' },
            }),
        ] } };
    }
    if ($regulatory eq 'matrixdiffers') {
        my $other = $maff->();
        $other->{elements}{2}{T} = 464;
        $data = { '21' => { MotifFeature => [
            $motif->('ENSM_FIRST', 8, 45009809, 45009823, 1, 1, 'TTTGCTGAGTCAGCA', $maff->()),
            $motif->('ENSM_SECOND', 9, 45009830, 45009844, 1, 1, 'ACGTACGTACGTACG', $other),
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

convert() { # convert <case> <info.txt lines> <with matrices> [regulatory mode] [mapper mode]
    local case=$1 lines=$2 matrices=$3 regulatory=${4:-none} mapper=${5:-none}
    run_perl make_cache.pl "$case" "$lines" "$matrices" "$regulatory" "$mapper" || { bad "$case: fixture cache built"; return 1; }
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
# says 0. A regulatory feature carries its activity as one code per name of the
# cell-type list (the sorted union of every feature's map; the input's own
# cell_types line is dropped), `-` where its map lacks the name; a motif carries
# its score, stranded sequence and the id of its matrix, and the matrices are
# written once, keyed by stable id, with their counts by position and the complex
# names TRANSCRIPTION_FACTORS prints.
if convert regulatory "${common}"$'regulatory\t0\ncell_types\tA549,HepG2\n' 0 full; then
    if python3 - "$WORK/regulatory.out" <<'PY'
import json
import pathlib
import sys
root = pathlib.Path(sys.argv[1])
info = json.loads((root / "info.json").read_text())
assert info["regulatory"] is True
assert info["regulatory_schema"] == 1
assert info["cell_types"] == ["A549", "HUVEC", "HepG2", "K562", "adrenal_gland_(96-97_dpf)"]
features = json.loads((root / "regulatory/21/45000001-46000000.json").read_text())
assert len(features) == 6
assert features[:3] == [
    dict(stable_id="ENSR21_BBX5K", dbid=613467, start=45000787, end=45001450,
         feature_type="RegulatoryFeature", biotype="enhancer", activity="APIRN"),
    dict(stable_id="ENSR_FIXTURE", dbid=2, start=45002000, end=45002100,
         feature_type="RegulatoryFeature", biotype="promoter", activity="I----"),
    dict(stable_id="ENSR_CTCF", dbid=3, start=45003000, end=45003100,
         feature_type="RegulatoryFeature", biotype="CTCF_binding_site", activity="-----"),
]
assert features[3] == dict(
    stable_id="ENSM00000587149", dbid=587149, start=45009809, end=45009823, strand=-1,
    feature_type="MotifFeature", score=5.46194289642971, seq="TTTGCTGAGTCAGCA",
    matrix="ENSPFM0363")
assert features[4] == dict(
    stable_id="ENSM_SHARED", dbid=4, start=45009830, end=45009844, strand=1,
    feature_type="MotifFeature", score=2.5, seq="ACGTACGTACGTACG", matrix="ENSPFM0363")
assert features[5] == dict(
    stable_id="ENSM_LENGTH_ONLY", dbid=5, start=45009900, end=45009905, strand=1,
    feature_type="MotifFeature", score=1.5, seq="ACGTAC", matrix="LENGTH_ONLY_MATRIX")
matrices = json.loads((root / "regulatory/binding_matrices.json").read_text())
assert matrices == {
    "ENSPFM0363": dict(name="MAFF_fixture", source="SELEX", threshold=4.1, length=2,
                       frequencies=[[137, 73, 62, 250], [22, 67, 1, 465]],
                       transcription_factors=["MAFK", "JUN::FOS"]),
    "LENGTH_ONLY_MATRIX": dict(length=6, frequencies=[], transcription_factors=[]),
}, matrices
assert not (root / "info.txt").exists()
PY
    then
        ok "regulatory objects, activity codes, motif sequences, the matrices table and metadata written"
    else
        bad "regulatory object conversion"
    fi
fi

# 5. Even an empty regulatory shard is a written tree, with an empty cell-type
# list and matrices table; an input marker without any source shards must not
# enable runtime regulatory lookup, and carries neither list nor schema.
if convert empty-regulatory "$common" 0 empty; then
    if python3 - "$WORK/empty-regulatory.out" <<'PY'
import json, pathlib, sys
root = pathlib.Path(sys.argv[1])
info = json.loads((root / "info.json").read_text())
assert info["regulatory"] is True and info["cell_types"] == [] and info["regulatory_schema"] == 1
assert json.loads((root / "regulatory/21/45000001-46000000.json").read_text()) == []
assert json.loads((root / "regulatory/binding_matrices.json").read_text()) == {}
PY
    then ok "empty regulatory shard enables the written tree"
    else bad "empty regulatory shard conversion"; fi
fi
if convert missing-regulatory "${common}"$'regulatory\t1\ncell_types\tA549\n' 0; then
    if python3 - "$WORK/missing-regulatory.out" <<'PY'
import json, pathlib, sys
root = pathlib.Path(sys.argv[1])
info = json.loads((root / "info.json").read_text())
assert not {"regulatory", "cell_types", "regulatory_schema"} & set(info), info
assert not (root / "regulatory").exists()
PY
    then ok "input regulatory marker without shards is dropped"
    else bad "input-only regulatory marker"; fi
fi

# 6. An activity value outside the five of the regulatory_activity table fails the
# conversion naming the feature and the value, as does a matrix that differs from
# an earlier copy under the same stable id; neither leaves metadata.
for case in 'badactivity:ENSR_SIXTH has activity .UNKNOWN. in HepG2' 'matrixdiffers:binding matrix ENSPFM0363 differs from its copy'; do
    mode=${case%%:*}; phrase=${case#*:}
    if run_perl make_cache.pl "$mode" "${common}" 0 "$mode"; then
        mkdir -p "$WORK/$mode.out"
        if run_perl storable_to_json.pl "$mode" "$mode.out" > "$WORK/$mode.log" 2>&1; then
            bad "$mode: conversion rejected"
        elif ! grep -q "$phrase" "$WORK/$mode.log"; then
            bad "$mode: refusal names the cause ($(tail -1 "$WORK/$mode.log"))"
        elif [[ -e "$WORK/$mode.out/info.json" || -e "$WORK/$mode.out/regulatory/binding_matrices.json" ]]; then
            bad "$mode: failed conversion leaves no metadata"
        else
            ok "$mode: rejected by name without metadata"
        fi
    else
        bad "$mode: fixture built"
    fi
done

# 7. A corrupt regulatory shard must fail conversion and leave no metadata,
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

# 8. A motif feature without a binding matrix is corrupt (the funcgen schema requires
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

# 9. edited_mapper_pairs_are_read_from_the_storable_mapper: the mapper pairs come from
# the Storable's Bio::EnsEMBL::Mapper object, which TranscriptMapper::_load_mapper wrote
# with the transcript's indel edits folded in, and from the exons only when the object
# carries no pair.
mapper_pairs() { # mapper_pairs <out dir>: the JSON pairs as `from_start-from_end:to_start-to_end` words
    python3 -c '
import json, sys
shard = json.load(open(sys.argv[1]))
pairs = shard[0]["variation_effect_feature_cache"]["mapper"]["pairs"]
print(" ".join("%d-%d:%d-%d" % (p["from_start"], p["from_end"], p["to_start"], p["to_end"]) for p in pairs))
' "$1/transcripts/21/1-1000000.json"
}
if convert folded-mapper "$common" 0 none folded; then
    pairs=$(mapper_pairs "$WORK/folded-mapper.out")
    if [[ $pairs == "1-50:1000-1049 51-98:1052-1099 99-198:1200-1299" ]]; then
        ok "edited_mapper_pairs_are_read_from_the_storable_mapper: the folded pairs, not the exons' ($pairs)"
    else
        bad "edited_mapper_pairs_are_read_from_the_storable_mapper: $pairs"
    fi
fi
if convert empty-mapper "$common" 0 none empty; then
    pairs=$(mapper_pairs "$WORK/empty-mapper.out")
    if [[ $pairs == "1-100:1000-1099 101-200:1200-1299" ]]; then
        ok "edited_mapper_pairs_are_read_from_the_storable_mapper: a mapper object without pairs falls back to the exons ($pairs)"
    else
        bad "edited_mapper_pairs_are_read_from_the_storable_mapper: fallback $pairs"
    fi
fi

# The Ensembl VEP test caches (t/testdata/cache of an ensembl-vep checkout, release
# 84 on GRCh38, chromosome 21 region 25000001-26000000): the RefSeq cache carries 86
# transcripts, three of them BAM-edited `ok` with four `_rna_edit` attributes, and
# the merged cache 145 transcripts stamped `_source_cache`. Read from
# $VEP_ENSEMBL_VEP_DIR, else from an ensembl-vep checkout beside this repository;
# without one the checks below are skipped, not failed.
ENSEMBL_VEP_DIR="${VEP_ENSEMBL_VEP_DIR:-$HERE/../../../ensembl-vep}"
TEST_CACHES="$ENSEMBL_VEP_DIR/t/testdata/cache"
REGION=25000001-26000000
cat > "$WORK/check_mapper.pl" <<'PERL'
# check_mapper.pl <Storable shard> <JSON shard>: every transcript's JSON pairs equal
# the Storable mapper object's cDNA pairs; the pair lengths sum to the exon lengths
# plus the indel edits TranscriptMapper::_load_mapper folds (those starting at or
# before the cDNA end of an exon, shifted by the edits before them); an unedited
# transcript's pairs equal the exon-derived ones; an `ok` transcript keeps its status
# and its cached cDNA. Prints the counts checked.
use strict;
use warnings;
use JSON;
use Storable qw(fd_retrieve);
my ($storable, $json_path) = @ARGV;
open my $fh, '-|', 'gzip', '-dc', $storable or die "gzip: $!";
my $data = fd_retrieve($fh);
close $fh;
open my $jfh, '<', $json_path or die "$json_path: $!";
# Keyed by stable id and dbID: two features of one gene can share a stable id.
my %json_by_id = map { ("$_->{stable_id}/$_->{dbID}", $_) } @{ JSON->new->decode(do { local $/; <$jfh> }) };
close $jfh;
my ($checked, $edited_ok, $unedited) = (0, 0, 0);
for my $chr (keys %$data) {
    for my $t (@{ $data->{$chr} }) {
        my $id = "$t->{stable_id}/$t->{dbID}";
        my $out = $json_by_id{$id} or die "$id missing from the JSON shard\n";
        my $vefc = $t->{_variation_effect_feature_cache};
        my @object = sort { $a->{from}{start} <=> $b->{from}{start} }
            @{ $vefc->{mapper}{exon_coord_mapper}{_pair_cdna}{CDNA} };
        my @json = @{ $out->{variation_effect_feature_cache}{mapper}{pairs} };
        die "$id: " . scalar(@json) . " JSON pairs, " . scalar(@object) . " in the Storable\n" unless @json == @object;
        for my $i (0 .. $#object) {
            my ($o, $j) = ($object[$i], $json[$i]);
            die "$id: pair $i differs\n" unless $j->{from_start} == $o->{from}{start}
                && $j->{from_end} == $o->{from}{end} && $j->{to_start} == $o->{to}{start}
                && $j->{to_end} == $o->{to}{end} && $j->{ori} == $o->{ori};
        }
        my @exons = @{ $t->{_trans_exon_array} };
        my @edits = sort { $a->[0] <=> $b->[0] }
            map { my ($s, $e, $alt) = split ' ', $_->{value}, 3; $alt //= ''; [$s, $e, length($alt) - ($e - $s + 1)] }
            grep { $_->{code} eq '_rna_edit' } @{ $t->{attributes} || [] };
        my $had_edits = scalar @edits;
        my ($exon_sum, $cdna_end, $shift) = (0, 0, 0);
        for my $exon (@exons) {
            my $length = $exon->{end} - $exon->{start} + 1;
            $exon_sum += $length;
            $cdna_end += $length;
            while (@edits && $edits[0][0] + $shift <= $cdna_end) {
                my $edit = shift @edits;
                next unless $edit->[2];
                $cdna_end += $edit->[2];
                $shift += $edit->[2];
            }
        }
        my $pair_sum = 0;
        $pair_sum += $_->{from_end} - $_->{from_start} + 1 for @json;
        die "$id: pair lengths sum to $pair_sum, exons to $exon_sum, folded edits shift $shift\n"
            unless $pair_sum == $exon_sum + $shift;
        if (($t->{_bam_edit_status} // '') eq 'ok') {
            die "$id: bam_edit_status not carried\n" unless ($out->{bam_edit_status} // '') eq 'ok';
            die "$id: spliced_seq not carried\n" unless defined $out->{variation_effect_feature_cache}{spliced_seq};
            $edited_ok++;
        }
        unless ($had_edits) {
            my @sorted = sort { $t->{strand} < 0 ? $b->{start} <=> $a->{start} : $a->{start} <=> $b->{start} } @exons;
            my ($pos, @derived) = (1);
            for my $exon (@sorted) {
                my $length = $exon->{end} - $exon->{start} + 1;
                push @derived, [$pos, $pos + $length - 1, $exon->{start}, $exon->{end}];
                $pos += $length;
            }
            die "$id: unedited pairs differ from the exon-derived ones\n"
                unless @derived == @json && !grep {
                    my ($d, $j) = ($derived[$_], $json[$_]);
                    $d->[0] != $j->{from_start} || $d->[1] != $j->{from_end} || $d->[2] != $j->{to_start} || $d->[3] != $j->{to_end}
                } 0 .. $#derived;
            $unedited++;
        }
        $checked++;
    }
}
print "checked=$checked edited_ok=$edited_ok unedited=$unedited\n";
PERL
if [[ -f "$TEST_CACHES/homo_sapiens_refseq/84_GRCh38/21/$REGION.gz" ]]; then
    for set in homo_sapiens_refseq homo_sapiens_merged homo_sapiens; do
        mkdir -p "$WORK/$set/84_GRCh38" "$WORK/$set/84_GRCh38.out"
        cp "$TEST_CACHES/$set/84_GRCh38/info.txt" "$WORK/$set/84_GRCh38/"
        cp -r "$TEST_CACHES/$set/84_GRCh38/21" "$WORK/$set/84_GRCh38/"
        if run_perl storable_to_json.pl "$set/84_GRCh38" "$set/84_GRCh38.out" 21 > "$WORK/$set.log" 2>&1; then
            ok "$set/84_GRCh38: converter exit 0"
        else
            bad "$set/84_GRCh38: converter failed ($(tail -1 "$WORK/$set.log"))"
        fi
    done
    counts=$(run_perl check_mapper.pl "homo_sapiens_refseq/84_GRCh38/21/$REGION.gz" "homo_sapiens_refseq/84_GRCh38.out/transcripts/21/$REGION.json" 2>&1)
    if [[ $counts == "checked=86 edited_ok=3 unedited=83" ]]; then
        ok "edited_mapper_pairs_are_read_from_the_storable_mapper: the RefSeq test region ($counts)"
    else
        bad "edited_mapper_pairs_are_read_from_the_storable_mapper: the RefSeq test region: $counts"
    fi
else
    echo "  skip edited_mapper_pairs_are_read_from_the_storable_mapper on the RefSeq test region: no $TEST_CACHES (set VEP_ENSEMBL_VEP_DIR)"
fi

# 10. refseq_and_merged_directories_convert_with_their_source_type: info.json declares the
# transcript set from the species directory's suffix and carries the `bam` line.
source_type_state() { # source_type_state <out dir>: `source_type=<value or ->` and `bam=<value or ->`
    python3 -c '
import json, sys
info = json.load(open(sys.argv[1]))
print("source_type=%s bam=%s" % (info.get("source_type", "-"), info.get("bam", "-")))
' "$1/info.json"
}
bam_line=$'bam\tinterim_GRCh38.p10_knownrefseq_alignments_2017-01-13.bam\n'
for set in homo_sapiens_refseq homo_sapiens_merged homo_sapiens; do
    lines=$common
    [[ $set == homo_sapiens_refseq ]] && lines="${common}${bam_line}"
    if convert "sets/$set/84_GRCh38" "$lines" 0; then
        state=$(source_type_state "$WORK/sets/$set/84_GRCh38.out")
        want="source_type=${set#homo_sapiens_} bam=-"
        [[ $set == homo_sapiens ]] && want="source_type=ensembl bam=-"
        [[ $set == homo_sapiens_refseq ]] && want="source_type=refseq bam=interim_GRCh38.p10_knownrefseq_alignments_2017-01-13.bam"
        if [[ $state == "$want" ]]; then
            ok "refseq_and_merged_directories_convert_with_their_source_type: $set ($state)"
        else
            bad "refseq_and_merged_directories_convert_with_their_source_type: $set: $state"
        fi
    fi
done
if [[ -f "$WORK/polyphen-only.out/info.json" ]]; then
    state=$(source_type_state "$WORK/polyphen-only.out")
    if [[ $state == "source_type=- bam=-" ]]; then
        ok "refseq_and_merged_directories_convert_with_their_source_type: a directory outside the layout declares none"
    else
        bad "refseq_and_merged_directories_convert_with_their_source_type: outside the layout: $state"
    fi
fi
if [[ -f "$WORK/homo_sapiens_merged/84_GRCh38.out/info.json" ]]; then
    if python3 - "$WORK" <<'PY'
import json, pathlib, sys
root = pathlib.Path(sys.argv[1])
assert json.loads((root / "homo_sapiens_refseq/84_GRCh38.out/info.json").read_text())["source_type"] == "refseq"
assert json.loads((root / "homo_sapiens_refseq/84_GRCh38.out/info.json").read_text())["bam"] == "interim_GRCh38.p10_knownrefseq_alignments_2017-01-13.bam"
merged_info = json.loads((root / "homo_sapiens_merged/84_GRCh38.out/info.json").read_text())
assert merged_info["source_type"] == "merged" and "bam" not in merged_info
ensembl_info = json.loads((root / "homo_sapiens/84_GRCh38.out/info.json").read_text())
assert ensembl_info["source_type"] == "ensembl" and "bam" not in ensembl_info
merged = json.loads((root / "homo_sapiens_merged/84_GRCh38.out/transcripts/21/25000001-26000000.json").read_text())
sources = [t.get("source_cache") for t in merged]
assert len(merged) == 145 and sources.count("Ensembl") == 70 and sources.count("RefSeq") == 75, (len(merged), sources.count("Ensembl"), sources.count("RefSeq"))
refseq = json.loads((root / "homo_sapiens_refseq/84_GRCh38.out/transcripts/21/25000001-26000000.json").read_text())
assert len(refseq) == 86 and not any("source_cache" in t for t in refseq)
assert sorted(t["bam_edit_status"] for t in refseq if "bam_edit_status" in t) == ["ok", "ok", "ok"]
PY
    then
        ok "refseq_and_merged_directories_convert_with_their_source_type: the test caches (merged 70 Ensembl, 75 RefSeq; RefSeq bam line and 3 ok)"
    else
        bad "refseq_and_merged_directories_convert_with_their_source_type: the test caches"
    fi
fi

# 11. real_regulatory_shard_converts_to_the_documented_values: the release 116 GRCh38
# shard 21/44000001-45000000_reg.gz (273 regulatory and 963 motif features, read from
# $VEP_REGULATORY_SHARD) converts to the values its Storable records hold: the first
# feature's activity over the 112 epigenomes of the build (77 of whose names contain a
# comma), the motif ENSM00000586186 and its matrix ENSPFM0364. Skipped, not failed,
# when the variable is unset.
if [[ -n ${VEP_REGULATORY_SHARD:-} && -f $VEP_REGULATORY_SHARD ]]; then
    if run_perl make_cache.pl real-shard "${common}"$'regulatory\t1\nsource_regbuild\t1.0\n' 0; then
        cp "$VEP_REGULATORY_SHARD" "$WORK/real-shard/21/44000001-45000000_reg.gz"
        mkdir -p "$WORK/real-shard.out"
        if run_perl storable_to_json.pl real-shard real-shard.out > "$WORK/real-shard.log" 2>&1 \
            && python3 - "$WORK/real-shard.out" <<'PY'
import collections, json, pathlib, sys
root = pathlib.Path(sys.argv[1])
info = json.loads((root / "info.json").read_text())
names = info["cell_types"]
assert len(names) == 112 and sum("," in name for name in names) == 77, (len(names), names)
assert names == sorted(names) and info["regulatory_schema"] == 1 and info["source_versions"]["regbuild"] == "1.0"
features = json.loads((root / "regulatory/21/44000001-45000000.json").read_text())
kinds = collections.Counter(feature["feature_type"] for feature in features)
assert kinds == {"RegulatoryFeature": 273, "MotifFeature": 963}, kinds
first = features[0]
assert first["stable_id"].startswith("ENSR") and len(first["activity"]) == 112
assert collections.Counter(first["activity"]) == {"A": 58, "I": 54}, first["activity"]
assert all(len(feature["activity"]) == 112 for feature in features if feature["feature_type"] == "RegulatoryFeature")
motif = next(feature for feature in features if feature["stable_id"] == "ENSM00000586186")
assert motif == dict(stable_id="ENSM00000586186", dbid=586186, start=44006279, end=44006290, strand=-1,
                     feature_type="MotifFeature", score=10.6784360161335, seq="AATAATGCTGAG",
                     matrix="ENSPFM0364"), motif
matrices = json.loads((root / "regulatory/binding_matrices.json").read_text())
assert len(matrices) == 44, len(matrices)
matrix = matrices["ENSPFM0364"]
assert matrix["length"] == 12 and matrix["threshold"] == 4.4 and matrix["source"] == "SELEX"
assert matrix["frequencies"][0] == [4881, 1573, 1789, 1847] and len(matrix["frequencies"]) == 12
assert matrix["transcription_factors"] == ["MAFK"], matrix
PY
        then
            ok "real_regulatory_shard_converts_to_the_documented_values"
        else
            bad "real_regulatory_shard_converts_to_the_documented_values ($(tail -1 "$WORK/real-shard.log"))"
        fi
    else
        bad "real_regulatory_shard_converts_to_the_documented_values: fixture built"
    fi
else
    echo "  skip real_regulatory_shard_converts_to_the_documented_values: VEP_REGULATORY_SHARD is unset"
fi

echo "$PASS ok, $FAIL failed"
[[ $FAIL -eq 0 ]]
