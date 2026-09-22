#!/usr/bin/perl
# Storable<->JSON cache parity audit (N=1000 transcripts).
#
# Reads N transcripts from both Storable + JSON caches, recursively diffs every
# key path, and classifies each Storable field against the conversion rules of
# scripts/data/storable_to_json.pl:
#   identical            same path, same container type on both sides
#   renamed              present on the JSON side under the converter's renamed path
#   converted            present but as a different container type (object -> string)
#   dropped_undef        undefined on the Storable side, so the converter writes nothing
#   dropped_by_converter absent by a named converter rule (back-references, the slice,
#                        the exon coordinate mapper, prediction matrices, ...)
#   unexplained          absent on the JSON side and matched by no rule: a parity gap
# JSON fields with no Storable counterpart are `derived` when a rule produces them and
# `only_in_json` otherwise.
#
# Rationale: a null `vefc.three_prime_utr` in the JSON cache changes which
# consequences vep-rs emits, so a field the converter silently drops shows up as
# a concordance difference rather than as an error. This confirms no other field
# is similarly missing, since field-coverage gaps invalidate F1 numbers.
#
# Usage:
#   perl cache_parity_audit.pl <storable_cache_dir> <json_cache_dir> [N] [SEED]
# Where:
#   storable_cache_dir = /opt/vep/.vep (or the Perl cache root, with homo_sapiens/<version>/)
#   json_cache_dir     = directory containing chr-keyed JSON files
#   N                  = number of transcripts to sample (default 1000)
#   SEED               = srand seed for the chunk order (default 1), recorded in the
#                        manifest so a run is reproducible
#
# Output: a JSON manifest on stdout carrying one count map per class above.

use strict;
use warnings;
use Storable qw(fd_retrieve);
use JSON::PP;
use File::Find;
use Scalar::Util qw(blessed reftype refaddr);

my $storable_root = $ARGV[0] or die "Usage: $0 storable_cache_dir json_cache_dir [N]\n";
my $json_root     = $ARGV[1] or die "Usage: $0 storable_cache_dir json_cache_dir [N]\n";
my $N             = $ARGV[2] // 1000;
my $SEED          = $ARGV[3] // 1;
srand($SEED);

# --- Discover Storable transcript chunks ---
my @storable_files;
find(
    sub { push @storable_files, $File::Find::name if -f && /\.gz$/ && !/_(?:reg|var)\.gz$/; },
    $storable_root
);
die "No .gz Storable chunks found under $storable_root\n" unless @storable_files;

# --- Sample N transcripts from the Storable chunks ---
# The chunk files are ordered by a randomising comparator, not a shuffle, and the
# quota is filled in that order, so the sample comes from however few leading chunks
# reach N and is not a uniform draw over the cache. That is sufficient for a
# field-presence audit, which is all this script claims.
my @sampled;       # list of [chr, transcript_hashref, source_file]
my $remaining = $N;
my @shuffled = sort { rand(1) <=> rand(1) } @storable_files;

for my $file (@shuffled) {
    last if $remaining <= 0;
    open(my $zfh, "-|", "gunzip -c '$file'") or next;
    binmode $zfh;
    my $chunk = eval { fd_retrieve($zfh) };
    close $zfh;
    next unless ref $chunk eq 'HASH';
    my ($chr) = ($file =~ m{/(\w+)/[^/]+\.gz$});
    $chr //= "unknown";
    for my $entry (values %$chunk) {
        my @transcripts;
        if (ref $entry eq 'ARRAY') {
            push @transcripts, grep { ref($_) && exists($_->{stable_id}) } @$entry;
        } elsif (ref $entry eq 'HASH' && exists $entry->{stable_id}) {
            push @transcripts, $entry;
        }
        for my $tr (@transcripts) {
            last if $remaining <= 0;
            push @sampled, [$chr, $tr, $file];
            $remaining--;
        }
        last if $remaining <= 0;
    }
}

# --- Load JSON side: lookup by stable_id ---
my %json_by_id;
for my $jfile (sort glob("$json_root/*.json")) {
    open(my $jfh, "<", $jfile) or next;
    my $payload = do { local $/; <$jfh> };
    close $jfh;
    my $data = eval { decode_json($payload) };
    next unless ref $data;
    # Expected structure: hash keyed by chr -> array of transcripts (or hash by region).
    my $iter = sub {
        my ($v) = @_;
        if (ref $v eq 'ARRAY') {
            for my $t (@$v) {
                if (ref $t eq 'HASH' && $t->{stable_id}) {
                    $json_by_id{ $t->{stable_id} } = $t;
                }
            }
        } elsif (ref $v eq 'HASH') {
            for my $sub (values %$v) {
                if (ref $sub eq 'ARRAY') {
                    for my $t (@$sub) {
                        if (ref $t eq 'HASH' && $t->{stable_id}) {
                            $json_by_id{ $t->{stable_id} } = $t;
                        }
                    }
                }
            }
        }
    };
    $iter->($data);
}

# --- Recursive field-set walker ---
#
# Storable transcripts are BLESSED hashes (Bio::EnsEMBL::Transcript and the objects
# nested under it), so `ref $v` returns the class name, never 'HASH', and a walker
# keyed on `ref $v eq 'HASH'` descends no Storable object: every JSON field lands in
# `parity_gaps` as `only_in_json:*`, `identical_field_counts` and `conversion_rules`
# come out empty, and `missing_in_json` is vacuously zero, which reads as a clean audit
# that compared nothing. `reftype` sees through the blessing; the type recorded for a
# field is the underlying container type plus the class where one exists, so a blessed
# object on the Storable side against a plain object on the JSON side is reported as a
# CONVERSION (`Bio::EnsEMBL::Slice(HASH)->HASH`), not as identical and not as a gap.
sub _kind {
    my ($x) = @_;
    my $rt = reftype($x);
    return(defined $x ? "scalar" : "undef") unless defined $rt;
    my $cls = blessed($x);
    return $cls ? "$cls($rt)" : $rt;
}
# Ensembl objects hold back-references (a Translation to its Transcript, a Feature to
# its Slice), so the Storable graph is cyclic and an unguarded descent never ends.
# Hashes on the current descent path are not re-entered; a field that closes a cycle
# is recorded with its class and the `cycle` marker, and JSON, which cannot carry a
# cycle, reports whatever the converter put there. Shared objects that are not
# ancestors are walked at every path they appear under, so each path is compared.
# Keys are visited in sorted order so the manifest is independent of hash order.
sub field_keys {
    my ($v, $prefix, $out, $path_refs) = @_;
    $path_refs //= {};
    if ((reftype($v) // '') eq 'HASH') {
        my $addr = refaddr($v);
        return $out if $path_refs->{$addr};
        $path_refs->{$addr} = 1;
        for my $k (sort keys %$v) {
            my $path  = $prefix ? "$prefix.$k" : $k;
            my $child = $v->{$k};
            my $child_is_hash = (reftype($child) // '') eq 'HASH';
            if ($child_is_hash && $path_refs->{ refaddr($child) }) {
                $out->{$path} = _kind($child) =~ s/\)$/,cycle)/r;
                next;
            }
            $out->{$path} = _kind($child);
            field_keys($child, $path, $out, $path_refs) if $child_is_hash;
        }
        delete $path_refs->{$addr};
    }
    return $out;
}

# --- The converter's rules, restated so a missing field can be attributed -----------
#
# storable_to_json.pl whitelists fields, renames the underscore-prefixed Perl internals,
# flattens the UTR objects to sequence strings, replaces the exon coordinate mapper with
# a pair list, and copies prediction matrices only from a `predictions_data` field that
# the Storable cache does not carry. Each entry names the behaviour that accounts for
# a Storable path being absent under its own name on the JSON side.
my %RENAME = (
    '_gene_stable_id'                 => 'gene_stable_id',
    '_gene_symbol'                    => 'gene_symbol',
    '_gene_symbol_source'             => 'gene_symbol_source',
    '_gene_hgnc_id'                   => 'gene_hgnc_id',
    '_gene_phenotype'                 => 'gene_phenotype',
    '_gene_version'                   => 'gene_version',
    '_ccds'                           => 'ccds',
    '_refseq'                         => 'refseq',
    '_protein'                        => 'protein',
    '_swissprot'                      => 'swissprot',
    '_trembl'                         => 'trembl',
    '_uniparc'                        => 'uniparc',
    '_trans_exon_array'               => 'exons',
    '_variation_effect_feature_cache' => 'variation_effect_feature_cache',
);
my @DROPPED = (
    [qr/^slice(?:\.|$)/,                       'slice: the chromosome is the shard directory and no slice field is read'],
    [qr/^_gene(?:\.|$)/,                       'gene object: the gene fields are carried by the _gene_* renames'],
    [qr/^_vep_lazy_loaded$/,                   'runtime flag set by the Perl cache loader'],
    [qr/^translation\.(?:start|end)_exon(?:\.|$)/, 'translation exons: derivable from exons and the translation start and end'],
    [qr/^translation\.transcript$/,            'back-reference to the enclosing transcript'],
    [qr/^_variation_effect_feature_cache\.seq_edits$/, 'seq_edits are not converted'],
    [qr/^_variation_effect_feature_cache\.mapper\.exon_coord_mapper(?:\.|$)/, 'exon coordinate mapper: replaced by mapper.pairs'],
    [qr/^_variation_effect_feature_cache\.protein_function_predictions\.\w+\.(?:matrix|matrix_compressed|translation_md5)$/,
        'prediction matrices are not converted: the converter reads predictions_data, which only the native cache builder writes'],
    [qr/^_variation_effect_feature_cache\.(?:five|three)_prime_utr\./, 'UTR object flattened to its sequence string'],
);
my %DERIVED = (
    'variation_effect_feature_cache.mapper.pairs'      => 'pair list taken from mapper.pairs or computed from sorted_exons',
    'variation_effect_feature_cache.mapper.pair_count' => 'length of mapper.pairs',
);

sub mapped_path {
    my ($k) = @_;
    my @parts = split /\./, $k;
    $parts[0] = $RENAME{ $parts[0] } if exists $RENAME{ $parts[0] };
    return join ".", @parts;
}

my (%identical, %renamed, %converted, %dropped_undef, %dropped_rule, %unexplained, %derived, %only_json);
my %gap_examples;   # gap path -> up to three stable ids, so a gap can be looked up
sub note_gap { my ($key, $sid) = @_; my $l = $gap_examples{$key} //= []; push @$l, $sid if @$l < 3; }
my ($n_matched) = (0);
# The converter processes only the primary-assembly directories (1-22, X, Y, MT), so a
# sampled transcript from a patch, haplotype, LRG or unplaced-contig directory has no
# JSON counterpart by design; the count is reported per class so the unmatched sample
# is attributable rather than a silent shortfall.
my %unmatched = (canonical_contig => 0, non_canonical_contig => 0);

for my $entry (@sampled) {
    my ($chr, $st, $file) = @$entry;
    my $sid = $st->{stable_id} or next;
    unless (exists $json_by_id{$sid}) {
        $unmatched{ $chr =~ /^(?:\d+|X|Y|MT)$/ ? 'canonical_contig' : 'non_canonical_contig' }++;
        next;
    }
    $n_matched++;
    my $st_keys = field_keys($st, "", {});
    my $j_keys  = field_keys($json_by_id{$sid}, "", {});
    my %accounted;   # JSON paths explained by a Storable path or a derivation rule

    # Shortest paths first, so the children of a container the converter flattened to
    # a string are attributed to that flattening rather than reported one by one.
    my %flattened_prefix;
    for my $k (sort { length($a) <=> length($b) || $a cmp $b } keys %$st_keys) {
        my ($parent) = ($k =~ /^(.*)\.[^.]+$/);
        next if defined $parent && $flattened_prefix{$parent};
        my $m = mapped_path($k);
        if (exists $j_keys->{$m}) {
            $accounted{$m} = 1;
            my ($sk, $jk) = ($st_keys->{$k}, $j_keys->{$m});
            if ($sk eq $jk) {
                ($m eq $k ? \%identical : \%renamed)->{$k}++;
            } else {
                $converted{"$k:$sk->$jk"}++;
                $flattened_prefix{$k} = 1 if $jk !~ /HASH|ARRAY/;
            }
            next;
        }
        if ($st_keys->{$k} eq 'undef') { $dropped_undef{$k}++; next; }
        my ($rule) = grep { $k =~ $_->[0] } @DROPPED;
        if ($rule) { $dropped_rule{"$k: $rule->[1]"}++; next; }
        $unexplained{"missing_in_json:$k"}++;
        note_gap("missing_in_json:$k", $sid);
    }
    for my $m (keys %$j_keys) {
        next if $accounted{$m};
        if (exists $DERIVED{$m}) { $derived{"$m: $DERIVED{$m}"}++; next; }
        $only_json{"only_in_json:$m"}++;
        note_gap("only_in_json:$m", $sid);
    }
}

print encode_json({
    sampled_transcripts         => scalar(@sampled),
    sample_seed                 => $SEED,
    parity_gap_examples         => \%gap_examples,
    matched_in_json             => $n_matched,
    unmatched_in_json_by_contig_class => \%unmatched,
    storable_files_scanned      => scalar(@shuffled),
    json_transcripts_indexed    => scalar(keys %json_by_id),
    identical_field_counts      => \%identical,
    renamed_field_counts        => \%renamed,
    conversion_rules            => \%converted,
    dropped_undef_counts        => \%dropped_undef,
    dropped_by_converter        => \%dropped_rule,
    derived_on_json_side        => \%derived,
    parity_gaps                 => { %unexplained, %only_json },
});

