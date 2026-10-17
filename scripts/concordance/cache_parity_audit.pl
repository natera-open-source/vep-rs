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
#                        the exon coordinate mapper's trees, ...)
#   unexplained          absent on the JSON side and matched by no rule: a parity gap
# JSON fields with no Storable counterpart are `derived` when a rule produces them and
# `only_in_json` otherwise. An array (exons, attributes, protein features, sequence
# edits) is compared as one container and then element by element under the path
# `<array>[]`, index-free, so every element's fields are classified.
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
#   json_cache_dir     = the converter's output directory (every *.json under it is read)
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
my @json_files;
find(
    sub { push @json_files, $File::Find::name if -f && /\.json$/ },
    $json_root
);
for my $jfile (sort @json_files) {
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
# Keys are visited in sorted order so the manifest is independent of hash order. The
# elements of an array are walked under one index-free path, `<array>[]`, so a field
# present on any element is classified once per transcript.
sub field_keys {
    my ($v, $prefix, $out, $path_refs) = @_;
    $path_refs //= {};
    my $rt = reftype($v) // '';
    if ($rt eq 'HASH') {
        my $addr = refaddr($v);
        return $out if $path_refs->{$addr};
        $path_refs->{$addr} = 1;
        for my $k (sort keys %$v) {
            _field_key($v->{$k}, ($prefix ? "$prefix.$k" : $k), $out, $path_refs);
        }
        delete $path_refs->{$addr};
    } elsif ($rt eq 'ARRAY') {
        for my $element (@$v) {
            _field_key($element, $prefix . '[]', $out, $path_refs);
        }
    }
    return $out;
}

sub _field_key {
    my ($child, $path, $out, $path_refs) = @_;
    my $child_rt = reftype($child) // '';
    if ($child_rt eq 'HASH' && $path_refs->{ refaddr($child) }) {
        $out->{$path} = _kind($child) =~ s/\)$/,cycle)/r;
        return;
    }
    $out->{$path} = _kind($child);
    field_keys($child, $path, $out, $path_refs) if $child_rt eq 'HASH' || $child_rt eq 'ARRAY';
}

# --- The converter's rules, restated so a missing field can be attributed -----------
#
# storable_to_json.pl whitelists fields, renames the underscore-prefixed Perl internals,
# flattens the UTR objects to sequence strings and each protein feature's analysis
# object to its `_display_label`, replaces the exon coordinate mapper with a pair list,
# and writes each prediction matrix's gzip blob base64-encoded under `predictions_data`.
# Each entry names the behaviour that accounts for a Storable path being absent under
# its own name on the JSON side.
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
    '_uniprot_isoform'                => 'uniprot_isoform',
    '_source_cache'                   => 'source_cache',
    '_gene_symbol_id'                 => 'gene_symbol_id',
    '_bam_edit_status'                => 'bam_edit_status',
    '_trans_exon_array'               => 'exons',
    '_variation_effect_feature_cache' => 'variation_effect_feature_cache',
);
my @DROPPED = (
    [qr/^slice(?:\.|$)/,                       'slice: the chromosome is the shard directory and no slice field is read'],
    [qr/^_gene(?:\.|$)/,                       'gene object: its stable_id and version are carried by the _gene_* renames; start, end and strand are read by nothing at run time'],
    [qr/^_vep_lazy_loaded$/,                   'runtime flag set by the Perl cache loader'],
    [qr/^translation\.(?:start|end)_exon(?:\.|$)/, 'translation exons: derivable from exons and the translation start and end'],
    [qr/^translation\.transcript$/,            'back-reference to the enclosing transcript'],
    [qr/^_variation_effect_feature_cache\.mapper\.exon_coord_mapper(?:\.|$)/, 'exon coordinate mapper: replaced by mapper.pairs'],
    [qr/^(?:_trans_exon_array|_variation_effect_feature_cache\.sorted_exons)\[\]\.strand$/, 'exon strand: the transcript strand'],
    [qr/^_variation_effect_feature_cache\.introns\[\]\.strand$/, 'intron strand: the transcript strand'],
    [qr/^_variation_effect_feature_cache\.introns\[\]\.(?:prev|next)(?:\.|$)/, 'intron exon references: the flanking exons are the exons array'],
    [qr/^_variation_effect_feature_cache\.introns\[\]\.slice(?:\.|$)/, 'intron slice: the chromosome is the shard directory and no slice field is read'],
    [qr/^_variation_effect_feature_cache\.seq_edits\[\]\.description$/, 'sequence edit description: read by nothing at run time'],
    [qr/^_variation_effect_feature_cache\.protein_function_predictions\.\w+\.matrix_compressed$/,
        'matrix_compressed: constant, every exported predictions_data is the gzip blob'],
    [qr/^_variation_effect_feature_cache\.(?:five|three)_prime_utr\./, 'UTR object flattened to its sequence string'],
);
my @DERIVED = (
    [qr/^variation_effect_feature_cache\.mapper\.pairs(?:\[\]|\.|$)/, 'pair list taken from mapper.pairs or computed from sorted_exons'],
    [qr/^variation_effect_feature_cache\.mapper\.pair_count$/,          'length of mapper.pairs'],
    [qr/^variation_effect_feature_cache\.seq_edits$/,                    'seq_edits: written as [] when the transcript has no translation'],
    [qr/^variation_effect_feature_cache\.seq_edits\[\]\.length_diff$/, 'seq_edits length_diff: computed from start, end and alt_seq as SeqEdit::length_diff does'],
);

sub mapped_path {
    my ($k) = @_;
    my @parts = split /\./, $k;
    my ($head, $suffix) = $parts[0] =~ /^(.*?)((?:\[\])?)$/;
    $parts[0] = $RENAME{$head} . $suffix if exists $RENAME{$head};
    my $path = join ".", @parts;
    $path =~ s/^(variation_effect_feature_cache\.protein_function_predictions\.\w+)\.matrix$/$1.predictions_data/;
    return $path;
}

my (%identical, %renamed, %converted, %dropped_undef, %dropped_rule, %unexplained, %derived, %only_json);
my %gap_examples;   # gap path -> up to three stable ids, so a gap can be looked up
sub note_gap { my ($key, $sid) = @_; my $l = $gap_examples{$key} //= []; push @$l, $sid if @$l < 3; }
my ($n_matched) = (0);
# The converter writes every chromosome directory, so a sampled transcript without a
# JSON counterpart is a gap; the count is reported per contig class (the assembled
# chromosomes against patches, haplotypes, LRG and unplaced contigs) so a conversion
# restricted to a chromosome list is attributable rather than a silent shortfall.
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
        my ($rule) = grep { $m =~ $_->[0] } @DERIVED;
        if ($rule) { $derived{"$m: $rule->[1]"}++; next; }
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

