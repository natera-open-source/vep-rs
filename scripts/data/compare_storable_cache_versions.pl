#!/usr/bin/env perl
# Dump every transcript of one VEP Storable cache directory as one TSV line:
#   stable_id  contig  start  end  md5_full  md5_content
# md5_full is the MD5 of the transcript's canonical Data::Dumper text (sorted keys, no
# indentation). md5_content is the same after two normalisations that remove what a
# cache build changes without changing what VEP annotates: the keys `dbID`, `adaptor`,
# `created_date`, `modified_date` and `_gene_version` are deleted at every depth (database
# identifiers, dates and a gene-version label only some builds carry), and every scalar
# that is an integer is coerced to a number, because Storable preserves Perl's
# string-or-number flag and Data::Dumper quotes the string form ('12345' against 12345)
# while VEP's arithmetic reads both alike. Only the canonical contigs
# (1-22, X, Y, MT) are read, and only the `<start>-<end>.gz` transcript shards; `_var.gz`
# and `_reg.gz` shards are variation and regulatory data, not transcripts. A transcript
# that spans a shard boundary appears in two shards; it is written once, and the run
# aborts if its two copies differ, since that would make the dump ambiguous.
#
# Runs inside ensemblorg/ensembl-vep:release_115.2, which carries Storable,
# IO::Uncompress::Gunzip, Digest::MD5, Data::Dumper and the Bio::EnsEMBL classes the
# shards are blessed into. `compare_storable_cache_versions.py` joins two of these dumps.
#
# Usage: compare_storable_cache_versions.pl <cache_dir e.g. /cache/homo_sapiens/113_GRCh37> <out.tsv>
use strict;
use warnings;
use Storable qw(fd_retrieve);
use IO::Uncompress::Gunzip qw(gunzip);
use Digest::MD5 qw(md5_hex);
use Data::Dumper;

$Data::Dumper::Sortkeys = 1;
$Data::Dumper::Indent   = 0;
$Data::Dumper::Terse    = 1;

my ($cache_dir, $out) = @ARGV;
die "usage: $0 <cache_dir> <out.tsv>\n" unless defined $cache_dir && defined $out;
my %canonical = map { $_ => 1 } (1 .. 22, 'X', 'Y', 'MT');
my %STRIP = map { $_ => 1 } qw(dbID adaptor created_date modified_date _gene_version);

sub strip_copy {
    my ($node, $seen) = @_;
    unless (ref $node) {
        return $node unless defined $node;
        return $node + 0 if $node =~ /^-?\d+\z/;
        return $node;
    }
    return $seen->{$node} if exists $seen->{$node};
    my $r = ref $node;
    if ($r eq 'HASH' || (Scalar::Util::blessed($node) && Scalar::Util::reftype($node) eq 'HASH')) {
        my %copy;
        $seen->{$node} = \%copy;
        for my $k (keys %$node) {
            next if $STRIP{$k};
            $copy{$k} = strip_copy($node->{$k}, $seen);
        }
        return \%copy;
    }
    if ($r eq 'ARRAY' || (Scalar::Util::blessed($node) && Scalar::Util::reftype($node) eq 'ARRAY')) {
        my @copy;
        $seen->{$node} = \@copy;
        push @copy, strip_copy($_, $seen) for @$node;
        return \@copy;
    }
    if ($r eq 'SCALAR' || (Scalar::Util::blessed($node) && Scalar::Util::reftype($node) eq 'SCALAR')) {
        my $v = $$node;
        return \$v;
    }
    return $node;
}
require Scalar::Util;

open(my $fh_out, '>', $out) or die "cannot write $out: $!";
my %seen_id;
my ($n_files, $n_records, $n_dupes, $n_conflicts) = (0, 0, 0, 0);
for my $contig (sort keys %canonical) {
    my $dir = "$cache_dir/$contig";
    next unless -d $dir;
    for my $f (sort glob("$dir/*.gz")) {
        next unless $f =~ m{/(\d+)-(\d+)\.gz$};
        my $buf;
        gunzip($f => \$buf) or die "gunzip $f failed: $IO::Uncompress::Gunzip::GunzipError";
        open(my $fh, '<', \$buf) or die "in-memory open failed for $f";
        my $data = fd_retrieve($fh);
        close $fh;
        $n_files++;
        my @txs = ref($data) eq 'HASH' ? (map { @{ $_ || [] } } values %$data) : @{ $data || [] };
        for my $t (@txs) {
            next unless ref $t;
            my $id = $t->{stable_id} // $t->{_stable_id};
            next unless defined $id;
            $id =~ s/\.\d+$//;
            $n_records++;
            my $full    = md5_hex(Dumper($t));
            my $content = md5_hex(Dumper(strip_copy($t, {})));
            if (exists $seen_id{$id}) {
                $n_dupes++;
                if ($seen_id{$id} ne $content) {
                    $n_conflicts++;
                    warn "WARN: $id differs between two shards of $cache_dir\n";
                }
                next;
            }
            $seen_id{$id} = $content;
            my $start = $t->{start} // '';
            my $end   = $t->{end}   // '';
            print {$fh_out} join("\t", $id, $contig, $start, $end, $full, $content), "\n";
        }
    }
}
close $fh_out;
print STDERR sprintf("[compare_storable_cache_versions] %s: %d shard files, %d transcript records, %d distinct transcripts, %d shard-boundary duplicates, %d conflicting duplicates\n",
    $cache_dir, $n_files, $n_records, scalar(keys %seen_id), $n_dupes, $n_conflicts);
exit($n_conflicts ? 1 : 0);
