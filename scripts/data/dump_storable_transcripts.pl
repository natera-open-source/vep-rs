#!/usr/bin/env perl
# Dump named transcripts of one VEP Storable cache as canonical Data::Dumper text, one
# block per transcript (`### <stable_id>` then the dump), for a field-level comparison of
# two cache builds. Prediction matrices are summarised as their byte length and
# translation md5 so the dump stays readable; everything else is verbatim. Runs inside
# ensemblorg/ensembl-vep:release_115.2 like compare_storable_cache_versions.pl, whose
# `differing_*.tsv` id list is the usual input; `classify_storable_transcript_differences.py`
# reads two of these dumps.
#
# Usage: dump_storable_transcripts.pl <cache_dir> <ids.txt (one unversioned id per line)> <out.txt>
use strict;
use warnings;
use Storable qw(fd_retrieve);
use IO::Uncompress::Gunzip qw(gunzip);
use Data::Dumper;

$Data::Dumper::Sortkeys = 1;
$Data::Dumper::Indent   = 1;
$Data::Dumper::Terse    = 1;

my ($dir, $idfile, $out) = @ARGV;
die "usage: $0 <cache_dir> <ids.txt> <out.txt>\n" unless defined $dir && defined $idfile && defined $out;
my %want;
open(my $ih, '<', $idfile) or die "cannot read $idfile: $!";
while (<$ih>) { chomp; next unless /^ENS/; s/\.\d+$//; $want{$_} = 1; }
close $ih;
open(my $oh, '>', $out) or die "cannot write $out: $!";
my %seen;
for my $contig (1 .. 22, 'X', 'Y', 'MT') {
    for my $f (sort glob("$dir/$contig/*.gz")) {
        next unless $f =~ m{/(\d+)-(\d+)\.gz$};
        my $buf;
        gunzip($f => \$buf) or die "gunzip $f failed";
        open(my $fh, '<', \$buf) or die;
        my $d = fd_retrieve($fh);
        my @t = ref($d) eq 'HASH' ? (map { @{ $_ || [] } } values %$d) : @{ $d || [] };
        for my $t (@t) {
            my $id = $t->{stable_id} // '';
            $id =~ s/\.\d+$//;
            next unless $want{$id} && !$seen{$id}++;
            my $c = {%$t};
            if (ref $c->{_variation_effect_feature_cache} eq 'HASH') {
                my %v = %{ $c->{_variation_effect_feature_cache} };
                if (ref $v{protein_function_predictions} eq 'HASH') {
                    my %pf;
                    for my $k (keys %{ $v{protein_function_predictions} }) {
                        my $m = $v{protein_function_predictions}{$k};
                        $pf{$k} = ref $m ? sprintf("matrix:%s md5:%s", defined $m->{matrix} ? length($m->{matrix}) : 'undef', $m->{translation_md5} // '') : 'none';
                    }
                    $v{protein_function_predictions} = \%pf;
                }
                $c->{_variation_effect_feature_cache} = \%v;
            }
            print {$oh} "### $id\n", Dumper($c), "\n";
        }
    }
}
close $oh;
print STDERR sprintf("[dump_storable_transcripts] %s: %d of %d requested transcripts dumped\n", $dir, scalar(keys %seen), scalar(keys %want));
