#!/usr/bin/env perl
# Builds a synthetic merged Ensembl VEP cache of one sequence, `contig`, 4,000,000 bases
# long, holding four single-exon non-coding transcripts on the forward strand, two of each
# set. Each set numbers its transcripts apart, as the two caches a merged cache is joined
# from do, so a `dbID` names one transcript of each set:
#
#   region file 1-1000000        set      dbID  span
#     XR_000001.1                RefSeq   7     100001-102000
#     ENSEXT00000000001          Ensembl  7     101001-103000
#     ENSEXT00000000002          Ensembl  9     500001-502000
#   region file 3000001-4000000
#     XR_000002.1                RefSeq   9     3200001-3202000
#
# The first pair overlaps and stands in one region file; the second pair lies 2.7 million
# bases apart in two region files. A region file lists its transcripts by start, as the
# merge of the two sets' files writes them.
# Usage: make_cache.pl <cache root>   (run with the perl of an Ensembl VEP installation)
use strict;
use warnings;
use Storable qw(nstore);
use IO::Compress::Gzip qw(gzip $GzipError);
use File::Path qw(make_path);

my ($root) = @ARGV;
die "usage: make_cache.pl <cache root>\n" unless defined $root;

my ($species, $assembly, $version) = ('two_set_example', 'EXAMPLE1', 116);
my ($name, $length, $region_size) = ('contig', 4_000_000, 1_000_000);

my $slice = bless({
    seq_region_name   => $name,
    seq_region_length => $length,
    start             => 1,
    end               => $length,
    strand            => 1,
    circular          => 0,
    coord_system      => bless({
        name => 'chromosome', version => $assembly, rank => 1, default => 1,
        sequence_level => 0, top_level => 0, dbID => 1, alias_to => undef,
    }, 'Bio::EnsEMBL::CoordSystem'),
}, 'Bio::EnsEMBL::Slice');

# One single-exon non-coding transcript as a cache dump holds it.
sub transcript {
    my (%t) = @_;
    my $bases = $t{end} - $t{start} + 1;
    my $exon = bless({
        stable_id => $t{exon}, start => $t{start}, end => $t{end}, strand => 1,
        phase => -1, end_phase => -1,
    }, 'Bio::EnsEMBL::Exon');
    my $pair = bless({
        from => bless({ id => 'cdna',   start => 1,         end => $bases }, 'Bio::EnsEMBL::Mapper::Unit'),
        to   => bless({ id => 'genome', start => $t{start}, end => $t{end} }, 'Bio::EnsEMBL::Mapper::Unit'),
        ori  => 1,
    }, 'Bio::EnsEMBL::Mapper::Pair');
    my $mapper = bless({
        cdna_coding_start => undef,
        cdna_coding_end   => undef,
        start_phase       => -1,
        exon_coord_mapper => bless({
            from => 'cdna', to => 'genomic', from_cs => undef, to_cs => undef,
            pair_count => 1, _is_sorted => 0,
            _pair_cdna    => { CDNA   => [$pair] },
            _pair_genomic => { GENOME => [$pair] },
        }, 'Bio::EnsEMBL::Mapper'),
    }, 'Bio::EnsEMBL::TranscriptMapper');
    return bless({
        stable_id => $t{id}, version => 1, dbID => $t{dbID}, biotype => 'lncRNA', source => $t{source},
        description => undef, attributes => [], translation => undef,
        start => $t{start}, end => $t{end}, strand => 1, slice => $slice,
        cdna_coding_start => undef, cdna_coding_end => undef,
        coding_region_start => undef, coding_region_end => undef,
        _gene => bless({ stable_id => $t{gene_object}, start => $t{start}, end => $t{end}, strand => 1 }, 'Bio::EnsEMBL::Gene'),
        _gene_stable_id => $t{gene}, _gene_version => 1, _gene_symbol => $t{symbol}, _gene_symbol_source => 'synthetic',
        _gene_phenotype => 0, _ccds => '-',
        _source_cache => $t{set},
        _trans_exon_array => [$exon],
        _vep_lazy_loaded  => 1,
        _variation_effect_feature_cache => {
            codon_table => 1, mapper => $mapper, introns => [], sorted_exons => [$exon],
            translateable_seq => '', five_prime_utr => undef, three_prime_utr => undef,
        },
    }, 'Bio::EnsEMBL::Transcript');
}

my @transcripts = (
    transcript(id => 'XR_000001.1', dbID => 7, set => 'RefSeq', source => 'Gnomon',
               gene => '900001', gene_object => 'gene1', symbol => 'exampleA', exon => 'id1',
               start => 100_001, end => 102_000),
    transcript(id => 'ENSEXT00000000001', dbID => 7, set => 'Ensembl', source => 'ensembl',
               gene => 'ENSEXG00000000001', gene_object => 'ENSEXG00000000001', symbol => 'exampleB', exon => 'ENSEXE00000000001',
               start => 101_001, end => 103_000),
    transcript(id => 'ENSEXT00000000002', dbID => 9, set => 'Ensembl', source => 'ensembl',
               gene => 'ENSEXG00000000002', gene_object => 'ENSEXG00000000002', symbol => 'exampleC', exon => 'ENSEXE00000000002',
               start => 500_001, end => 502_000),
    transcript(id => 'XR_000002.1', dbID => 9, set => 'RefSeq', source => 'Gnomon',
               gene => '900002', gene_object => 'gene2', symbol => 'exampleD', exon => 'id2',
               start => 3_200_001, end => 3_202_000),
);

my $dir = "$root/${species}_merged/${version}_$assembly";
make_path("$dir/$name");
open(my $info, '>', "$dir/info.txt") or die "$dir/info.txt: $!\n";
print $info "species\t$species\nassembly\t$assembly\nsift\t\npolyphen\t\nsource_assembly\t$assembly\nsource_genebuild\tsynthetic\nsource_refseq\tsynthetic\n";
close($info) or die "$dir/info.txt: $!\n";

# A transcript goes into the region file of each bin its span reaches, listed by start.
my %by_region;
for my $t (sort { $a->{start} <=> $b->{start} } @transcripts) {
    for my $bin (int(($t->{start} - 1) / $region_size) .. int(($t->{end} - 1) / $region_size)) {
        push @{ $by_region{$bin} }, $t;
    }
}
for my $bin (sort { $a <=> $b } keys %by_region) {
    my $file  = sprintf('%s/%s/%d-%d.gz', $dir, $name, $bin * $region_size + 1, ($bin + 1) * $region_size);
    my $plain = "$dir/$name/region.tmp";
    nstore({ $name => $by_region{$bin} }, $plain) or die "$plain: $!\n";
    gzip($plain => $file) or die "gzip: $GzipError\n";
    unlink($plain) or die "$plain: $!\n";
}
