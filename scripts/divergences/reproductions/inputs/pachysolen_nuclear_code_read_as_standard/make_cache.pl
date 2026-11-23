#!/usr/bin/env perl
# Builds a synthetic Ensembl VEP cache of one sequence, `contig`, 3,000 bases long, holding
# one forward protein-coding gene of six codons from base 1,001 whose record names genetic
# code 26, the nuclear code of Pachysolen tannophilus. The peptide is written as Ensembl's
# cache dump writes it, through `Bio::Tools::CodonTable->new(-id => 26)`: the module holds
# no table 26, warns, and translates under its last table, which gives the standard residues.
# Usage: make_cache.pl <cache root>   (run with the perl of an Ensembl VEP installation)
use strict;
use warnings;
use Storable qw(nstore);
use IO::Compress::Gzip qw(gzip $GzipError);
use File::Path qw(make_path);
use Bio::Tools::CodonTable;

my ($root) = @ARGV;
die "usage: make_cache.pl <cache root>\n" unless defined $root;

my ($species, $assembly, $version) = ('nuclear_code_example', 'EXAMPLE1', 116);
my ($name, $length) = ('contig', 3_000);
my ($start, $end) = (1_001, 1_018);
my $genetic_code = '26';

# The coding sequence: a start codon, an alanine codon, the codon genetic code 26 reads as
# alanine and the standard code as leucine, a leucine codon one base from it, a glutamate
# codon and a stop codon.
my @codons = qw(ATG GCT CTG CTA GAA TAA);
my $coding = join('', @codons);
my $bases  = length($coding);
die "the gene must span its coding bases\n" unless $bases == $end - $start + 1;

# The peptide as Ensembl's Transcript::translate gives it, without the final stop.
my $table   = Bio::Tools::CodonTable->new(-id => $genetic_code);
my $peptide = $table->translate($coding);
$peptide =~ s/\*$//;

my $slice = bless({
    seq_region_name   => $name,
    seq_region_length => $length,
    start             => 1,
    end               => $length,
    strand            => 1,
    coord_system      => bless({
        name => 'chromosome', version => $assembly, rank => 1, default => 1,
        sequence_level => 0, top_level => 0, dbID => 1, alias_to => undef,
    }, 'Bio::EnsEMBL::CoordSystem'),
}, 'Bio::EnsEMBL::Slice');
my $exon = bless({
    stable_id => 'CODEE1', start => $start, end => $end, strand => 1,
    phase => 0, end_phase => 0,
}, 'Bio::EnsEMBL::Exon');
my $pair = bless({
    from => bless({ id => 'cdna',   start => 1,               end => $bases },        'Bio::EnsEMBL::Mapper::Unit'),
    to   => bless({ id => 'genome', start => $start, end => $end }, 'Bio::EnsEMBL::Mapper::Unit'),
    ori  => 1,
}, 'Bio::EnsEMBL::Mapper::Pair');
my $mapper = bless({
    cdna_coding_start => 1,
    cdna_coding_end   => $bases,
    start_phase       => 0,
    exon_coord_mapper => bless({
        from => 'cdna', to => 'genomic', from_cs => undef, to_cs => undef,
        pair_count => 1, _is_sorted => 0,
        _pair_cdna    => { CDNA   => [$pair] },
        _pair_genomic => { GENOME => [$pair] },
    }, 'Bio::EnsEMBL::Mapper'),
}, 'Bio::EnsEMBL::TranscriptMapper');
my $transcript = bless({
    stable_id => 'CODET1', version => undef, dbID => 1, biotype => 'protein_coding', source => 'ensembl',
    description => undef, is_canonical => 1, attributes => [],
    start => $start, end => $end, strand => 1, slice => $slice,
    cdna_coding_start => 1, cdna_coding_end => $bases,
    coding_region_start => undef, coding_region_end => undef,
    _gene => bless({ stable_id => 'CODEG1', start => $start, end => $end, strand => 1 }, 'Bio::EnsEMBL::Gene'),
    _gene_stable_id => 'CODEG1', _gene_version => undef, _gene_symbol => 'codeA', _gene_symbol_source => 'synthetic',
    _gene_phenotype => 0, _protein => 'CODEP1',
    _ccds => '-', _refseq => '-', _swissprot => '-', _trembl => '-', _uniparc => '-', _uniprot_isoform => '-',
    _trans_exon_array => [$exon],
    _vep_lazy_loaded  => 1,
    _variation_effect_feature_cache => {
        codon_table => $genetic_code, mapper => $mapper, introns => [], sorted_exons => [$exon],
        translateable_seq => $coding, peptide => $peptide,
        five_prime_utr => undef, three_prime_utr => undef, protein_features => [], seq_edits => [],
    },
}, 'Bio::EnsEMBL::Transcript');
$transcript->{translation} = bless({
    stable_id => 'CODEP1', version => undef, dbID => 1, start => 1, end => $bases, seq => undef,
    start_exon => $exon, end_exon => $exon, transcript => $transcript,
}, 'Bio::EnsEMBL::Translation');

my $dir = "$root/$species/${version}_$assembly";
make_path("$dir/$name");
open(my $info, '>', "$dir/info.txt") or die "$dir/info.txt: $!\n";
print $info "species\t$species\nassembly\t$assembly\nsift\t\npolyphen\t\nsource_assembly\t$assembly\nsource_genebuild\tsynthetic\n";
close($info) or die "$dir/info.txt: $!\n";
my $plain = "$dir/$name/region.tmp";
nstore({ $name => [$transcript] }, $plain) or die "$plain: $!\n";
gzip($plain => "$dir/$name/1-1000000.gz") or die "gzip: $GzipError\n";
unlink($plain) or die "$plain: $!\n";
