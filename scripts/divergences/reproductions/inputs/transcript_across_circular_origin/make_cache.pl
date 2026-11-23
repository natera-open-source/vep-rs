#!/usr/bin/env perl
# Builds a synthetic Ensembl VEP cache of one circular sequence, `ring`, 12,000 bases long,
# holding one forward protein-coding gene of 150 bases that lies across the origin: its
# first 90 bases are the sequence's last 90 (11,911 to 12,000) and its last 60 bases the
# sequence's first 60. The transcript is written as Ensembl's cache dump stores such a
# transcript: its start one sequence length below the database start with the database
# end unchanged (-89 to 60), a mapper pair that keeps the database coordinates (11,911
# to 60), and a slice flagged circular.
# Usage: make_cache.pl <cache root>   (run with the perl of an Ensembl VEP installation)
use strict;
use warnings;
use Storable qw(nstore);
use IO::Compress::Gzip qw(gzip $GzipError);
use File::Path qw(make_path);
use Bio::Tools::CodonTable;

my ($root) = @ARGV;
die "usage: make_cache.pl <cache root>\n" unless defined $root;

my ($species, $assembly, $version) = ('circular_example', 'RING1', 116);
my ($name, $length) = ('ring', 12_000);
my ($database_start, $database_end) = (11_911, 60);

# The coding sequence, ten codons a line: a start codon, 48 codons and a stop codon.
my @codons = qw(
    ATG TCC TGT CTG CAA ATC AAT AGG GCA GGC
    TCT CTA CAC ATT ACG AGA GCC GGT TTG CTC
    CAT CGG ACA AGC GCT GAG TTA TAC CTT CCG
    CGA ACC AGT GTG GAA TTC TAT TGG CCA CGC
    ACT AAG GTA GAC TTT TCG CCC CGT ATG TAA
);
my $coding = join('', @codons);
my $bases  = length($coding);
die "the gene must be 150 bases\n" unless $bases == 150 && $length - $database_start + 1 + $database_end == $bases;

# The peptide as Ensembl's Transcript::translate gives it: the bacterial code, no final stop.
my $table   = Bio::Tools::CodonTable->new(-id => 11);
my $peptide = $table->translate($coding);
$peptide =~ s/\*$//;

my $stored_start = $database_start - $length;
my $slice = bless({
    seq_region_name   => $name,
    seq_region_length => $length,
    start             => 1,
    end               => $length,
    strand            => 1,
    circular          => 1,
    coord_system      => bless({
        name => 'chromosome', version => $assembly, rank => 1, default => 1,
        sequence_level => 0, top_level => 0, dbID => 1, alias_to => undef,
    }, 'Bio::EnsEMBL::CoordSystem'),
}, 'Bio::EnsEMBL::Slice');
my $exon = bless({
    stable_id => 'RINGE1', start => $stored_start, end => $database_end, strand => 1,
    phase => 0, end_phase => 0,
}, 'Bio::EnsEMBL::Exon');
my $pair = bless({
    from => bless({ id => 'cdna',   start => 1,               end => $bases },        'Bio::EnsEMBL::Mapper::Unit'),
    to   => bless({ id => 'genome', start => $database_start, end => $database_end }, 'Bio::EnsEMBL::Mapper::Unit'),
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
    stable_id => 'RINGT1', version => undef, dbID => 1, biotype => 'protein_coding', source => 'ensembl',
    description => undef, is_canonical => 1, attributes => [],
    start => $stored_start, end => $database_end, strand => 1, slice => $slice,
    cdna_coding_start => 1, cdna_coding_end => $bases,
    coding_region_start => undef, coding_region_end => undef,
    _gene => bless({ stable_id => 'RINGG1', start => $stored_start, end => $database_end, strand => 1 }, 'Bio::EnsEMBL::Gene'),
    _gene_stable_id => 'RINGG1', _gene_version => undef, _gene_symbol => 'ringA', _gene_symbol_source => 'synthetic',
    _gene_phenotype => 0, _protein => 'RINGP1',
    _ccds => '-', _refseq => '-', _swissprot => '-', _trembl => '-', _uniparc => '-', _uniprot_isoform => '-',
    _trans_exon_array => [$exon],
    _vep_lazy_loaded  => 1,
    _variation_effect_feature_cache => {
        codon_table => 11, mapper => $mapper, introns => [], sorted_exons => [$exon],
        translateable_seq => $coding, peptide => $peptide,
        five_prime_utr => undef, three_prime_utr => undef, protein_features => [], seq_edits => [],
    },
}, 'Bio::EnsEMBL::Transcript');
$transcript->{translation} = bless({
    stable_id => 'RINGP1', version => undef, dbID => 1, start => 1, end => $bases, seq => undef,
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
