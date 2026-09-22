#!/usr/bin/env perl
# trace_perl_intermediates.pl: dump Perl VEP intermediate values for discordant variants.
#
# Runs inside the Perl VEP Docker container (ensemblorg/ensembl-vep:release_115.2).
# Requires the VEP cache to be mounted.
#
# Usage:
#   docker run --rm \
#     -v "${VEP_PERL_CACHE_DIR}:/opt/vep/.vep" \
#     -v $(pwd)/scripts:/scripts \
#     -v "${VEP_VCF_DIR}:/vcf" \
#     ensemblorg/ensembl-vep:release_115.2 \
#     perl /scripts/concordance/trace_perl_intermediates.pl \
#       --input /vcf/discordant.vcf \
#       --cache_dir /opt/vep/.vep \
#       --assembly GRCh37
#
# Output: TSV to stdout with columns:
#   location  allele  transcript_id  cds_start  cds_end  translation_start  translation_end
#   ref_codon  alt_codon  ref_peptide  alt_peptide  consequence_set

use strict;
use warnings;
use Getopt::Long;

# Add VEP module paths
use lib '/opt/vep/ensembl-vep';
use lib '/opt/vep/ensembl-vep/modules';

use Bio::EnsEMBL::VEP::Config;
use Bio::EnsEMBL::VEP::Runner;
use Bio::EnsEMBL::VEP::Parser::VCF;

my ($input, $cache_dir, $assembly, $fasta);
GetOptions(
    'input=s'     => \$input,
    'cache_dir=s' => \$cache_dir,
    'assembly=s'  => \$assembly,
    'fasta=s'     => \$fasta,
) or die "Usage: $0 --input VCF --cache_dir DIR --assembly GRCh37 [--fasta FA]\n";

die "--input required\n" unless $input;
$cache_dir //= '/opt/vep/.vep';
$assembly  //= 'GRCh37';

# Build VEP config
my %config_args = (
    input_file  => $input,
    dir_cache   => $cache_dir,
    assembly    => $assembly,
    offline     => 1,
    cache       => 1,
    format      => 'vcf',
    no_stats    => 1,
    buffer_size => 1000,
    quiet       => 1,
);
$config_args{fasta} = $fasta if $fasta;

my $runner = Bio::EnsEMBL::VEP::Runner->new(\%config_args);

# Initialize cache / adaptors
$runner->init();

# Print header
print join("\t", qw(
    location allele transcript_id
    cds_start cds_end translation_start translation_end
    ref_codon alt_codon ref_peptide alt_peptide
    consequence_set
)), "\n";

# Process each annotated input buffer
my $input_buffer = $runner->get_InputBuffer();
while (my $batch = $input_buffer->next()) {
    last unless $batch && scalar @$batch;

    for my $as (@{$runner->get_all_AnnotationSources}) {
        $as->annotate_InputBuffer($input_buffer);
    }

    $input_buffer->finish_annotation();
    next unless scalar @{$input_buffer->buffer};

    for my $vf (@{$input_buffer->buffer}) {
        my $loc = $vf->location_string();

        for my $tv (@{$vf->get_all_TranscriptVariations}) {
            my $tr_id = $tv->transcript->stable_id;

            for my $tva (@{$tv->get_all_alternate_TranscriptVariationAlleles}) {
                my $allele = $tva->variation_feature_seq // '-';

                # CDS coordinates
                my $cds_start = $tv->cds_start // '';
                my $cds_end   = $tv->cds_end   // '';
                my $tr_start  = $tv->translation_start // '';
                my $tr_end    = $tv->translation_end   // '';

                # Codon and peptide alleles
                my $ref_codon   = '';
                my $alt_codon   = '';
                my $ref_peptide = '';
                my $alt_peptide = '';

                eval {
                    my $codons = $tva->codon;
                    if ($codons && $codons =~ /(.+)\/(.+)/) {
                        $ref_codon = $1;
                        $alt_codon = $2;
                    }
                };

                eval {
                    my $peps = $tva->pep_allele_string;
                    if ($peps && $peps =~ /(.+)\/(.+)/) {
                        $ref_peptide = $1;
                        $alt_peptide = $2;
                    }
                };

                # Consequence terms
                my @cons = map { $_->SO_term } @{$tva->get_all_OverlapConsequences};
                my $cons_str = join(',', sort @cons);

                print join("\t",
                    $loc, $allele, $tr_id,
                    $cds_start, $cds_end, $tr_start, $tr_end,
                    $ref_codon, $alt_codon, $ref_peptide, $alt_peptide,
                    $cons_str,
                ), "\n";
            }
        }
    }
}

$runner->finish();
