# An Ensembl VEP plugin that fails on purpose: its first run returns a value that
# is not a hash reference, its second dies, and every later run returns nothing,
# so a reference run with it loaded records how Ensembl VEP keeps the line and
# the fields of the other plugins when a plugin fails.
package Faulty;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

my $runs = 0;

sub feature_types {
  return ['Transcript'];
}

sub get_header_info {
  return {};
}

sub run {
  $runs++;
  return 1 if $runs == 1;
  die "Faulty on purpose\n" if $runs == 2;
  return {};
}

1;
