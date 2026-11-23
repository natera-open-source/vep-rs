#!/usr/bin/env perl
# Print every genetic code of Bio::Tools::CodonTable, the module Ensembl VEP
# translates with, as tab-separated text: one line a codon in the module's own
# codon order, holding the table id, the codon, its amino acid, and 1 or 0 for
# each of is_start_codon and is_ter_codon. The ids 0 to 40 are asked in turn
# through `Bio::Tools::CodonTable->new(-id => <id>)`, the call Ensembl's
# `Transcript::translate` and BioPerl's `Bio::Seq::translate` make. An id holds
# a table when the module takes it without a warning and reads each of the 64
# codons as one amino acid letter or the stop sign. The `#` lines that open the
# output name the BioPerl release, the ids that hold no table and each id the
# module answers under another id.
#
# `crates/vep-core/tests/data/genetic_codes.tsv` is this output from the
# Ensembl VEP container image, where the script runs with this directory
# mounted read-only.
#
# Usage: dump_genetic_codes.pl > genetic_codes.tsv
use strict;
use warnings;
use File::Basename qw(dirname);
use Bio::Root::Version;
use Bio::Tools::CodonTable;

my ($first_id, $last_id) = (0, 40);

# The release string describes the module only when both load from one tree.
my $module_tree  = dirname(dirname($INC{'Bio/Tools/CodonTable.pm'}));
my $version_tree = dirname(dirname($INC{'Bio/Root/Version.pm'}));
die "ERROR: [dump_genetic_codes] Bio::Tools::CodonTable and Bio::Root::Version load from two trees\n"
    unless $module_tree eq $version_tree;

my $order = $Bio::Tools::CodonTable::TRCOL;
die "ERROR: [dump_genetic_codes] the module's codon order does not hold 64 codons\n"
    unless ref($order) eq 'HASH' && scalar(keys %$order) == 64;
my @codons = map { uc $order->{$_} } sort { $a <=> $b } keys %$order;

my (@rows, @without_a_table, @answered_as_another);
for my $id ($first_id .. $last_id) {
    my $warnings = 0;
    my $table;
    {
        local $SIG{__WARN__} = sub { $warnings++ };
        $table = Bio::Tools::CodonTable->new(-id => $id);
    }
    my @amino_acids = map { $table->translate($_) } @codons;
    if ($warnings || grep { !/^[A-Z*]$/ } @amino_acids) {
        push @without_a_table, $id;
        next;
    }
    push @answered_as_another, $id . '=' . $table->id if $table->id != $id;
    for my $i (0 .. $#codons) {
        push @rows, join("\t",
            $id, $codons[$i], $amino_acids[$i],
            $table->is_start_codon($codons[$i]) ? 1 : 0,
            $table->is_ter_codon($codons[$i])   ? 1 : 0);
    }
}

print "# module: Bio::Tools::CodonTable\n";
print "# bioperl: $Bio::Root::Version::VERSION\n";
print "# ids_asked: $first_id-$last_id\n";
print "# ids_without_a_table: @without_a_table\n";
print "# ids_answered_as_another: @answered_as_another\n";
print "# columns: table codon amino_acid start stop\n";
print "$_\n" for @rows;
