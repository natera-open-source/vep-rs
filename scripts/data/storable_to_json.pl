#!/usr/bin/perl
# Converts an Ensembl VEP Storable transcript cache into the JSON cache vep-rs
# reads: transcripts/<chr>/<start>-<end>.json for every chromosome directory,
# the variation files copied as they are into variations/<chr>/ (all_vars.gz with
# its index, or the <start>-<end>_var.gz region files of a cache dumped without
# the tabix conversion), regulatory/<chr>/<start>-<end>.json from *_reg.gz with
# the binding matrices they name stored once in regulatory/binding_matrices.json,
# info.json from the cache's info.txt, and chr_synonyms.txt copied beside them.
# Every field the cache dumper keeps on a transcript is written under its own or
# a renamed key, flattened (a UTR object
# to its sequence, a protein feature's analysis to its display label) or dropped
# under a rule the parity audit names (the DROPPED table of
# scripts/concordance/cache_parity_audit.pl).
use strict;
use warnings;
use Storable qw(fd_retrieve);
use JSON;
use File::Copy qw(copy);
use File::Path qw(make_path);
use MIME::Base64 qw(encode_base64);
use Scalar::Util qw(blessed reftype);

my $cache_dir = $ARGV[0] or die "Usage: $0 cache_dir output_dir [chr1,chr2,...]\n";
my $output_dir = $ARGV[1] or die "Usage: $0 cache_dir output_dir [chr1,chr2,...]\n";
my @chrs = $ARGV[2] ? split(/,/, $ARGV[2]) : ();

sub fail {
    my ($message) = @_;
    print STDERR "ERROR: [storable_to_json] $message\n";
    exit 1;
}

# Transcript keys copied under their own name.
my %TRANSCRIPT_FIELDS = map { $_ => 1 } qw(
    stable_id version dbID start end strand biotype source is_canonical
    cdna_coding_start cdna_coding_end
);
# The identifiers VEP sets on the transcript under underscore-prefixed keys
# (AnnotationSource/Database/Transcript.pm:272-284 and its prefetch_*_ids subs),
# written unprefixed.
my %TRANSCRIPT_RENAME = (
    '_gene_stable_id'    => 'gene_stable_id',
    '_gene_symbol'       => 'gene_symbol',
    '_gene_symbol_source'=> 'gene_symbol_source',
    '_gene_hgnc_id'      => 'gene_hgnc_id',
    '_gene_phenotype'    => 'gene_phenotype',
    '_ccds'              => 'ccds',
    '_refseq'            => 'refseq',
    '_gene_version'      => 'gene_version',
    '_protein'           => 'protein',
    '_swissprot'         => 'swissprot',
    '_trembl'            => 'trembl',
    '_uniparc'           => 'uniparc',
    '_uniprot_isoform'   => 'uniprot_isoform',
    '_source_cache'      => 'source_cache',
    '_gene_symbol_id'    => 'gene_symbol_id',
    '_bam_edit_status'   => 'bam_edit_status',
);

# Special compound fields handled by dedicated functions
# exons, translation, attributes, variation_effect_feature_cache

sub extract_exon {
    my ($exon) = @_;
    return undef unless ref $exon;
    my %out;
    for my $k (qw(stable_id start end phase end_phase rank)) {
        $out{$k} = $exon->{$k} if defined $exon->{$k};
    }
    return \%out;
}

sub extract_translation {
    my ($tr) = @_;
    return undef unless ref $tr;
    my %out;
    for my $k (qw(stable_id version dbID start end)) {
        $out{$k} = $tr->{$k} if defined $tr->{$k};
    }
    # seq may be a Bio::Seq blessed object: extract the string
    my $seq = $tr->{seq};
    if (defined $seq) {
        if (ref $seq) {
            # Bio::Seq object: extract sequence string
            $out{seq} = eval { $seq->seq() } // $seq->{seq} // $seq->{primary_seq}{seq} // undef;
        } else {
            $out{seq} = $seq;
        }
    }
    return \%out;
}

sub extract_attribute {
    my ($attr) = @_;
    return undef unless ref $attr;
    my %out;
    for my $k (qw(code name value)) {
        $out{$k} = $attr->{$k} if defined $attr->{$k};
    }
    return \%out;
}

sub extract_intron {
    my ($intron) = @_;
    return undef unless ref $intron;
    my %out;
    for my $k (qw(start end)) {
        $out{$k} = $intron->{$k} if defined $intron->{$k};
    }
    return \%out;
}

sub extract_protein_feature {
    my ($pf) = @_;
    return undef unless ref $pf;
    my %out;
    for my $k (qw(start end hseqname)) {
        $out{$k} = $pf->{$k} if defined $pf->{$k};
    }
    # The dumper keeps one key of the feature's analysis object, _display_label
    # (AnnotationSource/Database/Transcript.pm:474-477); it becomes the string
    # the DOMAINS column prints before the feature's own identifier.
    my $analysis = $pf->{analysis};
    if ((reftype($analysis) // '') eq 'HASH') {
        $out{analysis} = $analysis->{_display_label} if defined $analysis->{_display_label};
    } elsif (defined $analysis && !ref $analysis) {
        $out{analysis} = $analysis;
    }
    return \%out;
}

sub extract_mapper_pair {
    my ($pair) = @_;
    return undef unless ref $pair;
    my %out;
    for my $k (qw(from_start from_end to_start to_end ori)) {
        $out{$k} = $pair->{$k} if defined $pair->{$k};
    }
    return \%out;
}

# The pairs of a Bio::EnsEMBL::TranscriptMapper's exon_coord_mapper, a
# Bio::EnsEMBL::Mapper whose `_pair_<from>` hash holds, under the upper-cased
# coordinate system name, one Bio::EnsEMBL::Mapper::Pair per pair added
# (Mapper.pm add_map_coordinates): a `from` and a `to` Unit (id, start, end)
# and the orientation. TranscriptMapper::_load_mapper writes them with every
# insertion or deletion SeqEdit of the transcript folded in, so a BAM-edited
# transcript's cDNA coordinates come from here and never from its exons.
# Returned in cDNA order; empty when the object carries none.
sub mapper_object_pairs {
    my ($mapper) = @_;
    my $coord_mapper = $mapper->{exon_coord_mapper};
    return () unless (reftype($coord_mapper) // '') eq 'HASH' && defined $coord_mapper->{from};
    my $lists = $coord_mapper->{ '_pair_' . $coord_mapper->{from} };
    return () unless (reftype($lists) // '') eq 'HASH';
    my @pairs;
    for my $list (values %$lists) {
        next unless (reftype($list) // '') eq 'ARRAY';
        for my $pair (@$list) {
            next unless (reftype($pair) // '') eq 'HASH'
                && (reftype($pair->{from}) // '') eq 'HASH'
                && (reftype($pair->{to}) // '') eq 'HASH';
            push @pairs, {
                from_start => $pair->{from}{start} + 0,
                from_end   => $pair->{from}{end} + 0,
                to_start   => $pair->{to}{start} + 0,
                to_end     => $pair->{to}{end} + 0,
                ori        => ($pair->{ori} // 1) + 0,
            };
        }
    }
    return sort { $a->{from_start} <=> $b->{from_start} } @pairs;
}

sub extract_mapper {
    my ($mapper, $sorted_exons, $strand) = @_;
    return undef unless ref $mapper;
    my %out;
    for my $k (qw(start_phase cdna_coding_start cdna_coding_end)) {
        $out{$k} = $mapper->{$k} if defined $mapper->{$k};
    }

    my @object_pairs = mapper_object_pairs($mapper);
    # Try to get pairs from the mapper directly (GRCh37 format)
    if (ref $mapper->{pairs} eq 'ARRAY' && @{$mapper->{pairs}}) {
        $out{pairs} = [map { extract_mapper_pair($_) } grep { defined $_ && ref $_ } @{$mapper->{pairs}}];
        $out{pair_count} = scalar @{$out{pairs}};
    }
    elsif (@object_pairs) {
        $out{pairs} = \@object_pairs;
        $out{pair_count} = scalar @object_pairs;
    }
    # A mapper object without pairs: the exons give the unedited cDNA coordinates.
    elsif (ref $sorted_exons eq 'ARRAY' && @{$sorted_exons}) {
        my @pairs;
        my $cdna_pos = 1;
        # Sort exons by start position (forward strand) or reverse (reverse strand)
        my @exons = sort { ($strand && $strand < 0) ? ($b->{start} <=> $a->{start}) : ($a->{start} <=> $b->{start}) } @{$sorted_exons};
        for my $exon (@exons) {
            next unless ref $exon && $exon->{start} && $exon->{end};
            my $exon_len = $exon->{end} - $exon->{start} + 1;
            push @pairs, {
                from_start => $cdna_pos,
                from_end   => $cdna_pos + $exon_len - 1,
                to_start   => $exon->{start},
                to_end     => $exon->{end},
                ori        => ($strand && $strand < 0) ? -1 : 1,
            };
            $cdna_pos += $exon_len;
        }
        $out{pairs} = \@pairs;
        $out{pair_count} = scalar @pairs;
    }
    else {
        $out{pairs} = [];
        $out{pair_count} = 0;
    }

    return \%out;
}

my %matrices_written;

sub extract_prediction_matrix {
    my ($pm, $tool) = @_;
    return undef unless ref $pm;
    my %out;
    for my $k (qw(analysis sub_analysis peptide_length translation_md5)) {
        $out{$k} = $pm->{$k} if defined $pm->{$k};
    }
    # A Bio::EnsEMBL::Variation::ProteinFunctionPredictionMatrix holds its
    # predictions in `matrix`, the gzip blob the variation database stores
    # (a `VEP` header, then one little-endian short per position and amino
    # acid), and prediction_from_matrix gunzips it in place, so a matrix read
    # before the dump is the expanded bytes and is gzipped again here.
    if (defined $pm->{matrix}) {
        my $blob = $pm->{matrix};
        if (substr($blob, 0, 2) ne "\x1f\x8b") {
            require Compress::Zlib;
            $blob = Compress::Zlib::memGzip($blob) or die "gzip of a prediction matrix failed\n";
        }
        $out{predictions_data} = encode_base64($blob, '');
        $matrices_written{$tool}++;
    }
    return \%out;
}

sub extract_seq_edit {
    # A Bio::EnsEMBL::SeqEdit is a blessed hash whose start, end and alt_seq
    # were split out of the attribute value as strings; the class itself is not
    # loaded here, so its fields are read directly and length_diff is computed
    # as SeqEdit::length_diff does (0 when either coordinate is undefined).
    my ($se) = @_;
    return undef unless ref $se;
    my %out;
    for my $k (qw(code name)) {
        $out{$k} = $se->{$k} if defined $se->{$k};
    }
    for my $k (qw(start end)) {
        $out{$k} = $se->{$k} + 0 if defined $se->{$k};
    }
    $out{alt_seq} = defined $se->{alt_seq} ? "$se->{alt_seq}" : '';
    $out{length_diff} = (defined $se->{start} && defined $se->{end})
        ? length($out{alt_seq}) - ($out{end} - $out{start} + 1)
        : 0;
    return \%out;
}

sub extract_seq_string {
    # Extract a literal sequence string from a value that may be a plain string,
    # a blessed Bio::EnsEMBL::Slice (has ->seq method), a Bio::Seq, or a HASH
    # with a {seq} or {primary_seq}{seq} field. Returns undef when no string
    # can be recovered.
    #
    # Ensembl VEP's _variation_effect_feature_cache->{three_prime_utr} stores a
    # Bio::Seq (BaseTranscriptVariation.pm:1111 caches the transcript's
    # three_prime_utr()), a blessed hash whose ->seq() reads its primary_seq.
    # Without this extraction, JSON encoding emits null and vep-rs falls back to
    # the clamped-peptide path, which emits `stop_lost` where Perl emits a 3'UTR
    # term. Extracting the seq string keeps vep-rs offline+cache-only mode
    # self-contained (no FASTA dependency).
    my ($val) = @_;
    return undef unless defined $val;
    return $val unless ref $val;
    # Try ->seq method (Bio::EnsEMBL::Slice, Bio::Seq, etc.)
    my $seq = eval { $val->seq() };
    if (defined $seq && !ref $seq && $seq ne '') {
        return $seq;
    }
    # Try direct hash field (when blessed object's underlying hash holds the seq)
    if (ref $val eq 'HASH' || (blessed($val) && reftype($val) eq 'HASH')) {
        return $val->{seq} if defined $val->{seq} && !ref $val->{seq};
        if (ref $val->{primary_seq}) {
            my $ps_seq = eval { $val->{primary_seq}->seq() };
            return $ps_seq if defined $ps_seq && !ref $ps_seq && $ps_seq ne '';
            if (ref $val->{primary_seq} eq 'HASH' || (blessed($val->{primary_seq}) && reftype($val->{primary_seq}) eq 'HASH')) {
                return $val->{primary_seq}{seq} if defined $val->{primary_seq}{seq} && !ref $val->{primary_seq}{seq};
            }
        }
    }
    return undef;
}

sub extract_vefc {
    my ($vefc, $transcript_strand) = @_;
    return undef unless ref $vefc;
    # Use the underscore-prefixed name from Storable
    $vefc = $vefc->{'_variation_effect_feature_cache'} || $vefc->{'variation_effect_feature_cache'} || $vefc
        if exists $vefc->{'_variation_effect_feature_cache'} || exists $vefc->{'variation_effect_feature_cache'};

    my %out;
    $out{codon_table} = $vefc->{codon_table} if defined $vefc->{codon_table};
    # five/three_prime_utr are blessed Bio::Seq objects in Storable;
    # extract the underlying sequence string so vep-rs can use cache-only mode
    # without --fasta. See extract_seq_string above for rationale.
    if (defined $vefc->{five_prime_utr}) {
        my $s = extract_seq_string($vefc->{five_prime_utr});
        $out{five_prime_utr} = $s if defined $s;
    }
    if (defined $vefc->{three_prime_utr}) {
        my $s = extract_seq_string($vefc->{three_prime_utr});
        $out{three_prime_utr} = $s if defined $s;
    }
    $out{translateable_seq} = $vefc->{translateable_seq} if defined $vefc->{translateable_seq};
    $out{peptide} = $vefc->{peptide} if defined $vefc->{peptide};
    # Only a BAM-edited cache stores the spliced cDNA (AnnotationType/Transcript.pm:617).
    $out{spliced_seq} = $vefc->{spliced_seq} if defined $vefc->{spliced_seq};
    # The translation's SeqEdits (initial_met, _selenocysteine, amino_acid_sub,
    # _stop_codon_rt), which VEP applies to the reference allele's peptide at
    # runtime (TranscriptVariationAllele::peptide). The cached peptide has them
    # applied already, and Transcript::translate writes M over any start codon
    # of the table with or without an edit, so the peptide alone cannot say
    # which residues are edited. Always an array: a translation with none, or
    # no translation, gets [] as BaseTranscriptVariation::_seq_edits does.
    $out{seq_edits} = ref $vefc->{seq_edits} eq 'ARRAY'
        ? [map { extract_seq_edit($_) } grep { defined $_ && ref $_ } @{$vefc->{seq_edits}}]
        : [];

    # introns
    if (ref $vefc->{introns} eq 'ARRAY') {
        $out{introns} = [map { extract_intron($_) } grep { defined $_ && ref $_ } @{$vefc->{introns}}];
    }

    # sorted_exons
    if (ref $vefc->{sorted_exons} eq 'ARRAY') {
        $out{sorted_exons} = [map { extract_exon($_) } grep { defined $_ && ref $_ } @{$vefc->{sorted_exons}}];
    }

    # mapper: pass sorted_exons and strand to compute pairs if needed
    if (ref $vefc->{mapper}) {
        $out{mapper} = extract_mapper($vefc->{mapper}, $vefc->{sorted_exons}, $transcript_strand);
    }

    # protein_features
    if (ref $vefc->{protein_features} eq 'ARRAY') {
        $out{protein_features} = [map { extract_protein_feature($_) } grep { defined $_ && ref $_ } @{$vefc->{protein_features}}];
    }

    # protein_function_predictions
    my $pfp = $vefc->{protein_function_predictions};
    if (ref $pfp eq 'HASH' || (ref $pfp && eval { keys %{$pfp}; 1 })) {
        my %pfp_out;
        for my $key (qw(sift polyphen_humvar polyphen_humdiv)) {
            if (ref $pfp->{$key}) {
                my ($tool) = split /_/, $key;
                $pfp_out{$key} = extract_prediction_matrix($pfp->{$key}, $tool);
            }
        }
        $out{protein_function_predictions} = \%pfp_out if keys %pfp_out;
    }

    return \%out;
}

sub extract_transcript {
    my ($t) = @_;
    return undef unless ref $t;
    my %out;

    for my $k (keys %TRANSCRIPT_FIELDS) {
        $out{$k} = $t->{$k} if defined $t->{$k};
    }
    for my $pk (keys %TRANSCRIPT_RENAME) {
        $out{ $TRANSCRIPT_RENAME{$pk} } = $t->{$pk} if defined $t->{$pk};
    }

    # exons: from _trans_exon_array or exons
    my $exons = $t->{exons} || $t->{'_trans_exon_array'};
    if (ref $exons eq 'ARRAY') {
        $out{exons} = [map { extract_exon($_) } grep { defined $_ && ref $_ } @{$exons}];
    }

    # translation
    $out{translation} = extract_translation($t->{translation}) if ref $t->{translation};

    # attributes
    if (ref $t->{attributes} eq 'ARRAY') {
        $out{attributes} = [map { extract_attribute($_) } grep { defined $_ && ref $_ } @{$t->{attributes}}];
    }

    # variation_effect_feature_cache
    my $vefc = $t->{'variation_effect_feature_cache'} || $t->{'_variation_effect_feature_cache'};
    if (ref $vefc) {
        $out{variation_effect_feature_cache} = extract_vefc($vefc, $t->{strand});
    }

    return \%out;
}
my $json = JSON->new->utf8->canonical->allow_blessed->convert_blessed;

# The shape of the regulatory shards and of regulatory/binding_matrices.json;
# info.json carries it as `regulatory_schema`, and the annotator refuses a
# regulatory tree whose info.json declares none.
my $REGULATORY_SCHEMA = 1;

# One code per activity value of the funcgen regulatory_activity table
# (sql/table.sql: INACTIVE, REPRESSED, POISED, ACTIVE, NA). A feature's
# `activity` string holds one code per name of info.json's cell_types list, in
# list order, and `-` where the feature's map lacks the name.
my %ACTIVITY_CODE = (ACTIVE => 'A', INACTIVE => 'I', POISED => 'P', REPRESSED => 'R', NA => 'N');

# The binding matrices the shards name, each stored once under its stable id as
# its canonical JSON, with the shard that first carried it: the Storable cache
# repeats the whole matrix on every motif, and a repeat that differs is corrupt.
my %binding_matrices;
my %binding_matrix_source;

sub read_regulatory_shard {
    my ($in) = @_;
    open my $fh, '-|', 'gzip', '-dc', $in or die "gzip pipe failed: $!";
    my $data = fd_retrieve($fh);
    close($fh) or die "gzip -dc failed (wait status $?)\n";
    die "invalid regulatory Storable image\n" unless (reftype($data) // '') eq 'HASH';
    for my $chr (keys %$data) {
        die "invalid regulatory chromosome entry\n" unless (reftype($data->{$chr}) // '') eq 'HASH';
        for my $type (keys %{$data->{$chr}}) {
            die "invalid regulatory feature list\n" unless (reftype($data->{$chr}{$type}) // '') eq 'ARRAY';
        }
    }
    return $data;
}

# The cell_types list of info.json: the sorted union of every regulatory
# feature's activity-map keys (Database/RegFeat.pm:196-203, epigenome short
# names with spaces as `_`) over the chromosomes being converted. The shards are
# read once for the list before any is written, because every feature's
# activity string is aligned to it; the cache's own info.txt line cannot be
# split back into names, since some names contain commas.
sub collect_cell_types {
    my (@chr_dirs) = @_;
    my %names;
    for my $chr (@chr_dirs) {
        my $chr_path = "$cache_dir/$chr";
        opendir(my $ch, $chr_path) or fail("cannot open $chr_path: $!");
        my @files = sort grep { /^\d+-\d+_reg\.gz$/ } readdir($ch);
        closedir($ch);
        for my $file (@files) {
            my $data = eval { read_regulatory_shard("$chr_path/$file") };
            fail("$chr_path/$file: $@") if $@;
            for my $c (keys %$data) {
                for my $rf (@{ $data->{$c}{RegulatoryFeature} // [] }) {
                    my $map = (reftype($rf) // '') eq 'HASH' ? $rf->{cell_types} : undef;
                    next unless (reftype($map) // '') eq 'HASH';
                    $names{$_} = 1 for keys %$map;
                }
            }
        }
    }
    return [sort keys %names];
}

sub activity_string {
    my ($feature, $cell_types) = @_;
    my $map = $feature->{cell_types};
    $map = {} unless (reftype($map) // '') eq 'HASH';
    my $codes = '';
    for my $name (@$cell_types) {
        my $value = $map->{$name};
        if (!defined $value) {
            $codes .= '-';
            next;
        }
        my $code = $ACTIVITY_CODE{$value};
        die "RegulatoryFeature $feature->{stable_id} has activity '$value' in $name; expected one of "
            . join(', ', sort keys %ACTIVITY_CODE) . "\n" unless defined $code;
        $codes .= $code;
    }
    return $codes;
}

# A motif's binding matrix as the table stores it (funcgen BindingMatrix.pm:
# `elements` by position then nucleotide, `length` the number of positions,
# the complexes' display names as TRANSCRIPTION_FACTORS prints them). Returns
# the stable id and the entry.
sub extract_binding_matrix {
    my ($matrix, $motif_id) = @_;
    die "MotifFeature $motif_id lacks binding_matrix\n" unless defined $matrix;
    die "invalid binding matrix\n" unless (reftype($matrix) // '') eq 'HASH';
    die "binding matrix lacks stable_id\n" unless defined $matrix->{stable_id};
    my $id = "$matrix->{stable_id}";
    my %entry;
    $entry{name} = "$matrix->{name}" if defined $matrix->{name};
    $entry{source} = "$matrix->{source}" if defined $matrix->{source};
    $entry{threshold} = $matrix->{threshold} + 0 if defined $matrix->{threshold};
    # The motif score arithmetic reads the elements as counts, as funcgen's
    # converter does (BindingMatrix/Converter.pm:108-114 throws on another unit).
    die "binding matrix $id has unit '$matrix->{unit}'; only Frequencies converts\n"
        if defined $matrix->{unit} && $matrix->{unit} ne 'Frequencies';
    my $complexes = $matrix->{associated_transcription_factor_complexes} // [];
    die "invalid transcription factor complexes\n" unless (reftype($complexes) // '') eq 'ARRAY';
    $entry{transcription_factors} = [
        map { "$_->{display_name}" }
        grep { (reftype($_) // '') eq 'HASH' && defined $_->{display_name} }
        @$complexes
    ];
    my $elements = $matrix->{elements};
    my $length = $matrix->{length};
    if (defined $elements) {
        die "invalid binding matrix elements\n" unless (reftype($elements) // '') eq 'HASH';
        my $positions = scalar keys %$elements;
        $length //= $positions;
        die "binding matrix $id declares length $length but has $positions positions\n"
            unless $length == $positions;
        my @rows;
        for my $position (1 .. $length) {
            my $row = $elements->{$position};
            die "binding matrix $id lacks position $position\n" unless (reftype($row) // '') eq 'HASH';
            push @rows, [map {
                die "binding matrix $id lacks $_ at position $position\n" unless defined $row->{$_};
                $row->{$_} + 0
            } qw(A C G T)];
        }
        $entry{frequencies} = \@rows;
    } else {
        die "binding matrix $id has neither elements nor a length\n" unless defined $length;
        $entry{frequencies} = [];
    }
    $entry{length} = $length + 0;
    return ($id, \%entry);
}

sub record_binding_matrix {
    my ($id, $entry, $file) = @_;
    my $encoded = $json->encode($entry);
    if (exists $binding_matrices{$id}) {
        die "binding matrix $id differs from its copy in $binding_matrix_source{$id}\n"
            unless $binding_matrices{$id} eq $encoded;
        return;
    }
    $binding_matrices{$id} = $encoded;
    $binding_matrix_source{$id} = $file;
}

# RegFeat.pm stores chromosome => feature class => array of blessed hashes.
# OutputFactory.pm uses feature_type (or its so_name) for regulatory BIOTYPE,
# and complex display_name, not component gene names, for TRANSCRIPTION_FACTORS.
# A motif's stranded sequence is the one the dumper stored (Regulation.pm:137).
sub extract_regulatory_feature {
    my ($feature, $type, $cell_types, $file) = @_;
    die "invalid $type object\n" unless (reftype($feature) // '') eq 'HASH';
    my %out = (feature_type => $type);
    for my $key (qw(stable_id start end dbID)) {
        die "$type lacks $key\n" unless defined $feature->{$key};
    }
    $out{stable_id} = "$feature->{stable_id}";
    $out{start} = $feature->{start} + 0;
    $out{end} = $feature->{end} + 0;
    $out{dbid} = $feature->{dbID} + 0;
    if ($type eq 'RegulatoryFeature') {
        my $biotype = $feature->{feature_type};
        $biotype = $biotype->{so_name} if (reftype($biotype) // '') eq 'HASH';
        $out{biotype} = "$biotype" if defined $biotype && !ref $biotype;
        $out{activity} = activity_string($feature, $cell_types);
    } elsif ($type eq 'MotifFeature') {
        die "MotifFeature $out{stable_id} lacks strand\n" unless defined $feature->{strand};
        $out{strand} = $feature->{strand} + 0;
        die "MotifFeature $out{stable_id} lacks score\n" unless defined $feature->{score};
        $out{score} = $feature->{score} + 0;
        my $vefc = $feature->{_variation_effect_feature_cache};
        my $seq = (reftype($vefc) // '') eq 'HASH' ? $vefc->{seq} : undef;
        die "MotifFeature $out{stable_id} lacks its sequence\n" unless defined $seq && !ref $seq;
        my $span = $out{end} - $out{start} + 1;
        die "MotifFeature $out{stable_id} has a sequence of " . length($seq) . " bases on a span of $span\n"
            unless length($seq) == $span;
        $out{seq} = "$seq";
        # Every motif feature carries a binding matrix: the funcgen schema requires
        # one and the GFF route blesses a stub, so a motif without one is corrupt.
        my ($matrix_id, $entry) = extract_binding_matrix($feature->{binding_matrix}, $out{stable_id});
        record_binding_matrix($matrix_id, $entry, $file);
        $out{matrix} = $matrix_id;
    } else {
        die "unknown regulatory feature type $type\n";
    }
    return \%out;
}

sub convert_regulatory_files {
    my ($chr_path, $out_chr, $entries, $cell_types) = @_;
    my @files = sort grep { /^\d+-\d+_reg\.gz$/ } @$entries;
    return 0 unless @files;
    make_path($out_chr);
    for my $file (@files) {
        my $in = "$chr_path/$file";
        (my $base = $file) =~ s/_reg\.gz$/.json/;
        my $out = "$out_chr/$base";
        eval {
            my $data = read_regulatory_shard($in);
            my @features;
            for my $chr (sort keys %$data) {
                # RegulatoryFeature before MotifFeature, as Cache/RegFeat.pm:212-213 lists them.
                for my $type (reverse sort keys %{$data->{$chr}}) {
                    push @features,
                        map { extract_regulatory_feature($_, $type, $cell_types, $in) }
                        @{$data->{$chr}{$type}};
                }
            }
            open my $ofh, '>', $out or die "write failed: $!";
            print $ofh $json->encode(\@features) or die "write failed: $!";
            close($ofh) or die "write failed: $!";
        };
        if ($@) {
            my $error = $@;
            unlink $out;
            fail("$in: $error");
        }
    }
    return scalar @files;
}

sub write_binding_matrices {
    my ($out_path) = @_;
    my $table = '{' . join(',', map { $json->encode($_) . ':' . $binding_matrices{$_} } sort keys %binding_matrices) . "}\n";
    open(my $ofh, '>', $out_path) or fail("cannot write $out_path: $!");
    print $ofh $table;
    unless (close($ofh)) {
        my $error = $!;
        unlink $out_path;
        fail("write of $out_path failed: $error");
    }
    print "Wrote $out_path (" . scalar(keys %binding_matrices) . " binding matrices)\n";
}


my $info = read_info_txt("$cache_dir/info.txt");
fail("$cache_dir uses the sereal serialiser; only a Storable cache converts")
    if ($info->{serialiser_type} // '') eq 'sereal';

opendir(my $dh, $cache_dir) or fail("cannot open $cache_dir: $!");
# Every directory is a chromosome, as VEP lists them (CacheDir.pm:431): the
# assembled chromosomes, patches, alternate loci, scaffolds and LRG regions.
my @chr_dirs = sort grep { -d "$cache_dir/$_" && $_ !~ /^\./ } readdir($dh);
closedir($dh);
if (@chrs) {
    my %present = map { $_ => 1 } @chr_dirs;
    my @absent = grep { !$present{$_} } @chrs;
    fail("no chromosome directory named @absent under $cache_dir") if @absent;
    my %wanted = map { $_ => 1 } @chrs;
    @chr_dirs = grep { $wanted{$_} } @chr_dirs;
}

print "Processing " . scalar(@chr_dirs) . " chromosomes\n";
# A run that stops early leaves none of these files, whatever an earlier run left there.
unlink "$output_dir/info.json", "$output_dir/chr_synonyms.txt", "$output_dir/regulatory/binding_matrices.json";
my $cell_types = collect_cell_types(@chr_dirs);

my $total = 0;
my $ok = 0;
my $transcripts_total = 0;
my $variation_files_total = 0;
my $regulatory_files_total = 0;

foreach my $chr (@chr_dirs) {
    my $chr_path = "$cache_dir/$chr";
    my $out_chr = "$output_dir/transcripts/$chr";

    opendir(my $ch, $chr_path) or fail("cannot open $chr_path: $!");
    my @entries = readdir($ch);
    closedir($ch);
    my @files = sort grep { /^\d+-\d+\.gz$/ } @entries;
    fail("$chr_path holds sereal shards; only a Storable cache converts")
        if !@files && grep { /^\d+-\d+\.sereal$/ } @entries;
    make_path($out_chr);

    $total += scalar @files;
    my $chr_tx = 0;

    foreach my $f (@files) {
        my $in = "$chr_path/$f";
        (my $base = $f) =~ s/\.gz$/.json/;
        my $out = "$out_chr/$base";

        eval {
            open my $fh, '-|', 'gzip', '-dc', $in or die "gzip pipe failed: $!";
            my $data = fd_retrieve($fh);
            # close on a piped open is false when gzip exited non-zero, and
            # fd_retrieve returns undef, not an exception, when the image ends early.
            close($fh) or die "gzip -dc failed (wait status $?)\n";
            die "fd_retrieve returned no data (truncated Storable image)\n" unless defined $data;

            my @transcripts;
            if (ref($data) && (reftype($data) eq 'HASH' || blessed($data))) {
                for my $key (keys %{$data}) {
                    my $arr = $data->{$key};
                    next unless ref($arr) && reftype($arr) eq 'ARRAY';
                    for my $t (@{$arr}) {
                        if (defined $t && ref $t) {
                            my $clean = extract_transcript($t);
                            push @transcripts, $clean if $clean;
                        }
                    }
                }
            } elsif (ref($data) && reftype($data) eq 'ARRAY') {
                for my $t (@{$data}) {
                    if (defined $t && ref $t) {
                        my $clean = extract_transcript($t);
                        push @transcripts, $clean if $clean;
                    }
                }
            }

            $chr_tx += scalar @transcripts;
            open(my $ofh, ">", $out) or die "write failed: $!";
            print $ofh $json->encode(\@transcripts);
            close($ofh) or die "write failed: $!";
            $ok++;
        };
        if ($@) {
            unlink $out;
            warn "WARN: $in: $@\n";
        }
    }
    $transcripts_total += $chr_tx;
    my $variation_files = copy_variation_files($chr_path, "$output_dir/variations/$chr", \@entries);
    $variation_files_total += $variation_files;
    my $regulatory_files = convert_regulatory_files($chr_path, "$output_dir/regulatory/$chr", \@entries, $cell_types);
    $regulatory_files_total += $regulatory_files;
    print "  chr$chr: " . scalar(@files) . " files, $chr_tx transcripts, $variation_files variation files, $regulatory_files regulatory files\n";
}
fail("no transcript shards (<start>-<end>.gz) under $cache_dir") unless $total;
my $failed_shards = $total - $ok;
# A failed shard stops the run here, before chr_synonyms.txt and info.json are
# written, so a partial output lacks both.
fail("$failed_shards of $total shards failed to convert (see the WARN lines above)") if $failed_shards;

# VEP reads the cache's chr_synonyms.txt when --synonyms is not given (CacheDir.pm:254).
my $synonyms_in  = "$cache_dir/chr_synonyms.txt";
my $synonyms_out = "$output_dir/chr_synonyms.txt";
if (-e $synonyms_in) {
    copy($synonyms_in, $synonyms_out) or fail("copy of chr_synonyms.txt failed: $!");
    print "Wrote $synonyms_out\n";
} else {
    warn "WARN: no $synonyms_in; the cache carries no chromosome synonyms\n";
}

$info->{contig_order} = [sort by_contig_order @chr_dirs];
# The dumper's `sift b` and `polyphen b` lines say the cache was dumped with
# predictions. VEP reads a tool's availability from its `source_<tool>` line
# alone (Cache/Transcript.pm:155 over the version data CacheDir.pm:468-470
# collects), so a marker whose version line info.txt lacks is dropped, as is
# one whose matrices this conversion did not write.
for my $tool (qw(sift polyphen)) {
    delete $info->{$tool} unless $matrices_written{$tool} && $info->{source_versions}{$tool};
}
# An input marker alone must not advertise a tree that was never converted; the
# cell-type list and the schema describe the written tree alone.
delete $info->{$_} for qw(regulatory cell_types regulatory_schema);
if ($regulatory_files_total) {
    write_binding_matrices("$output_dir/regulatory/binding_matrices.json");
    $info->{regulatory} = JSON::true;
    $info->{cell_types} = $cell_types;
    $info->{regulatory_schema} = $REGULATORY_SCHEMA;
}
write_info_json($info, "$output_dir/info.json");

print "\nDone. Converted: $ok / $total files, $transcripts_total total transcripts, $variation_files_total variation files copied, $regulatory_files_total regulatory files\n";

# The variation files are read by vep-rs in the dumper's own format
# (Pipeline/DumpVEP/Dumper/Variation.pm): all_vars.gz with the .csi or .tbi
# tabix built on it, else every <start>-<end>_var.gz region file. Copied byte
# for byte, never re-encoded. An all_vars.gz without its index stops the run:
# the reader queries the index and never scans the file.
sub copy_variation_files {
    my ($chr_path, $out_chr, $entries) = @_;
    my %present = map { $_ => 1 } @{$entries};
    my @names_to_copy;
    if ($present{'all_vars.gz'}) {
        my ($index) = grep { $present{$_} } ('all_vars.gz.csi', 'all_vars.gz.tbi');
        fail("$chr_path/all_vars.gz has no all_vars.gz.csi or all_vars.gz.tbi beside it") unless $index;
        @names_to_copy = ('all_vars.gz', $index);
    } else {
        @names_to_copy = sort grep { /^\d+-\d+_var\.gz$/ } @{$entries};
    }
    return 0 unless @names_to_copy;
    make_path($out_chr);
    for my $name (@names_to_copy) {
        copy("$chr_path/$name", "$out_chr/$name") or fail("copy of $chr_path/$name failed: $!");
    }
    return scalar @names_to_copy;
}

# info.txt is translated the way VEP's own reader does (CacheDir.pm
# read_cache_info_file): every `source_<key>` line becomes a source version,
# `variation_cols` becomes a list, `-` means unset, and the `cell_types` line is
# dropped (its comma-joined names cannot be split back; collect_cell_types
# rebuilds the list from the features). Regulatory availability is set only
# after its shards have been written.
sub read_info_txt {
    my ($info_txt) = @_;
    my %info = (source_versions => {});
    if (open my $ifh, "<", $info_txt) {
        while (my $line = <$ifh>) {
            chomp $line;
            next if $line =~ /^#/ || $line !~ /\t/;
            my ($key, $value) = split /\t/, $line, 2;
            next if $key eq 'cell_types';
            if ($key =~ s/^source_//) {
                $info{source_versions}{$key} = $value;
            } elsif ($key eq 'variation_cols') {
                $info{$key} = [grep { length } split /,/, $value];
            } elsif ($key eq 'regulatory') {
                next;
            } elsif ($value ne '-') {
                $info{$key} = $value;
            }
        }
        close $ifh;
    } else {
        warn "WARN: no $info_txt; info.json carries no source versions\n";
    }
    # The cache version is the trailing path component of a VEP cache directory
    # (<species>[_<assembly>]/<version>), never a line of info.txt.
    if ($cache_dir =~ m{/(\d+)(?:_[^/]+)?/?$}) {
        $info{cache_version} = $1 + 0;
    }
    # The transcript set is the species directory's suffix (CacheDir.pm:312
    # appends `_refseq` or `_merged` to the species under --refseq or --merged);
    # Ensembl VEP's info.txt has no line for it, and a cache outside that
    # layout declares none.
    if ($cache_dir =~ m{(?:^|/)([^/]+)/\d+(?:_[^/]+)?/?$}) {
        my $species_dir = $1;
        $info{source_type} = $species_dir =~ /_refseq$/ ? 'refseq'
            : $species_dir =~ /_merged$/ ? 'merged'
            : 'ensembl';
    }
    return \%info;
}

# Autosomes by number, then X, Y and MT, then every other contig by name. The
# numeric comparison runs on the copies contig_sort_key returns: sort aliases
# its two operands to the contig list's own elements, and a `<=>` on those
# marks them as numbers, which JSON::PP then writes unquoted in contig_order.
sub contig_sort_key {
    my ($contig) = @_;
    return (0, $contig) if $contig =~ /^\d+$/;
    return (1, $contig) if $contig eq 'X';
    return (2, $contig) if $contig eq 'Y';
    return (3, $contig) if $contig eq 'MT';
    return (4, $contig);
}

sub by_contig_order {
    my ($tier_a, $key_a) = contig_sort_key($a);
    my ($tier_b, $key_b) = contig_sort_key($b);
    return $tier_a <=> $tier_b if $tier_a != $tier_b;
    return $tier_a == 0 ? $key_a <=> $key_b : $key_a cmp $key_b;
}

sub write_info_json {
    my ($info, $out_path) = @_;
    open(my $ofh, ">", $out_path) or fail("cannot write $out_path: $!");
    print $ofh JSON->new->canonical->pretty->encode($info);
    unless (close($ofh)) {
        my $error = $!;
        unlink $out_path;
        fail("write of $out_path failed: $error");
    }
    print "Wrote $out_path\n";
}
