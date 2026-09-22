#!/usr/bin/perl
# Storable-to-JSON converter that produces clean JSON matching vep-rs json_cache.rs.
# Whitelists ONLY the fields the Rust parser expects, stripping all Perl internals.
use strict;
use warnings;
use Storable qw(fd_retrieve);
use JSON;
use File::Path qw(make_path);
use Scalar::Util qw(blessed reftype refaddr);

my $cache_dir = $ARGV[0] or die "Usage: $0 cache_dir output_dir [chr1,chr2,...]\n";
my $output_dir = $ARGV[1] or die "Usage: $0 cache_dir output_dir [chr1,chr2,...]\n";
my @chrs = $ARGV[2] ? split(/,/, $ARGV[2]) : ();

# --- Field whitelists (matching json_cache.rs structs) ---

my %TRANSCRIPT_FIELDS = map { $_ => 1 } qw(
    stable_id version dbID start end strand biotype source description
    gene_symbol gene_symbol_source gene_hgnc_id gene_phenotype gene_stable_id
    is_canonical mane_select mane_plus_clinical tsl appris ccds
    protein refseq swissprot trembl uniparc chr
    cdna_coding_start cdna_coding_end coding_region_start coding_region_end
);
# These also map to transcript fields with underscored Perl names
my %TRANSCRIPT_RENAME = (
    '_gene_stable_id'    => 'gene_stable_id',
    '_gene_symbol'       => 'gene_symbol',
    '_gene_symbol_source'=> 'gene_symbol_source',
    '_gene_hgnc_id'      => 'gene_hgnc_id',
    '_gene_phenotype'    => 'gene_phenotype',
    '_ccds'              => 'ccds',
    '_refseq'            => 'refseq',
    '_gene_version'      => 'gene_version',
    # Protein-database cross-refs are stored on the Perl transcript with
    # underscore-prefixed keys. Map them to the unprefixed name vep-rs expects.
    '_protein'           => 'protein',
    '_swissprot'         => 'swissprot',
    '_trembl'            => 'trembl',
    '_uniparc'           => 'uniparc',
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
    # analysis: extract display_label if it's an object
    my $analysis = $pf->{analysis};
    if (ref $analysis eq 'HASH' || (ref $analysis && eval { $analysis->{display_label} })) {
        $out{analysis} = $analysis->{display_label};
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

sub extract_mapper {
    my ($mapper, $sorted_exons, $strand) = @_;
    return undef unless ref $mapper;
    my %out;
    for my $k (qw(start_phase cdna_coding_start cdna_coding_end)) {
        $out{$k} = $mapper->{$k} if defined $mapper->{$k};
    }

    # Try to get pairs from the mapper directly (GRCh37 format)
    if (ref $mapper->{pairs} eq 'ARRAY' && @{$mapper->{pairs}}) {
        $out{pairs} = [map { extract_mapper_pair($_) } grep { defined $_ && ref $_ } @{$mapper->{pairs}}];
        $out{pair_count} = scalar @{$out{pairs}};
    }
    # Otherwise, compute pairs from sorted exons (GRCh38 Storable format)
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

sub extract_prediction_matrix {
    my ($pm) = @_;
    return undef unless ref $pm;
    my %out;
    for my $k (qw(analysis sub_analysis peptide_length predictions_data)) {
        $out{$k} = $pm->{$k} if defined $pm->{$k};
    }
    return \%out;
}

sub extract_seq_string {
    # Extract a literal sequence string from a value that may be a plain string,
    # a blessed Bio::EnsEMBL::Slice (has ->seq method), a Bio::Seq, or a HASH
    # with a {seq} or {primary_seq}{seq} field. Returns undef when no string
    # can be recovered.
    #
    # Perl VEP's _variation_effect_feature_cache->{three_prime_utr} stores a
    # Bio::EnsEMBL::Slice ref blessed object whose ->seq() reads from the
    # cache's primary_seq region. Without this extraction, JSON encoding emits
    # null and vep-rs falls back to the clamped-peptide path, which emits
    # `stop_lost` where Perl emits a 3'UTR term. Extracting the seq string keeps
    # vep-rs offline+cache-only mode self-contained (no FASTA dependency).
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
    # five/three_prime_utr are blessed Bio::EnsEMBL::Slice objects in Storable;
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
                $pfp_out{$key} = extract_prediction_matrix($pfp->{$key});
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

    # Direct fields
    for my $k (keys %TRANSCRIPT_FIELDS) {
        $out{$k} = $t->{$k} if defined $t->{$k};
    }
    # Renamed fields (underscore-prefixed Perl internals)
    for my $pk (keys %TRANSCRIPT_RENAME) {
        my $rk = $TRANSCRIPT_RENAME{$pk};
        $out{$rk} = $t->{$pk} if defined $t->{$pk} && !defined $out{$rk};
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

opendir(my $dh, $cache_dir) or die "Cannot open $cache_dir: $!\n";
my @chr_dirs = sort grep { -d "$cache_dir/$_" && $_ !~ /^\./ } readdir($dh);
closedir($dh);

# Only standard chromosomes
@chr_dirs = grep { /^(\d+|X|Y|MT)$/ } @chr_dirs;
if (@chrs) {
    my %wanted = map { $_ => 1 } @chrs;
    @chr_dirs = grep { $wanted{$_} } @chr_dirs;
}

print "Processing " . scalar(@chr_dirs) . " chromosomes\n";

my $total = 0;
my $ok = 0;
my $transcripts_total = 0;

foreach my $chr (@chr_dirs) {
    my $chr_path = "$cache_dir/$chr";
    my $out_chr = "$output_dir/transcripts/$chr";
    make_path($out_chr);

    opendir(my $ch, $chr_path) or next;
    my @files = sort grep { /^\d+-\d+\.gz$/ } readdir($ch);
    closedir($ch);

    $total += scalar @files;
    my $chr_tx = 0;

    foreach my $f (@files) {
        my $in = "$chr_path/$f";
        (my $base = $f) =~ s/\.gz$/.json/;
        my $out = "$out_chr/$base";

        eval {
            open my $fh, "gzip -dc $in |" or die "gzip pipe failed: $!";
            my $data = fd_retrieve($fh);
            close($fh);

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
            close($ofh);
            $ok++;
        };
        if ($@) {
            warn "WARN: $in: $@\n";
        }
    }
    $transcripts_total += $chr_tx;
    print "  chr$chr: " . scalar(@files) . " files, $chr_tx transcripts\n";
}
# info.json: the cache metadata vep-rs prints in its output headers, translated
# from the VEP cache's info.txt the way VEP's own reader does (CacheDir.pm
# read_cache_info_file): every `source_<key>` line becomes a source version,
# `variation_cols` becomes a list, `-` means unset, `cell_types` is not carried.
write_info_json("$cache_dir/info.txt", "$output_dir/info.json");

print "\nDone. Converted: $ok / $total files, $transcripts_total total transcripts\n";

sub write_info_json {
    my ($info_txt, $out_path) = @_;
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
                $info{$key} = ($value ne '' && $value ne '0' && $value ne '-') ? JSON::true : JSON::false;
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
    open(my $ofh, ">", $out_path) or die "write failed: $!";
    print $ofh JSON->new->canonical->pretty->encode(\%info);
    close($ofh);
    print "Wrote $out_path\n";
}
