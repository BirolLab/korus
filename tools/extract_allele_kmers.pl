#!/usr/bin/env perl
#RLW2026
use strict;
use warnings;

if ($#ARGV < 1) {
    die "Usage: $0 <GRCh38.fa(.gz)> <k>\n";
}

my $refg = $ARGV[0];
my $k    = $ARGV[1];
my $half = int($k / 2);

# =========================================================
# Layer 3: Canonical k-mer (INDEXING ONLY)
# =========================================================
sub canonical {
    my ($kmer) = @_;

    my $rc = reverse($kmer);
    $rc =~ tr/ACGT/TGCA/;

    return ($kmer lt $rc) ? $kmer : $rc;
}

# =========================================================
# Load genome (unchanged, forward strand reference)
# =========================================================
my %genome;
my $chr = "";

#open(my $FA, "unpigz -c $refg |") or die $!;
my $cmd;

if ($refg =~ /\.gz$/) {
    $cmd = "unpigz -c $refg";
}
elsif ($refg =~ /\.bz2$/) {
    $cmd = "bunzip2 -c $refg";
}
elsif ($refg =~ /\.(fa|fasta)(\.gz|\.bz2)?$/ || $refg !~ /\.gz$|\.bz2$/) {
    # try plain file first (no decompression)
    $cmd = "cat $refg";
}
else {
    die "Unsupported reference format: $refg";
}

open(my $FA, "$cmd |") or die "Failed to open reference with command: $cmd";


while (<$FA>) {
    chomp;

    if (/^>(\S+)/) {
        $chr = $1;
        $genome{$chr} = "";
    } else {
        $genome{$chr} .= uc($_);
    }
}
close $FA;

# =========================================================
# State
# =========================================================
my $cur_chr = "";
my $cur_pos = 0;

my $ref_kmer;
my %alt_seen;

# =========================================================
# Layer 1: Extract fixed window (TRUTH)
# =========================================================
sub extract_window {
    my ($chr, $pos) = @_;

    return unless defined $genome{$chr};

    my $start = $pos - $half - 1;

    return if $start < 0;
    return if ($start + $k) > length($genome{$chr});

    my $seq = substr($genome{$chr}, $start, $k);

    return ($seq, $start);
}

# =========================================================
# Layer 2: Apply variant (BIOLOGICAL embedding)
# =========================================================
sub apply_variant {
    my ($seq, $ref, $alt, $is_snv) = @_;

    die "k-mer length mismatch"
        unless length($seq) == $k;

    # strict REF validation (only once per site)
    if (!defined $ref_kmer) {
        die "REF mismatch at $cur_chr:$cur_pos (expected $ref, got "
            . substr($seq, $half, length($ref)) . ")"
            unless substr($seq, $half, length($ref)) eq $ref;
    }

    # -----------------------------------------------------
    # SNV
    # -----------------------------------------------------
    if ($is_snv) {
        substr($seq, $half, 1) = $alt;
        return $seq;
    }

    my $ref_len = length($ref);
    my $alt_len = length($alt);

    # split reference window around anchor
    my $left  = substr($seq, 0, $half);
    #my $right = substr($seq, $half + $ref_len);
    my $right = "";

    if (($half + $ref_len) < length($seq)) {
        $right = substr($seq, $half + $ref_len);
    }

    # -----------------------------------------------------
    # INSERTION (ALT longer than REF)
    # -----------------------------------------------------
    if ($alt_len > $ref_len) {

        my $ins = substr($alt, 1);   # inserted sequence (after anchor base)

        my $new = $left . substr($alt, 0, 1) . $ins . $right;

        # re-encode to fixed k by truncating right side
        return substr($new, 0, $k);
    }

    # -----------------------------------------------------
    # DELETION (ALT shorter than REF)
    # -----------------------------------------------------
    if ($alt_len < $ref_len) {


        # deletion removes REF bases at center,
        # so collapse REF region and refill from genomic downstream sequence

        my $new = $left . substr($alt, 0, 1) . $right;

        # genomic coordinate immediately AFTER deleted REF
        my $fill_start = ($cur_pos - 1) + $ref_len + length($right);

        my $fill_needed = $k - length($new);

        if ($fill_needed > 0) {
            $new .= substr($genome{$cur_chr}, $fill_start, $fill_needed);
        }

        return substr($new, 0, $k);

    }

    # -----------------------------------------------------
    # complex substitution-like INDEL
    # -----------------------------------------------------
    my $new = $left . substr($alt, 0, 1) . $right;

    return substr($new, 0, $k);
}


# =========================================================
# Flush site output
# =========================================================
sub flush_site {
    return unless defined $ref_kmer;

    #my $ref_can = canonical($ref_kmer);
    my $ref_can = $ref_kmer;

    print join("\t", $cur_chr, $cur_pos, "REF", $ref_can), "\n";

    for my $alt_can (values %alt_seen) {
        print join("\t", $cur_chr, $cur_pos, "ALT", $alt_can), "\n";
    }
}

# =========================================================
# Main VCF stream
# =========================================================
while (<STDIN>) {

    next if /^#/;
    chomp;

    my @f = split /\t/;
    my ($chr, $pos, undef, $ref, $alt) = @f[0..4];

    # site boundary
    if ($chr ne $cur_chr || $pos != $cur_pos) {
        flush_site();

        $cur_chr = $chr;
        $cur_pos = $pos;

        %alt_seen = ();
        $ref_kmer  = undef;
    }

    my ($seq, $start) = extract_window($chr, $pos);
    next unless defined $seq;

    my @alts = split /,/, $alt;

    foreach my $a (@alts) {

        my $is_snv = (length($ref) == 1 && length($a) == 1);

        # Build REF once (only SNVs used as anchor)
	#if (!defined $ref_kmer && $is_snv) {
	#    next unless substr($seq, $half, 1) eq $ref;
	#    $ref_kmer = $seq;
	#}


	 # Build REF k-mer ONCE (but now safe for ALL variant types)
        if (!defined $ref_kmer) {
            next unless substr($seq, $half, length($ref)) eq $ref;
            $ref_kmer = $seq;
        }

        next unless defined $ref_kmer;

        # Apply ALT in forward space
        my $alt_seq = apply_variant($seq, $ref, $a, $is_snv);
	#$alt_seq = substr($alt_seq,0,$k);### deal with pesky problem
        die "ALT k-mer wrong length"
            unless length($alt_seq) == $k;


        # Canonicalize ONLY for indexing/output
	#my $alt_can = canonical($alt_seq);
        my $alt_can = $alt_seq;

	#$alt_seen{$a . "|" . $alt_can} = $alt_can;
	
	$alt_seen{$ref . ">" . $a . "|" . $alt_can} = $alt_can;
    }
}

flush_site();
