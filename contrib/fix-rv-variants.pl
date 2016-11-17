#!/usr/bin/perl -w

# perl packages needed
use strict;
use Compress::Zlib;

# usage
sub usage()
{
print STDERR << "EOF";

usage: perl $0 -r -i infile > outfile

processes a gzipped dbSNP file looking for "RV" code lines and reversing their bases
download it from ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/00-All.vcf.gz

by Jorge Amigo Lechuga

OPTIONS
    -i    gzipped input file
    -r    remove "RV" codes
    -h    this (help) message
EOF
exit;
}

# get options
use Getopt::Std;
my %opt;
my $opt_string = 'i:rh';
getopts( "$opt_string", \%opt ) or usage();
usage() if $opt{h};
my $inputfile = $opt{i} ? $opt{i}: usage();
my $rvremoval = $opt{r} ? 1 : 0;

# define the base translations
my %trans = qw(A T T A C G G C);
my $check = join '|', keys %trans;

# process the gzipped dbSNP file
my $gz = gzopen($inputfile, "rb");
while ( $gz->gzreadline(my $line) > 0 ) {
    if ($line =~ /RV;/) {
        my @cols = split /\t/, $line;
        # reverse allele columns
        $cols[3] =~ s/($check)/$trans{$1}/g;
        $cols[4] =~ s/($check)/$trans{$1}/g;
        # remove "RV" code
        $cols[7] =~ s/RV;//g if $rvremoval;
        print join "\t", @cols;
    } else {
        print $line;
    }
}
$gz->gzclose();

