#!/usr/bin/perl
use strict;
use Getopt::Long;

## config
my $bedtools_path="/lustre/tools_centos6.3/tools";
my $VERSION = '1.0';

## option
my $variant_file;
my $unphase; ## force for unphased input variants
my $process;
my $gain; ## option default off
my $hg38; ## version of reference genome

sub execute
{
    my $cmd = shift;
    system($cmd) == 0 or die "[cope_tfbs.pl] ERROR executing command: $cmd\n.";
}

GetOptions ("variant=s" => \$variant_file,
	    "unphase!" => \$unphase,
	    "gain!" => \$gain,
	   );

unless($gain)
{
    $gain=0; ## default off
}

my $usage="perl cope_tfbs.pl --variant variant.vcf [--unphased] [--gain]";
unless($variant_file)
{
    die $usage;
}

my $phased_snp_indel_file;
if($unphase==0)
{
    print "Phased input variants.\n";
    $phased_snp_indel_file=$variant_file;
}
else
{
    print "Unphase input variants. Variants are forced to be '0|1'.\n";
    execute("perl script/force_phase_information.pl $variant_file $variant_file.forced.vcf");
    $phased_snp_indel_file="$variant_file.forced.vcf";
}

## filter conflict SNPs and Indels
print "Filtering conflict SNPs and Indels\n";
execute("perl script/filter_conflict_variant_after_phasing.pl $phased_snp_indel_file");

print "Running TFBS breaking prediction\n";
execute("perl script/TFBS_Breaking_annotation.pl $phased_snp_indel_file.remove.Conflicts.vcf data/ENCODE.tf.bound.union.sorted.bed.gz > TF_Breaking.txt");

if($gain==1)
{
    print "Running Gain of TFBS prediction\n";
    execute("perl script/promoter_annotation.pl $phased_snp_indel_file.remove.Conflicts.vcf data/hg19_refseq_NMgene_promoter.sorted.bed.gz > TF_Gain.txt");
}
unlink "$phased_snp_indel_file.remove.Conflicts.vcf";
unlink "$variant_file.forced.vcf";
