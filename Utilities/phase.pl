#!/usr/bin/perl
use strict;
use Getopt::Long;
## usage:

## config for shapeit
my $PATH_for_extractPIRs="/lustre/user/chengsj/15_Haplotype_Assembly/07_shapeit_static/extractPIRs/extractPIRs/bin";
my $PATH_for_shapeit2="/lustre/user/chengsj/15_Haplotype_Assembly/07_shapeit_static";
my $REF_panel="/lustre/user/chengsj/15_Haplotype_Assembly/07_shapeit_static/1000GP_Phase3/1000GP_Phase3";

## config for hapcut
my $path_for_hapcut="/lustre/user/chengsj/15_Haplotype_Assembly/03_HapCut/HAPCUT-v0.7";

## option
my $sample;
my $variant_file;
my $bam_file;
my $gender;
my $reference;
sub execute
{
    my $cmd = shift;
    system($cmd) == 0 or die "[corect.pl] ERROR executing command: $cmd\n.";
}

GetOptions ("sample|s=s" => \$sample,
	    "variant|v=s" => \$variant_file,
	    "bam|b:s"  => \$bam_file,
	    "gender|g:s" => \$gender,
	    "reference|r:s" => \$reference,
	   );

my $usage="perl phase.pl -sample sample -variant VCF.file -bam bamfile -gender [male|female]] -reference hg19.fa\n";
unless(-d "$sample" )
{
    mkdir $sample;
}

unless($variant_file && $bam_file && $gender && $reference && $sample)
{
    die $usage;
}

my $phased_snp_indel_file; ## final vcf format file
## haplotype reconstruction
print "Run haplotype inference\n";
execute("./script/hapcut.sh $path_for_hapcut $reference input/$variant_file input/$bam_file");
##run shapeit2
execute("perl script/change_vcf_format.pl input/$variant_file $gender > input/snp_indel_merged_for_shapeit.vcf");
if($gender eq 'male')
{
    execute("./script/run_shapeit2.sh $PATH_for_extractPIRs $PATH_for_shapeit2 $REF_panel input/$bam_file input/snp_indel_merged_for_shapeit.vcf $gender");
    `rm -f *mm`;
    `rm -f *log`;
}
elsif($gender eq 'female')
{
    execute("./script/run_shapeit2_for_female.sh $PATH_for_extractPIRs $PATH_for_shapeit2 $REF_panel input/$bam_file input/snp_indel_merged_for_shapeit.vcf $gender");
    `rm -f *mm`;
    `rm -f *log`;
}
else
{
    die "gender is wrong!";
}
## use hapcut to fix the result
execute("cat ./tmp_PIRs/*haplotype.txt.haps > ./tmp_PIRs/all.shapeit.result");
execute("perl ./script/use_hapcut_to_add_phase.pl input/$variant_file input/$variant_file.output_haplotype_file $gender > ./tmp_PIRs/add.haps.result");
execute("cat ./tmp_PIRs/all.shapeit.result ./tmp_PIRs/add.haps.result > $sample/final.phased.result");
## change final phased result into VCF format
print "change final phased result into VCF format\n";
open(VCF,"$sample/final.phased.result") or die;
open(FORMAT,">$sample/final.phased.vcf") or die; ##vcf format #chr pos id ref alt qual filter info format gt
print FORMAT "##fileformat=VCFv4.1\n";
while(my $line=<VCF>)
{
    chomp $line;
    my @info=split /\s+/,$line;
    print FORMAT "$info[0]\t$info[2]\t$info[1]\t$info[3]\t$info[4]\t.\t.\t.\tGT\t$info[5]|$info[6]\n";
}
close VCF;
close FORMAT;
$phased_snp_indel_file="final.phased.vcf";
## add X and Y
if($gender eq 'male')
{
    open(X,"./tmp_VCF/X.snp.vcf") or die;
    open(Y,"./tmp_VCF/Y.snp.vcf") or die;
    open(TO,">>$sample/final.phased.vcf");
    while(my $line=<X>)
    {
	chomp $line;
	if($line=~/^#/)
	{
	    next;
	}
	else
	{
	    my @info=split /\t/,$line;
	    print TO "$info[0]\t$info[1]\t$info[2]\t$info[3]\t$info[4]\t$info[5]\t$info[6]\t$info[7]\t$info[8]\t0|1\n";
	}
    }
    close X;
    while(my $line=<Y>)
    {
	chomp $line;
	if($line=~/^#/)
	{
	    next;
	}
	else
	{
	    my @info=split /\t/,$line;
	    print TO "$info[0]\t$info[1]\t$info[2]\t$info[3]\t$info[4]\t$info[5]\t$info[6]\t$info[7]\t$info[8]\t0|1\n";
	}
    }
    close Y;
	close TO;
}

print "final.phased.vcf is your final phased result.";
