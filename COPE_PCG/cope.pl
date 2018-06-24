#!/usr/bin/perl
#
# cope.pl
#
# Copyright (C) Center for Bioinformatics, Peking University <cope@mail.cbi.pku.edu.cn>
#
# The software provided herein is free for ACADEMIC INSTRUCTION AND
# RESEARCH USE ONLY. You are free to download, copy, compile, study,
# and refer to the source code for any personal use of yours. Usage by
# you of any work covered by this license should not, directly or
# indirectly, enable its usage by any other individual or
# organization.
#
# You are free to make any modifications to the source covered by this
# license. You are also free to compile the source after modifying it
# and using the compiled product obtained thereafter in compliance
# with this License.  You may NOT under any circumstance copy,
# redistribute and/or republish the source or a work based on it
# (which includes binary or object code compiled from it) in part or
# whole without the permission of the authors.
#
# If you intend to incorporate the source code, in part or whole, into
# any free or proprietary program, you need to explicitly write to the
# original author(s) to ask for permission via e-mail at
# cope@mail.cbi.pku.edu.cn.
#
# Commercial licenses are available to legal entities, including
# companies and organizations (both for-profit and non-profit),
# requiring the software for general commercial use. To obtain a
# commercial license please, contact us via e-mail at
# cope@mail.cbi.pku.edu.cn.
#
# DISCLAIMER
#
# This software is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# THE ORIGINAL AUTHOR OF THE PROGRAM IS NOT LIABLE TO YOU FOR DAMAGES,
# INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL OR CONSEQUENTIAL DAMAGES
# ARISING OUT OF THE USE OR INABILITY TO USE THE PROGRAM (INCLUDING
# BUT NOT LIMITED TO LOSS OF DATA OR DATA BEING RENDERED INACCURATE OR
# LOSSES SUSTAINED BY YOU OR THIRD PARTIES OR A FAILURE OF THE PROGRAM
# TO OPERATE WITH ANY OTHER PROGRAMS)

use strict;
use Getopt::Long;
use Bio::DB::Fasta;

## config
my $bedtools_path="/lustre/tools_centos6.3/tools";
my $needle_path="/rd/build/EMBOSS/install/bin";
my $VERSION = '1.0';

## option
my $sample;
my $variant_file;
my $unphase; ## force for unphased input variants
my $process;
my $context; ## local region,default -/+100bp
my $splicegain; ## option default off
my $ensembl; ## ensembl gene model
my $hg38; ## version of reference genome
my $bed;
my $score;
my $confident;

sub execute
{
    my $cmd = shift;
    system($cmd) == 0 or die "[cope.pl] ERROR executing command: $cmd\n.";
}

GetOptions ("sample=s" => \$sample,
	    "variant=s" => \$variant_file,
	    "unphase!" => \$unphase,
	    "thread:i" => \$process,
	    "context:i" => \$context,
	    "splicegain!" => \$splicegain,
	    "ensembl!" => \$ensembl,
	    "hg38!" => \$hg38,
	    "bed!" => \$bed,
	    "score!" => \$score,
	    "splicing_confident_score!" => \$confident,
	   );
unless($process)
{
    $process=10; ##default using 10 cpus
}
unless($context)
{
    $context=100;
}
unless($splicegain)
{
    $splicegain=0; ## default off
}
unless($bed)
{
    $bed=0; ## default off
}
unless($score)
{
    $score=0;
}
unless($confident)
{
    $confident=0;
}

my $gene_resource="refseq";
my $gene_version="hg19";
my $reference="hg19.fa";
if($ensembl ==1) ## use ensembl gene model instead of RefSeq gene model
{
    $gene_resource="ensembl";
}
if($hg38 == 1)
{
    $gene_version="hg38";
    $reference="hg38.fa";
}
my $gene_model_bed=$gene_resource."_".$gene_version."_filtered.bed";
my $gene_model_fa=$gene_resource."_".$gene_version."_filtered.fa";

my $usage="perl cope.pl --sample sample --variant VCF.file [--unphase] [-thread int] [-c int] [--splicegain] [--ensembl] [--hg38] [--bed] [--score] [--splicing_confident_score]\n";
unless(-d "$sample" )
{
    mkdir $sample;
}

unless($variant_file && $sample)
{
    die $usage;
}

my $phased_snp_indel_file; ## final vcf format file
## haplotype information
if($unphase==0)
{
    print "Phased input variants.\n";
    #execute("cp input/$variant_file $sample/$variant_file");
    $phased_snp_indel_file=$variant_file;
}
else
{
    print "Unphase input variants. Variants are forced to be '0|1'.\n";
    execute("perl script/force_phase_information.pl input/$variant_file input/$variant_file.forced.vcf");
    $phased_snp_indel_file="$variant_file.forced.vcf";
}
## filter conflict SNPs and Indels
print "Filtering conflict SNPs and Indels\n";
execute("perl script/filter_conflict_variant_after_phasing.pl input/$phased_snp_indel_file");

## annotate protein coding variants
print "Running protein-coding genes annotation\n";
execute("$bedtools_path/intersectBed -a input/$phased_snp_indel_file.remove.Conflicts.vcf -b resource/$gene_model_bed -wa|uniq > $sample/final.phased.vcf.remove.Conflicts.WithinGene.vcf");

execute("perl script/extract_novel_transcript_and_sequence_alignment.pl -needle $needle_path -aminoref $gene_model_fa -r $reference -v $sample/final.phased.vcf.remove.Conflicts.WithinGene.vcf -g resource/$gene_model_bed -p $process -context $context -splicegain $splicegain -bed $bed -score $score -confident $confident > $sample/final.protein.annotate.txt");

