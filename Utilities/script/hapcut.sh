#!/bin/bash
path_for_hapcut=$1
ref=$2
vcf=$3
bam=$4
vcf-sort -c $vcf > $vcf.sorted.vcf

$path_for_hapcut/extractHAIRS --VCF $vcf.sorted.vcf --bam $bam  --ref $ref --indels 1 > $vcf.fragment_matrix_file
$path_for_hapcut/HAPCUT --fragments $vcf.fragment_matrix_file --VCF $vcf.sorted.vcf --output $vcf.output_haplotype_file --maxmem 12000 > hapcut.log;printf "HapCut is finished!\n"

printf "Running SHAPEIT2\n"

