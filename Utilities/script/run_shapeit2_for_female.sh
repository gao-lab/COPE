#!/bin/bash

#useage: ./run_shapeit2.sh bam_file vcf_file
##this bash script is used to run shapeit for all chromosome
PATH_for_extractPIRs=$1 #"/lustre/user/chengsj/15_Haplotype_Assembly/07_shapeit_static/extractPIRs/extractPIRs/bin"
PATH_for_shapeit2=$2 #"/lustre/user/chengsj/15_Haplotype_Assembly/07_shapeit_static"
REF_panel=$3 #"/lustre/user/chengsj/15_Haplotype_Assembly/07_shapeit_static/1000GP_Phase3"

bam=$4
vcf=$5
sampleid=$6 ## same as sampleID in VCF file
##vcf format: CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsampleID
##22      16053444        .       A       T       451.01  PASS    .       GT      0/1
if [ ! -d "./tmp_VCF" ]; then
  mkdir ./tmp_VCF
fi

if [ ! -d "./tmp_PIRs" ]; then
  mkdir ./tmp_PIRs
fi

if [ ! -d "./tmp_BAM" ]; then
  mkdir ./tmp_BAM
fi

for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
do
{
    echo "$sampleid	$bam	$chr" > ./tmp_BAM/$chr.bamlist
    echo "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	$sampleid" >./tmp_VCF/$chr.snp.vcf
    grep "^$chr\b" $vcf|grep -v ",">> ./tmp_VCF/$chr.snp.vcf
     vcf-sort -c ./tmp_VCF/$chr.snp.vcf > ./tmp_VCF/$chr.snp.sorted.vcf
    gzip ./tmp_VCF/$chr.snp.sorted.vcf
    $PATH_for_extractPIRs/extractPIRs --bam ./tmp_BAM/$chr.bamlist --vcf ./tmp_VCF/$chr.snp.sorted.vcf.gz --out ./tmp_PIRs/$chr.PIRs 
    $PATH_for_shapeit2/shapeit -check --input-vcf ./tmp_VCF/$chr.snp.sorted.vcf.gz -R $REF_panel/1000GP_Phase3_chr$chr.hap.gz $REF_panel/1000GP_Phase3_chr$chr.legend.gz $REF_panel/1000GP_Phase3.sample --output-log ./tmp_PIRs/$chr.alignmentChecks.log 
    $PATH_for_shapeit2/shapeit -assemble --input-vcf ./tmp_VCF/$chr.snp.sorted.vcf.gz --input-pir ./tmp_PIRs/$chr.PIRs -R $REF_panel/1000GP_Phase3_chr$chr.hap.gz $REF_panel/1000GP_Phase3_chr$chr.legend.gz $REF_panel/1000GP_Phase3.sample -O ./tmp_PIRs/$chr.haplotype.txt --exclude-snp ./tmp_PIRs/$chr.alignmentChecks.snp.strand.exclude --no-mcmc
}
done

#REF_panel_X="/lustre/user/chengsj/15_Haplotype_Assembly/07_shapeit_static/ALL_1000G_phase1integrated_v3_impute"
echo "$sampleid\t$bam\tX" > ./tmp_BAM/X.bamlist
echo "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$sampleid" >./tmp_VCF/X.snp.vcf
grep "^X\b" $vcf|grep -v ",">> ./tmp_VCF/X.snp.vcf
vcf-sort -c ./tmp_VCF/X.snp.vcf > ./tmp_VCF/X.snp.sorted.vcf
gzip ./tmp_VCF/X.snp.sorted.vcf
$PATH_for_extractPIRs/extractPIRs --bam ./tmp_BAM/X.bamlist --vcf ./tmp_VCF/X.snp.sorted.vcf.gz --out ./tmp_PIRs/X.PIRs
$PATH_for_shapeit2/shapeit -check --input-vcf ./tmp_VCF/X.snp.sorted.vcf.gz -R $REF_panel/1000GP_Phase3_chrX_NONPAR.hap.gz $REF_panel/1000GP_Phase3_chrX_NONPAR.legend.gz $REF_panel/1000GP_Phase3.sample --output-log ./tmp_PIRs/chrX_nonPAR.log 
$PATH_for_shapeit2/shapeit -assemble --input-vcf ./tmp_VCF/X.snp.sorted.vcf.gz --input-pir ./tmp_PIRs/X.PIRs -R $REF_panel/1000GP_Phase3_chrX_NONPAR.hap.gz $REF_panel/1000GP_Phase3_chrX_NONPAR.legend.gz $REF_panel/1000GP_Phase3.sample -O ./tmp_PIRs/X_nonPAR.haplotype.txt --exclude-snp ./tmp_PIRs/chrX_nonPAR.snp.strand.exclude --no-mcmc

$PATH_for_shapeit2/shapeit -check --input-vcf ./tmp_VCF/X.snp.sorted.vcf.gz -R $REF_panel/1000GP_Phase3_chrX_PAR1.hap.gz $REF_panel/1000GP_Phase3_chrX_PAR1.legend.gz $REF_panel/1000GP_Phase3.sample --output-log ./tmp_PIRs/chrX_PAR1.log
$PATH_for_shapeit2/shapeit -assemble --input-vcf ./tmp_VCF/X.snp.sorted.vcf.gz --input-pir ./tmp_PIRs/X.PIRs -R $REF_panel/1000GP_Phase3_chrX_PAR1.hap.gz $REF_panel/1000GP_Phase3_chrX_PAR1.legend.gz $REF_panel/1000GP_Phase3.sample  -O ./tmp_PIRs/X_PAR1.haplotype.txt --exclude-snp ./tmp_PIRs/chrX_PAR1.snp.strand.exclude --no-mcmc

$PATH_for_shapeit2/shapeit -check --input-vcf ./tmp_VCF/X.snp.sorted.vcf.gz -R $REF_panel/1000GP_Phase3_chrX_PAR2.hap.gz $REF_panel/1000GP_Phase3_chrX_PAR2.legend.gz $REF_panel/1000GP_Phase3.sample --output-log ./tmp_PIRs/chrX_PAR2.log
$PATH_for_shapeit2/shapeit -assemble --input-vcf ./tmp_VCF/X.snp.sorted.vcf.gz --input-pir ./tmp_PIRs/X.PIRs -R $REF_panel/1000GP_Phase3_chrX_PAR2.hap.gz $REF_panel/1000GP_Phase3_chrX_PAR2.legend.gz $REF_panel/1000GP_Phase3.sample -O ./tmp_PIRs/X_PAR2.haplotype.txt --exclude-snp ./tmp_PIRs/chrX_PAR2.snp.strand.exclude --no-mcmc

cat ./tmp_PIRs/X_nonPAR.haplotype.txt.haps ./tmp_PIRs/X_PAR1.haplotype.txt.haps ./tmp_PIRs/X_PAR2.haplotype.txt.haps > ./tmp_PIRs/X.haplotype.txt.haps
cat ./tmp_PIRs/chrX_nonPAR.snp.strand.exclude ./tmp_PIRs/chrX_PAR1.snp.strand.exclude ./tmp_PIRs/chrX_PAR2.snp.strand.exclude > ./tmp_PIRs/merged.X.snp.strand.exclude
sort ./tmp_PIRs/merged.X.snp.strand.exclude |uniq -c |cut -d ' ' -f 7-8|grep ^3 |cut -d ' ' -f 2 > ./tmp_PIRs/X.alignmentChecks.snp.strand.exclude
##perl /lustre/user/chengsj/15_Haplotype_Assembly/07_shapeit_static/script/extract_miss_variant_in_X.pl ./tmp_PIRs/merged.X.snp.strand.exclude> ./tmp_PIRs/X.alignmentChecks.snp.strand.exclude
