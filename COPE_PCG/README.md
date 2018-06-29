# Manual for COPE-PCG

## Attention
In the folder "resource", the file "hg38.fa" is required. The download link will be attached here later. 

## Requrements

COPE is mainly written by Perl. Several Perl modules are required as followed.

    Bio::DB::Fasta
    Set::IntervalTree
    Getopt::Long
    Parallel::ForkManager

Besides, a function in COPE is implemented by a Python module NetworkX. Needle and bedtools are also included in the pipeline of COPE.

For variants without phase information, we recommend our phasing script which bases on HapCUT and SHAPEIT2.

## Configure

This will configure serveral paths for different requrements. You should replace the following paths in COPE.pl with the real destination in your own system.

  
    my $bedtools_path="/lustre/tools_centos6.3/tools";
    my $needle_path="/rd/build/EMBOSS/install/bin";

## Command Line

#### Usage:

    perl COPE.pl --sample sampleName --variant example.vcf [--unphase] [--thread <int>] [--context <int>] [--splicegain] [--ensembl] [--hg38]

There should be three subdirectories in your working directory: input, resource and script. All input files should be stored in the subdirectory named input. The annotation resources used by COPE are all stored in resource subdirectory. Particularly, users should put a fasta file with reference genome sequence named hg19.fa or hg38.fa into resource subdirectory. The scripts are stored in script subdirectory.  

#### Main Options  

 Flag 	 Alternate 	 Description  
 --sample 	-s 	The name of the Project. COPE will make a directory based on the option and the output of annotation will be stored in this fold.  
--variant 	-v 	The VCF must be sorted by chromosome and position. Besides, the VCF file should be put in the input fold.  --unphase 	-u 	This option force COPE to handle variants without phase information. Warning: we don't recommend this option. Default off  
--thread 	-t 	The number of CPUs used for prediction novel protein sequences. Default 10  
--context 	-c 	The region used for searching cryptic splicing sites. Default 100  
--splicegain 		The option is used to decide whether to predict novel splice site gain in the whole gene region.  
Default off  
--ensembl 		The option is used to choose the Ensembl gene model. (default RefSeq)  
--hg38 		The option is used to choose the hg38 version of human genome assembly. Default hg19  
--bed 		Write output in BED format. Default off  
--score 	Report the splice score change for splice-disrupting variant. Default off  
--splicing_confident_score 	Report the absolute splice motif score for identified novel potential splice site. Default off


## Result
The annotation result is stored in file final.protein.annotate.txt. The default output format is a tab-delimited text file. This is an example:

  
 NM_001005484	OR4F5	0|1	S	305	113-113:F->C;	1:69428:69428:T:G:1|1  
 NM_001005484	OR4F5	1|0	S	305	113-113:F->C;	1:69428:69428:T:G:1|1  
 NM_152228	TAS1R3	0|1	O	852	757-757:C->R;	1:1269554:1269554:T:C:1|1,1:1268847:1268847:T:G:1|1  
 NM_152228	TAS1R3	1|0	O	852	757-757:C->R;	1:1269554:1269554:T:C:1|1,1:1268847:1268847:T:G:1|1  
 
 
Each line contains the following data:

   1. Transcript ID.
   2. Gene name.
   3. Haplotype information. "0|1" means hap1 and "1|0" means hap2.
   4. State for splicing pattern.“O” stands for original splicing pattern; “R” stands for rescued splicing pattern; “G” means the transcript with novel splice site gain; “N” means the transcript without novel splice site. "S" indicates single exon transcprit. “G” and “N” are used with the option --splicegain.
   5. The length of reference protein sequence.
   6. Amino acid changes.
   7. Variants. Each variant is delimited by comma. Each variant six parts (chromosome,start,end,reference allele,alternative allele and haplotype) delimited by colon.

Besides, COPE also supports reporting in BED format by using --bed option. This is an example:

  
1	69090	70008 Transcript=NM_001005484;Gene=OR4F5;Haplotype=0|1;Splicing_Code=S;Protein_Length=305;Amino_Acid=113-113:F->C;Variant=1:69428:69428:T:G:1|1;  
1	69090	70008	Transcript=NM_001005484;Gene=OR4F5;Haplotype=1|0;Splicing_Code=S;Protein_Length=305;Amino_Acid=113-113:F->C;Variant=1:69428:69428:T:G:1|1;  
1	1266725	1269844	Transcript=NM_152228;Gene=TAS1R3;Haplotype=0|1;Splicing_Code=O;Protein_Length=852;Amino_Acid=757-757:C->R;Variant=1:1269554:1269554:T:C:1|1,1:1268847:1268847:T:G:1|1;  
1	1266725	1269844	Transcript=NM_152228;Gene=TAS1R3;Haplotype=1|0;Splicing_Code=O;Protein_Length=852;Amino_Acid=757-757:C->R;Variant=1:1269554:1269554:T:C:1|1,1:1268847:1268847:T:G:1|1;  

 
Each line contains the following data:

   1. Chromosome.
   2. The starting position of the transcript.
   3. The ending position of the transcript.
   4. The annotation result as shown before.

