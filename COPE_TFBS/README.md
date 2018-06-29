# Manual for COPE-TFBS

## Attention
In the folder "data", the file "ENCODE.tf.bound.union.sorted.bed.gz" is required. The download link will be attached here later.

## Requrements

COPE-TFBS is written by Perl. Several Perl modules and other tools are required as followed.

Perl modules:  

    Bio::DB::Fasta
    Getopt::Long
    
Tools:  

    tabix
    bedtools
    
For variants without phase information, we recommend our phasing script which bases on HapCUT and SHAPEIT2.

## Configure

This will configure serveral paths for different requrements. You should replace the following paths in cope_tfbs.pl with the real destination in your own system.

     my $bedtools_path="/lustre/tools_centos6.3/tools";


## Command Line

#### Usage:

	 perl cope_tfbs.pl --variant variant.vcf [--unphased] [--gain]

There should be three subdirectories in your working directory: pwm, resource and script. All PWM files for TFs are stored in the subdirectory named pwm. The annotation resources used by COPE-tfbs are all stored in data subdirectory. Particularly, users should put a fasta file with reference genome sequence named hg19.fa into data subdirectory. The scripts are stored in script subdirectory.

#### Main Options
Flag	 Alternate	 Description  
--variant 		 -v	The VCF must be sorted by chromosome and position. Besides, the VCF file should be put in the input fold.  
--unphase    -u		 This option force COPE to handle variants without phase information. Warning: we don't recommend this option. Default off.  
--gain	     		 The option is used to decide whether to predict novel TFBS in the whole promoter region. Default off.


# Result
The annotation result is stored in file TF_annotation.txt. As our method is designed to annotate multiple variant effects on TFBS accurately, variant-centric format is no longer suitable for our output, such as VCF format. The default output format is a tab-delimited text file.


***
* Summary Statistics for TFBS Gain
* Number of variants uploaded: 4
* Number of variants creating novel TFBSs: 2
* Number of Novel TFBSs: 1
* Number of multi-variant novel TFBSs: 1
* Number of transcripts affected: 2

Chr	 Start		End	  TF	Sequence	PWM_score	Transcript	Strand	Variant  
5	 1295224	1295229	  ETS_known8		TTTCCG		5.74321215572648	NM_198253(Name=TERT)	-	5_1295228_G/A_0|1;5_1295229_G/A_0|1  
5	 1295224	1295229	  ETS_known8		TTTCCG		5.74321215572648	NM_001193376(Name=TERT)	-	5_1295228_G/A_0|1;5_1295229_G/A_0|1  

***
* Summary Statistics for TFBS Breaking
* Number of variants uploaded: 4
* Number of variants within TFBSs: 2
* Number of TFBSs affected: 3
* Number of TFBSs with multi-variants: 2

Chr	 Start	  End  TF	       Ref_PWM_score	Alt_PWM_score	Strand	Nearest_gene	Variant  
10	 100091757     100091764       HIF1A::ARNT_1	6.028847	5.048077		-	NM_138469(Name=R3HCC1L;Distance=142243bp)	10_100091761_G/A_0|1  
10	 100091753     100091766       HIF1A_2		9.043477	7.913043		-	NM_138469(Name=R3HCC1L;Distance=142242bp)	10_100091754_G/A_0|1;10_100091761_G/A_0|1  
10	 100091754     100091765       HIF1A_1		7.916666	6.916666		-	NM_138469(Name=R3HCC1L;Distance=142242bp)	10_100091754_G/A_0|1;10_100091761_G/A_0|1  

