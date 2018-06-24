# Haplotype Inference
## Configure

This will configure serveral paths for different requrements. You should replace the following paths in phase.pl with the real destination in your own system.

    my $PATH_for_extractPIRs="PATH_for_extractPIRs/extractPIRs/extractPIRs/bin";
    my $PATH_for_shapeit2="PATH_for_shapeit2";
    my $REF_panel="PATH_for_shapeit2/1000GP_Phase3";
    my $path_for_hapcut="path_for_hapcut";


## Command Line

#### Usage:

    perl phase.pl --sample sampleName --variant example.vcf --bam bamFile --gender [male|female] --reference human_genome.fa

There should be two subdirectories in your working directory: input and script. All input files should be stored in the subdirectory named input. The scripts are stored in script subdirectory.

#### Main Options
Flag	 Alternate	 Description  
--sample	-s	The name of the Project. phase.pl will make a directory based on the option and the output will be stored in this fold.  
--variant	-v	The VCF must be sorted by chromosome and position. Besides, the VCF file should be put in the input fold.  
--bam		-b	The BAM file used for HapCUT. The file should be put in the input fold.  
--gender	-g	The gender of your sample [male|female].  
--reference	-r	The reference genome sequence in FASTA format.  

