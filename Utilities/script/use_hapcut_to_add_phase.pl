#!/usr/bin/perl
use warnings;

## usage:
##this script is used to add SNV/Indel hapcut handled into shapeit result
##load all hom SNV into memory(some hom SNVs are took into excluded snvs)
my $vcf_file=$ARGV[0]; ##all variants in VCF format
my $gender=$ARGV[2];
my %hom_variant;
open(VCF,$vcf_file) or die "Cannot open it";
while(my $vcf=<VCF>)
{
    if($vcf=~/1\/1/)
    {
	my($chr,$posi,$ref,$alt)=(split /\t/,$vcf)[0,1,3,4];
	if($gender eq 'male' && $chr eq 'X')
	{
	    next;
	}
	elsif($gender eq 'male' && $chr eq 'Y')
	{
	    next;
	}
	$hom_variant{$chr ."_" .$posi}=$ref ." ".$alt;
    }
}
##loading all excluded snv into memory
my $i=1;
my %miss_snps;
unless($gender eq 'male')
{
    if(-e "./tmp_PIRs/X.alignmentChecks.snp.strand.exclude")
    {
	open(CHRX,"X.alignmentChecks.snp.strand.exclude") or die "cannot open it";
	while(<CHRX>)
	{
	    chomp;
	    my $variant="X_". $_;
	    if(defined($hom_variant{$variant}))
	    {
		print "X . $_ $hom_variant{$variant} 1 1\n";
	    }
	    else
	    {
		$miss_snps{$variant}=$_;
	    }
	}
    }
}
while($i<23)
{
    open(FILE,"./tmp_PIRs/$i.alignmentChecks.snp.strand.exclude") or die "cannot open it";
    my @exclude; #remove duplicate
    while(<FILE>)
    {
	chomp;
	unless($_ ~~ @exclude)
	{
	    my $variant=$i ."_". $_;
	    if(defined($hom_variant{$variant}))
	    {
		print "$i . $_ $hom_variant{$variant} 1 1\n";
	    }
	    else
	    {
		$miss_snps{$variant}=$_;
	    }
	    push @exclude,$_;
	}
    }
    close FILE;
    $i++;
}

##loading shapeit result into memory(all chromosome)
my %shapeit;
open(SHAPEIT,"./tmp_PIRs/all.shapeit.result") or die "cannot open it";
while(my $line=<SHAPEIT>)
{
    chomp $line;
    my($chr,$pos,$ref_nt,$alt_nt,$ref_type,$alt_type)=(split / /,$line)[0,2,3,4,5,6];
    $shapeit{$chr ."_".$pos."_".$ref_nt ."_". $alt_nt}=$ref_type ."|" .$alt_type;
}

## loading hapcut result (base on each block)
my $hapcut_result=$ARGV[1];
open(HAPCUT,$hapcut_result) or die "cannot open the hapcut file";
$/="******** \n";
my @add_phase;
while(my $line=<HAPCUT>)
{
    chomp $line;
    my @info=(split /\n/,$line);
    shift @info;
    my @var_hap_num;
    my $state_for_miss_snp=0;
    #my $state_for_indel=0;
    my @useful_variant;
    my @nearby_variant;
    foreach my $var ( @info )
    {
	my ($ref,$alt,$chr,$pos,$ref_alle,$alt_alle)=(split /\t/,$var)[1,2,3,4,5,6];
	my $key=$chr ."_".$pos;
	my $tmp=$chr ."_". $pos ."_".$ref_alle ."_". $alt_alle . "_". $ref ."|".$alt;
	if(defined($miss_snps{$key}) || $tmp=~/1\|2/ ||$tmp=~/2\|1/) ## consider hapcut can handle 1/2; while shapeit cannot.
	{
	    $state_for_miss_snp=1;
	    push @useful_variant,$tmp;
	}
	#elsif( length($ref_alle)>1 || length($alt_alle)>1 )
	#{
	 #   $state_for_indel=1;
	  #  push @useful_variant,$tmp;
	#}
	else
	{
	    push @nearby_variant,$tmp;
	}
    }
    if($state_for_miss_snp  >0 && @nearby_variant)
    {
	##find the nearby variant
	foreach my $variant ( @useful_variant )
	{
	    my $near_by_variant=nearest_variant($variant,@nearby_variant);
	    my $state_of_haplotype=0;# 0 stands for not same
	    if($variant=~/1\|2/ ||$variant=~/2\|1/)
	    {
		my $haplotype=(split /_/,$variant)[4];
		my $haplotype_for_variant;
		my($chromsome,$position,$reference,$alter)=(split /_/,$near_by_variant)[0,1,2,3];
		my $haplotype_for_nearby=$shapeit{$chromsome ."_".$position."_".$reference ."_". $alter};
		if($near_by_variant=~/$haplotype_for_nearby/)
		{
		    $haplotype_for_variant=$haplotype;
		}
		else
		{
		    $haplotype_for_variant= reverse $haplotype;
		}
		my($chrom,$pos,$ref_nt,$alter_nt)=(split /_/,$variant)[0,1,2,3];
		my($first,$second)=(split /\|/,$haplotype_for_variant)[0,1];
		print "$chrom . $pos $ref_nt $alter_nt $first $second\n";
	    }
	    else
	    {
		$state_of_haplotype =1 if $near_by_variant =~/1\|0/ && $variant =~/1\|0/;
		$state_of_haplotype =1 if $near_by_variant =~/0\|1/ && $variant =~/0\|1/;
		my($chromsome,$position,$reference,$alter)=(split /_/,$near_by_variant)[0,1,2,3];
		my $haplotype_for_nearby=$shapeit{$chromsome ."_".$position."_".$reference ."_". $alter};
		my $haplotype;
		if($state_of_haplotype ==1)
		{
		    $haplotype=$haplotype_for_nearby;
		}
		else
		{
		    $haplotype=reverse($haplotype_for_nearby);
		}
		my($chrom,$pos,$ref_nt,$alter_nt)=(split /_/,$variant)[0,1,2,3];
		my($first,$second)=(split /\|/,$haplotype)[0,1];
		print "$chrom . $pos $ref_nt $alter_nt $first $second\n";
	    }
	}
    }
}
close HAPCUT;
$/="\n";

sub nearest_variant
{
    my($variant,@vars)=@_;
    my($chr,$pos)=(split /_/,$variant)[0,1];
    my $distance="";
    my $near_var="";
    foreach my $temp ( @vars )
    {
	my($chrom,$position)=(split /_/,$temp)[0,1];
	my $tmp=abs($pos-$position);
	if( $distance eq "" )
	{
	    $distance = $tmp;
	    $near_var=$temp;
	    next;
	}
	else
	{
	    if($tmp < $distance)
	    {
		$distance=$tmp;
		$near_var=$temp;
		next;
	    }
	    else
	    {
		last;
	    }
	}
    }
    return $near_var;
}
