#!/usr/bin/perl
use warnings;
use Bio::DB::Fasta;

## fix some coordinate bugs in version 1 and 2
## fix haplotype information, only support 0|1 1|1 1|0

## usage:
my $input_file=$ARGV[0]; ##BED chr\tpos\tid\tref\talt\tscore\tpass\tother  (VCF format 10-col)
my $promoter_file=$ARGV[1];
my $fasta_library = 'data/hg19.fa';
my $database = Bio::DB::Fasta->new("$fasta_library") or die "Failed to creat Fasta DP object on fasta library\n";
my $path="script/TFM-Pvalue";
my $p_cut=4e-8;

my $pfm_file="data/motif.PFM";
# Read in PFM and transform it to PWM using 20 sequences.
my $A = 1;
my $C = 2;
my $G = 3;
my $T = 4;
my %motif;
my $prev_name_id;
open MOTIF, "$pfm_file" or die;
while(my $line = <MOTIF>)
{
    chomp($line);
    if($line =~ /^>/)
    {
	$prev_name_id = (split/>|\s+/,$line)[1];
    }
    else
    {
	my @info = split/\s+/,$line;
	if(not exists $motif{$prev_name_id})
	{
	    $motif{$prev_name_id}->[0] = {(A=>log(($info[$A]*20+0.25)/21)-log(0.25), T=>log(($info[$T]*20+0.25)/21)-log(0.25), C=>log(($info[$C]*20+0.25)/21)-log(0.25), G=>log(($info[$G]*20+0.25)/21)-log(0.25))};
	}
	else
	{
	    my $temp = $motif{$prev_name_id};
	    $motif{$prev_name_id}->[scalar(@$temp)] = {(A=>log(($info[$A]*20+0.25)/21)-log(0.25), T=>log(($info[$T]*20+0.25)/21)-log(0.25), C=>log(($info[$C]*20+0.25)/21)-log(0.25), G=>log(($info[$G]*20+0.25)/21)-log(0.25))};
	}
    }
}
close MOTIF;

my $score_file="data/motif.score.cut";
my %score_lower;
my %score;
my %score_upper;

open SCORE,"$score_file";
while(<SCORE>){
    my ($prev_name,$cut_off,$p) = (split /\s+/,$_)[0..2];
    if ($p < $p_cut){
	$score_upper{$prev_name} = $cut_off;
    }elsif ($p == $p_cut){
	$score{$prev_name} = $cut_off;
    }elsif($p > $p_cut){
	$score_lower{$prev_name} = $cut_off;
    }
}
close SCORE;

## statistic summary
my $number_of_input_variant=0;
my @variant_create_tfbs;
my @novel_tfbs;
my @transcript_affected;
my @multi_variant_tfbs;

open(VCF,"$input_file");
my @region;
while(<VCF>)
{
    chomp;
    next if /^#/;
    $number_of_input_variant++;
    my($chr,$pos,$ref,$alt)=(split /\t/,$_)[0,1,3,4];
    my $end=$pos+length($ref)-1;
    push @region,"$chr:$pos-$end";
}
my $region_s=join(" ",@region);
my @promoter=split /\n+/,`tabix $promoter_file $region_s |sort|uniq|intersectBed -a $input_file -b stdin -wa -wb|sort -k 2,3`;

my %haplotype;
foreach my $line ( @promoter )
{
    my($var_chr,$var_pos,$var_id,$var_ref,$var_alt,$var_dot,$var_pass,$a,$b,$hap,$promoter_chr,$promoter_start,$promoter_end,$promoter_gene,$promoter_strand)=split /\t/,$line;
    my $promoter_sympol=join("\t",$promoter_chr,$promoter_start,$promoter_end,$promoter_gene,$promoter_strand);
    my $variant_sympol=join("\t",$var_chr,$var_pos,$var_ref,$var_alt,$hap);
    if($hap=~/[123456789]\|/) #hap1
    {
	if(not exists $haplotype{$promoter_sympol}{'hap1'})
	{
	    $haplotype{$promoter_sympol}{'hap1'}=$variant_sympol;
	}
	else
	{
	    $haplotype{$promoter_sympol}{'hap1'} .= ";" . $variant_sympol;
	}
    }
    if($hap=~/\|[123456789]/) ## hap2
    {
	if(not exists $haplotype{$promoter_sympol}{'hap2'})
	{
	    $haplotype{$promoter_sympol}{'hap2'}=$variant_sympol;
	}
	else
	{
	    $haplotype{$promoter_sympol}{'hap2'} .= ";" . $variant_sympol;
	}
    }
}

my @key = keys %haplotype;
my %result;
foreach my $i ( @key )
{
    my($pro_chrom,$pro_start_site,$pro_end_site,$pro_name,$strand)=split /\t/,$i;
    my $length=$pro_end_site-$pro_start_site;
    foreach my $hap (keys %{$haplotype{$i}})
    {
	my $var= $haplotype{$i}{$hap};
	my @pos=split /;/,$var; ##each element is "chr pos ref alt"
	my($ref_tmp,$alt_tmp)=rebuild_sequence($pro_chrom,$pro_start_site,$pro_end_site,\@pos,$database);
	my @ref=@$ref_tmp;
	my @alt=@$alt_tmp;
	my $length_for_ref=@ref;
	foreach my $variation ( @pos )
	{
	    my($chr,$pos,$ref,$alt,$varhap)=split /\t/,$variation;
	    my $var_logo=$chr."_".$pos."_".$ref."/".$alt."_".$varhap;
	    my $first=$pos-1;## 0-based
	    my $end=$pos-1+length($ref)-1;
	    my $relative_start=$first-$pro_start_site; ## 0-based
	    my $relative_end=$end-$pro_start_site;
	    my $ref_seq_for_variant="";
	    my $alt_seq_for_variant="";
	    my @ref_seq_var;
	    my @alt_seq_var;
	    my($motif_name,$motif_seq,$score_motif, $novel_tfbs_start);
	    if($relative_start-29>=0 && $relative_end+29 <= $length-1)
	    {
		foreach my $ii ($relative_start-29..$relative_end+29)
		{
		    $ref_seq_for_variant.=$ref[$ii];
		    $alt_seq_for_variant.=$alt[$ii];
		}
		my $left_distance=29;
		#print "$ref_seq_for_variant\n$alt_seq_for_variant\n";
		($motif_name,$motif_seq,$score_motif,$novel_tfbs_start)=seq_scan($ref_seq_for_variant,$alt_seq_for_variant,\%motif,$strand,$left_distance,$ref,$alt);
		$novel_tfbs_start=$first-29+$novel_tfbs_start+1 if $novel_tfbs_start; ## the reference start site of novel tfbs 1 based
	    }
	    if($relative_start-29 < 0 && $relative_end+29 <= $length-1)
	    {
		foreach my $ii (0..$relative_end+29)
		{
		    $ref_seq_for_variant.=$ref[$ii];
		    $alt_seq_for_variant.=$alt[$ii];
		}
		my $left_distance=$relative_start;
		($motif_name,$motif_seq,$score_motif,$novel_tfbs_start)=seq_scan($ref_seq_for_variant,$alt_seq_for_variant,\%motif,$strand,$left_distance,$ref,$alt);
		$novel_tfbs_start=$pro_start_site+$novel_tfbs_start+1 if $novel_tfbs_start;
	    }
	    if($relative_start-29 >=0 && $relative_end+29 > $length-1)
	    {
		foreach my $ii ($relative_start-29..$length-1)
		{
		    $ref_seq_for_variant.=$ref[$ii];
		    $alt_seq_for_variant.=$alt[$ii];
		}
		my $left_distance=29;
		($motif_name,$motif_seq,$score_motif,$novel_tfbs_start)=seq_scan($ref_seq_for_variant,$alt_seq_for_variant,\%motif,$strand,$left_distance,$ref,$alt);
		$novel_tfbs_start=$first-29+$novel_tfbs_start+1 if $novel_tfbs_start;
	    }
	    if($relative_start-29 < 0 && $relative_end+29 > $length-1) ## it is impossible
	    {
	    }
	    if(defined($motif_seq))
	    {
		my $novel_tfbs_end=$novel_tfbs_start+length($motif_seq)-1; #1-based
		my $key=join("\t",$chr,$novel_tfbs_start,$novel_tfbs_end,$motif_name,$motif_seq,$score_motif,$pro_name,$strand);
		push @{$result{$key}},$var_logo;
		push @variant_create_tfbs,$var_logo unless $var_logo ~~@variant_create_tfbs;
		my $tmp1=$chr.":".$novel_tfbs_start.":".$novel_tfbs_end;
		push @novel_tfbs,$tmp1 unless $tmp1 ~~ @novel_tfbs;
		push @transcript_affected,$i unless $i ~~ @transcript_affected;
	    }
	    #print "$variation\t$motif_name\t$motif_seq\t$score_motif\t$i\n" if defined($motif_seq);
	}
    }
}
my @output;
foreach my $k (keys %result)
{
    if(scalar(@{$result{$k}})>=2)
    {
	my($a,$b,$c)=(split /\t/,$k)[0,1,2];
	my $tmp=$a.":".$b.":".$c;
	push @multi_variant_tfbs,$tmp unless $tmp ~~ @multi_variant_tfbs;
    }
    push @output,$k."\t".join(";",@{$result{$k}});
}
print "# Number of variants uploaded: ",$number_of_input_variant,"\n";
print "# Number of variants creating novel TFBSs: ",scalar(@variant_create_tfbs),"\n";
print "# Number of Novel TFBSs: ",scalar(@novel_tfbs),"\n";
print "# Number of multi-variant novel TFBSs: ",scalar(@multi_variant_tfbs),"\n";
print "# Number of transcripts affected: ",scalar(@transcript_affected),"\n\n";

print "#Chr\tStart\tEnd\tTF\tSequence\tPWM_score\tTranscript\tStrand\tVariant\n";
print join("\n",@output),"\n";


sub rebuild_sequence
{
    my($chr,$start,$end,$position,$reference_database)=@_;
    $start++;
    my @pos=@$position; ##format: ##each element is "chr pos ref alt" 
    my @pro_ref=split //,$reference_database->seq($chr,$start,$end);
    my @pro_alt=@pro_ref;
    foreach my $location (@pos)
    {
	my @info=split /\t/,$location;
	my $vstart=$info[1];
	my $vend=$info[1]+length($info[2])-1;
	$vstart= $vstart-$start; #to change the start to relative position
	$vend= $vend-$start; # to change the end position
	for my $i ($vstart..$vend)
	{
	    $pro_alt[$i]="";
	}
        $pro_alt[$vstart]= $info[3];
    }
    return(\@pro_ref,\@pro_alt);
}

sub seq_scan
{
    my ($ref_seq,$alt_seq,$motif_pwm,$strand,$left_distance,$ref,$alt)=@_;
    my %motif=%$motif_pwm;
    my @ref;
    my @alt;
    @ref=split //,$ref_seq;
    @alt=split //,$alt_seq;
    my $length_ref=@ref;
    my $length_alt=@alt;
    undef my %refpwm;
    my $prev_name;
    foreach $prev_name(sort keys %motif)
    {
	#print "start $prev_name\n";
	my $motif_length = scalar @{$motif{$prev_name}};
	unless (-d "pwm")
	{
	    mkdir "pwm";
	}
	if (-e "pwm/$prev_name")
	{
	}
	else
	{
	    open(TMP,">pwm/$prev_name")||die;
	    my $pwm = "";
	    for $i(0 .. $motif_length-2)
	    {
		$pwm .= join("",$motif{$prev_name}->[$i]->{"A"},"\t");
	    }
	    $pwm .= join("",$motif{$prev_name}->[$motif_length-1]->{"A"},"\n");
	    for $i(0 .. $motif_length-2)
	    {
		$pwm .= join("",$motif{$prev_name}->[$i]->{"C"},"\t");
	    }
	    $pwm .= join("",$motif{$prev_name}->[$motif_length-1]->{"C"},"\n");
	    for $i(0 .. $motif_length-2)
	    {
		$pwm .= join("",$motif{$prev_name}->[$i]->{"G"},"\t");
	    }
	    $pwm .= join("",$motif{$prev_name}->[$motif_length-1]->{"G"},"\n");
	    for $i(0 .. $motif_length-2)
	    {
		$pwm .= join("",$motif{$prev_name}->[$i]->{"T"},"\t");
	    }
	    $pwm .= join("",$motif{$prev_name}->[$motif_length-1]->{"T"},"\n");
	    print TMP $pwm;
	    close TMP;
	}

	#print "$prev_name!\n";# sequence scanning .... reference
	my $scan_begain;
	my $scan_end;
	$scan_end=$left_distance + length($ref) if $left_distance + length($ref) <= $length_ref-$motif_length;
	$scan_end=$length_ref-$motif_length if $length_ref-$motif_length < $left_distance + length($ref);
	$scan_begain=$left_distance+1-$motif_length if $left_distance+1-$motif_length>=0;
	$scan_begain=0 if $left_distance+1-$motif_length < 0;
	for (my $i=$scan_begain; $i < $scan_end; $i ++)
	{
	    my $ref_pwm_score = 0;
	    # calculate reference  PWM score;
	    my $ref_seq_motif="";
	    for (my $j = 0; $j < $motif_length; $j++){
		if($strand eq "+")
		{
		    print $i+$j,"\n",$i,"\t",$j,"\n" unless($ref[$i+$j]);

		    $ref_pwm_score += $motif{$prev_name}->[$j]->{$ref[$i+$j]};
		    $ref_seq_motif.=$ref[$i+$j];
		}
		else
		{
		    my $nt=$ref[$i+$j];
		    $nt=~tr/ATGC/TACG/;
		    $ref_pwm_score += $motif{$prev_name}->[$motif_length-$j-1]->{$nt};
		    $ref_seq_motif=$nt.$ref_seq_motif;
		}
	    }
	    #print "$ref_seq_motif\t$ref_pwm_score\n" if $ref_pwm_score >10;
	    if (defined $score_upper{$prev_name} && $ref_pwm_score >= $score_upper{$prev_name})
	    {
		$refpwm{$prev_name}=1;
		#print $i,"\n";
	    }
	    elsif(defined $score_lower{$prev_name} && $ref_pwm_score < $score_lower{$prev_name})
	    {
		next;
	    }
	    elsif(defined $score{$prev_name} && $ref_pwm_score >= $score{$prev_name})
	    {
		$refpwm{$prev_name}=1;
		#print $i,"\n";
	    }
	    elsif(defined $score{$prev_name} && $ref_pwm_score < $score{$prev_name})
	    {
		next;
	    }
	    else
	    {
		#print "aaaa\n";
		my $ref_p = (split /\s+/,`$path/TFMpvalue-sc2pv -a 0.25 -c 0.25 -g 0.25 -t 0.25 -s $ref_pwm_score -m pwm/$prev_name -w`)[2];
		#print "chengsj\n";
		if ($ref_p < $p_cut)
		{
		    $refpwm{$prev_name}=1;
		}
	    }
		#print "here!\n";
	}
	#print "finished $prev_name!\n";
    }
	
    #print keys %refpwm,"\n";

    foreach $prev_name(sort keys %motif)
    {
	if (defined $refpwm{$prev_name})
	{
	}
	else
	{
	    $motif_length = scalar @{$motif{$prev_name}};

	    # sequence scanning....  alternative
	    my $final_position; ## output the start of the novel TFBS 0-based
	    my $scan_start;
	    $scan_start=$left_distance+1 - $motif_length if $left_distance+1-$motif_length>=0;
	    $scan_start=0 if $left_distance+1-$motif_length <0;
	    my $scan_end;
	    $scan_end=$left_distance + length($alt) if $left_distance + length($alt) <= $length_alt-$motif_length;
	    $scan_end=$length_alt-$motif_length if $length_alt-$motif_length < $left_distance + length($alt);
	    
	    for (my $i=$scan_start; $i < $scan_end; $i ++)
	    {
		my $alt_pwm_score = 0;
		my $alt_motif_seq='';
		$final_position=$i;
		# calculate reference & alternative PWM score;
		for (my $j = 0; $j < $motif_length; $j++)
		{
		    if($strand eq '+')
		    {
			$alt_pwm_score += $motif{$prev_name}->[$j]->{$alt[$i+$j]};
			$alt_motif_seq .=$alt[$i+$j];
		    }
		    else
		    {
			my $nt=$alt[$i+$j];
			$nt=~tr/ATGC/TACG/;
			$alt_pwm_score += $motif{$prev_name}->[$motif_length-$j-1]->{$nt};
			$alt_motif_seq=$nt .$alt_motif_seq;
		    }
		}
		if (defined $score_upper{$prev_name} && $alt_pwm_score >= $score_upper{$prev_name})
		{
		    return($prev_name,$alt_motif_seq,$alt_pwm_score,$final_position);
		}
		elsif(defined $score_lower{$prev_name} && $alt_pwm_score < $score_lower{$prev_name}){
		}
		elsif(defined $score{$prev_name} && $alt_pwm_score>= $score{$prev_name})
		{
			return($prev_name,$alt_motif_seq,$alt_pwm_score,$final_position);
		}
		else{
		    my $alt_p = (split /\s+/,`$path/TFMpvalue-sc2pv -a 0.25 -c 0.25 -g 0.25 -t 0.25 -s $alt_pwm_score -m pwm/$prev_name -w`)[2];
		    if ($alt_p < $p_cut){
			return($prev_name,$alt_motif_seq,$alt_pwm_score,$final_position);
		    }
		}
	    }
	}
    }
}

