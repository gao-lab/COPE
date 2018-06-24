#!/usr/bin/perl
use warnings;
use Bio::DB::Fasta;

## input file: variants(VCF,10-tab) only support 0|1 1|0 1|1
my $input_file=$ARGV[0];
my $bound_motif=$ARGV[1];
#my $bound_motif="data/ENCODE.tf.bound.union.bed"; ##7 col
my $motif_pfm="data/motif.PFM";

##reading in PFM
my %motif;
my $A = 1;
my $C = 2;
my $G = 3;
my $T = 4;
open MOTIF,$motif_pfm or die;
#my %motif;
my $motif_name_id;
while(<MOTIF>)
{
    chomp $_;
    if(/^>/){
	$motif_name_id = (split/>|\s+/,$_)[1];
    }else{
	my @info = split/\s+/,$_;
	if(not exists $motif{$motif_name_id}){
	    $motif{$motif_name_id}->[0] = {(A=>$info[$A], T=>$info[$T], C=>$info[$C], G=>$info[$G])};
	}else{
	    my $temp = $motif{$motif_name_id};
	    $motif{$motif_name_id}->[scalar(@$temp)] = {(A=>$info[$A], T=>$info[$T], C=>$info[$C], G=>$info[$G])};
	}
    }
}
close MOTIF;

## statistic summary
my $number_of_input_variant=0;
my @variant_within_tfbs;
my @tfbs_affected;
my @multi_variant_tfbs;
my @output;


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
my $tfbs= `tabix $bound_motif  $region_s |sort|uniq|intersectBed -a stdin -b $input_file -wa -wb`;
my @line = split /\n+/, `printf "$tfbs"|gawk -F '\t' '{print \$1"\t"\$2"\t"\$3"\t"\$4"\t"\$6"\t"\$8"\t"\$9"\t"\$11"\t"\$12"\t"\$17}'`;
my %tfbs_variant;
foreach my $line (@line)
{
    my($chr,$tfbs_start,$tfbs_end,$tfbs_name,$strand,$chrom,$var_pos,$ref,$alt,$hap)=split /\t/,$line;
    my $var_end=$var_pos+length($ref)-1;
    my $tfbs_sympol=join("",$chr,"\t",$tfbs_start,"\t",$tfbs_end,"\t",$tfbs_name,"\t",$strand);
    my $variant_sympol=join("",$chrom,"\t",$var_pos,"\t",$var_end,"\t",$ref,"\t",$alt,"\t",$hap);
    if($hap=~/[123456789]\|/) #hap1
    {
	if(not exists $tfbs_variant{$tfbs_sympol}{'hap1'})
	{
	    $tfbs_variant{$tfbs_sympol}{'hap1'}=$variant_sympol;
	}
	else
	{
	    $tfbs_variant{$tfbs_sympol}{'hap1'} .= ";" .$variant_sympol;
	}
    }
    if($hap=~/\|[123456789]/) ## hap2
    {
	if(not exists $tfbs_variant{$tfbs_sympol}{'hap2'})
	{
	    $tfbs_variant{$tfbs_sympol}{'hap2'}=$variant_sympol;
	}
	else
	{
	    $tfbs_variant{$tfbs_sympol}{'hap2'} .= ";" .$variant_sympol;
	}
    }
}

##rebuild tfbs sequence
my $fasta_library = 'data/hg19.fa';
my $database = Bio::DB::Fasta->new("$fasta_library") or die "Failed to creat Fasta DP object on fasta library\n";
foreach $key (keys %tfbs_variant)
{
    my($chr_tfbs,$start_tfbs,$end_tfbs,$id_tfbs,$strand)=split /\t/,$key;
    #my($chr_var,$start_var,$end_var,$ref_var,$alt_var)=split /\t/,$tfbs_variant{$key};
    $start_tfbs++;
    foreach my $hap (keys %{$tfbs_variant{$key}})
    {
	my @pos=split /;/,$tfbs_variant{$key}{$hap};
	my $ref_seq=$database->seq($chr_tfbs, $start_tfbs, $end_tfbs);
	my @alt_seq=(split //,$ref_seq);
	my @variant_tag;
	foreach my $postion (@pos)
	{
	    my @info=split /\t/,$postion;
	    my $tmp=$info[0]."_".$info[1]."_".$info[3]."/".$info[4]."_".$info[5];
	    push @variant_tag,$tmp;
	    $info[1]= $info[1]-$start_tfbs+1; #to change the start to relative position
	    $info[2]= $info[2]-$start_tfbs+1; # to change the end position
	    $info[3]="" if $info[3] eq "-";
	    $info[4]="" if $info[4] eq "-";
	    #$alt_seq=$database->subseq($id,1,$info[0]-1). $info[3] . $database->subseq($id,$info[1]+1,$len) if $info[0]>1;
	    if($info[1]>0)
	    {
		if($info[3] eq "")
		{
		    $alt_seq[$info[1]-1].= $info[4];
		}
		else{$alt_seq[$info[1]-1]= $info[4];}
		for my $i ($info[1]..$info[2]-1)
		{
		    if($alt_seq[$i])
		    {
			$alt_seq[$i]="";
		    }
		}
	    }
	    else
	    {
		$info[1]=1;
		$alt_seq[$info[1]-1]= $info[4];
		for my $i ($info[1]..$info[2]-1)
		{
		    if($alt_seq[$i])
		    {
			$alt_seq[$i]="";
		    }
		}
	    }
	}
	my @ref_final_seq;
	my @alt_final_seq;
	if($strand eq "-")
	{
	    $ref_seq=~ tr/ACGTacgt/TGCAtgca/;
	    @ref_final_seq=reverse(split //,$ref_seq);
	    my @tmp=reverse(@alt_seq);
	    my $temp=join "",@tmp;
	    $temp=~ tr/ACGTacgt/TGCAtgca/;
	    @alt_final_seq=split //,$temp;
	}
	else
	{
	    @ref_final_seq=split //,$ref_seq;
	    my $temp=join("",@alt_seq);
	    @alt_final_seq=split //,$temp;
	}
	##calculate score
	my $ref_length=@ref_final_seq;
	my $alt_length=@alt_final_seq;
	my $id=(split /_\dmer/,$id_tfbs)[0];
	if(not exists $motif{$id}){
	    #print "Motif_Name_Not_Found\t$line\n";
	    next;
	}
	else
	{
	    my $ref_score=0;
	    my $alt_score=0;
	    foreach my $b (0..$ref_length-1)
	    {
		$ref_score += $motif{$id}->[$b]->{$ref_final_seq[$b]};
	    }
	    if($ref_length == $alt_length)
	    {
		#my $ref_score=0;
		#my $alt_score=0;
		foreach my $a (0..$ref_length-1)
		{
		    #$ref_score += $motif{$id}->[$a]->{$ref_final_seq[$a]};
		    $alt_score += $motif{$id}->[$a]->{$alt_final_seq[$a]};
		}
	    }
	    elsif($ref_length < $alt_length)
	    {
		$alt_score=max_matrix_score(\@alt_final_seq,\%motif,$ref_length,$id);
	    }
	    else
	    {
		my $l=$ref_length-$alt_length;
		my $left=$database->seq($chr_tfbs, $start_tfbs-$l, $start_tfbs-1);
		my $right=$database->seq($chr_tfbs, $end_tfbs+1, $end_tfbs+$l);
		my @change_alt_seq;
		if($strand eq "+")
		{
		    push @change_alt_seq,$left;
		    push @change_alt_seq,@alt_final_seq;
		    push @change_alt_seq,$right;
		}
		else
		{
		    $left=~ tr/ATCGatcg/TAGCTAGC/;
		    $right=~ tr/ATCGatcg/TAGCTAGC/;
		    push @change_alt_seq,reverse($right);
		    push @change_alt_seq,@alt_final_seq;
		    push @change_alt_seq,reverse($left);
		}
		my @final_alt_sequence=split //,join("",@change_alt_seq);
		$alt_score=max_matrix_score(\@final_alt_sequence,\%motif,$ref_length,$id);
	    }
	    if($ref_score>$alt_score) ##TFBS breaking
	    {
		#print $chr_tfbs,"\t", $start_tfbs,"\t",$end_tfbs,"\t",$id,"\n";
		my $tmp2=$chr_tfbs."\t". $start_tfbs."\t".$end_tfbs."\t".$id;
		push @tfbs_affected,$tmp2 unless $tmp2 ~~ @tfbs_affected;
		if(scalar(@variant_tag)>=2)
		{
		    push @multi_variant_tfbs,$tmp2 unless $tmp2 ~~@multi_variant_tfbs;
		}
		foreach my $a (@variant_tag)
		{
		    push @variant_within_tfbs,$a unless $a ~~ @variant_within_tfbs;
		}
		 my $tmp1=$chr_tfbs."\t".$start_tfbs."\t".$end_tfbs."\t".$id."\t".$ref_score."\t".$alt_score."\t".$strand."\t".join(";",@variant_tag);
		 push @output,$tmp1;
	     }
	}
    }
}

print "# Number of variants uploaded: ",$number_of_input_variant,"\n";
print "# Number of variants within TFBSs: ",scalar(@variant_within_tfbs),"\n";
print "# Number of TFBSs affected: ",scalar(@tfbs_affected),"\n";
print "# Number of TFBSs with multi-variants: ",scalar(@multi_variant_tfbs),"\n\n";

print "#Chr\tStart\tEnd\tTF\tRef_PWM_score\tAlt_PWM_score\tStrand\tVariant\n";
print join("\n",@output),"\n";


sub max_matrix_score
{
    my($sequence,$motif,$motif_length,$prev_name)=@_;
    my @seq=@$sequence;
    my %motif_pfm=%$motif;
    my $matrix_score=0;
    my $max_score=0;
    my $length=@seq-$motif_length;
    foreach my $i (0..$length)
    {
	for ($j = 0; $j < $motif_length; $j++)
	{
	    $matrix_score += $motif{$prev_name}->[$j]->{$seq[$i+$j]};
	}
	    if( $matrix_score >$max_score)
	    {
		$max_score=$matrix_score;
	    }
    }
    return($max_score);
}
