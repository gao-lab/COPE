#!/usr/bin/perl

# extract_novel_transcript_and_sequence_alignment.pl 
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



use warnings;
use Set::IntervalTree;
use Bio::DB::Fasta;
use Parallel::ForkManager;
use Getopt::Long;

##useage
my $usage="perl find_new_transcript.pl -n needle -aminoref protein.fa -r ref.fa -v variant.vcf -g gene.bed -p int -c 100 -splicegain [0|1] -bed [0|1] -score [0|1] -confident [0|1]";

##input
my $reference_genome;
my $variant_file; ##VCF format without "chr"
my $gene_model_file=$ARGV[1]; ##bed format
my $max_process;
my $context;
my $splicegain;
my $needle_path;
my $aminoseq;
my $bed;
my $score;
my $confident;

GetOptions ("needle|n=s" => \$needle_path,
	    "aminoref|a=s" => \$aminoseq,
	    "reference|r=s" => \$reference_genome,
	    "variant|v=s" => \$variant_file,
	    "gene|g=s"  => \$gene_model_file,
	    "process|p=i"  => \$max_process,
	    "context|c=i" => \$context,
	    "splicegain|s=i" => \$splicegain,
	    "bed|b=i" => \$bed,
	    "score=i" => \$score,
	    "confident=i" => \$confident,
	   );
##global variable
unless($variant_file && $gene_model_file && $max_process && $context )
{
    die "$usage\n";
}
my %trans_id;
open(GE,"./resource/Refseq2Gene.NM.prefix.txt") or die;
while(<GE>)
{
    chomp;
    my($trans,$genename)=split /\t/,$_;
    $trans_id{$trans}=$genename;
}

my $pm = new Parallel::ForkManager($max_process);
my @total_seq;
$pm -> run_on_finish( sub {
			  my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data_structure_reference) = @_;
			  if($data_structure_reference)
			  {
			      push @total_seq,@$data_structure_reference;
			  }
		      }
		    );

## config for maxentscan
##3ass
my $dir="./resource/maxentscan/splicemodels/";
my @list = ('me2x3acc1','me2x3acc2','me2x3acc3','me2x3acc4','me2x3acc5','me2x3acc6','me2x3acc7','me2x3acc8','me2x3acc9');
my @metables; ##used for score 3UTR
my $num = 0 ;
foreach my $file (@list)
{
    my $n = 0;
    open (SCOREF,"<".$dir.$file) || die "Can't open $file!\n";
    while(<SCOREF>)
    {
	chomp;
	$_=~ s/\s//;
	$metables[$num]{$n} = $_;
	$n++;
    }
    close(SCOREF);
    #print STDERR $file."\t".$num."\t".$n."\n";
    $num++;
}
##5dss
my %me2x5 = dssmakescorematrix('./resource/maxentscan/me2x5');
my %seq = dssmakesequencematrix('./resource/maxentscan/splicemodels/splice5sequences');
sub dssmakesequencematrix{
    my $file = shift;
    my %matrix;my $n=0;
    open(SCOREF, $file) || die "Can't open $file!\n";
    while(<SCOREF>) {
	chomp;
	$_=~ s/\s//;
	$matrix{$_} = $n;
	$n++;
    }
    close(SCOREF);
    return %matrix;
}
sub dssmakescorematrix{
    my $file = shift;
    my %matrix;my $n=0;
    open(SCOREF, $file) || die "Can't open $file!\n";
    while(<SCOREF>) {
	chomp;
	$_=~ s/\s//;
	$matrix{$n} = $_;
	$n++;
    }
    close(SCOREF);
    return %matrix;
}

my $fasta_library = "./resource/$reference_genome";
my $pro_lib="./resource/$aminoseq";
my $pro_database = Bio::DB::Fasta->new("$pro_lib");
my $database = Bio::DB::Fasta->new("$fasta_library") or die "Failed to creat Fasta DP object on fasta library\n";
## splicing branch
open(PFM,"./resource/maxentscan/splicing_branch.pfm") or die "cannot open the file!";
my %motif;
my $A = 1;
my $C = 2;
my $G = 3;
my $T = 4;
while(my $line=<PFM>)
{
    chomp $line;
    next if $line=~/^>/;
    my $name='branch';
    my @info = split/\s+/,$line;
    if(not exists $motif{$name})
    {
	$motif{$name}->[0] = {(A=>log($info[$A])-log(0.25), T=>log($info[$T])-log(0.25), C=>log($info[$C])-log(0.25), G=>log($info[$G])-log(0.25))};
    }
    else
    {
	my $temp = $motif{$name};
	$motif{$name}->[scalar(@$temp)]={(A=>log($info[$A])-log(0.25), T=>log($info[$T])-log(0.25), C=>log($info[$C])-log(0.25), G=>log($info[$G])-log(0.25))};
    }
}
close PFM;


##loading variant with phasing result regardless of variants that without phasing result
my $tree = Set::IntervalTree->new;
open(VAR,$variant_file) or die "cannot open the file!";
while(<VAR>)
{
    next if /^#/;
    chomp;
    my($chr,$start,$ref,$alt,$haplotype)=(split /\t/,$_)[0,1,3,4,9];
    my $end = $start+length($ref)-1;
    my $id= $chr .":".$start.":".$end.":".$ref.":".$alt.":".$haplotype;
    if(length($ref) > length($alt))
    {
	$tree->insert($id,$start+1,$end+1);##like bed format [);
    }
    else
    {
	$tree->insert($id,$start,$start+1);
    }
}
close VAR;

##deal with each gene
open(GENEMODEL,$gene_model_file) or die "cannot open the gene model file";
while(my $line=<GENEMODEL>)
{
    $pm->start and next;
    my @novel_seq;
    chomp $line;
    my($chrom,$start,$end,$gene_id,$strand,$cds_start,$cds_end,$block_length,$block_start)=(split /\t/,$line)[0,1,2,3,5,6,7,10,11];
    my $GeneName="";
    my @splice_score;
    my @confident_score; ## motif score of new splice site
    if(defined($trans_id{$gene_id}))
    {
	$GeneName=$trans_id{$gene_id};
    }
    my $t_start=$start;
    my $t_end=$end;
    my @length=split /,/,$block_length;
    my @start=split /,/,$block_start;
    my $overlaps = $tree->fetch($start+1,$end+1); #like bed format [);
    my @variant_in_gene;
    if(scalar(@$overlaps)>0)
    {
	foreach my $region (@$overlaps)
	{
	    my $r_chr=(split /:/,$region)[0];
	    next if ($chrom ne $r_chr);  # ignore if overlap is on different chromosome
	    push @variant_in_gene,$region;
	}
	if(@variant_in_gene)
	{
	    my @zero_one_gene;
	    my @one_zero_gene;
	    foreach my $var ( @variant_in_gene )
	    {
		push @zero_one_gene,$var if $var =~/0\|1/ || $var =~/1\|1/;
		push @one_zero_gene,$var if $var =~/1\|0/ || $var =~/1\|1/;
	    }
	    if(@zero_one_gene)
	    {
		my $splice_state=0;
		my $state_for_splice_donor=0;
		my $state_for_splice_acceptor=0;## 0 strands for unbroken
		my($wild_gene_seq,@alter_seq_array)=reconstruct_seq($database,$chrom,$start+1,$end,$strand,@zero_one_gene);
		$pm->finish(0) if $wild_gene_seq=~/N/;
		my $alter_gene_seq=join("",@alter_seq_array);
		my @dss_zero_one;
		my @ass_zero_one;
		my $relative_start;
		if($strand eq "+")
		{
		    my $start_codon_altered=relative_position_altered($cds_start+1,@zero_one_gene);
		    $relative_start=$cds_start+1-$start+$start_codon_altered;
		    my $start_codon_first=$cds_start+1;
		    my @dss;
		    my @ass;
		    my $n=@length;
		    foreach my $i (0..$n-2)
		    {
			push @dss,$start+$start[$i]+$length[$i]+1; ##1-based
			push @ass,$start+$start[$i+1];##1-based
		    }
		    foreach my $j (0..$n-2)
		    {
			##for dss
			my $dss_overlaps = $tree->fetch($dss[$j]-3,$dss[$j]+5+1);
			my @dss_variant;
			if(scalar(@$dss_overlaps)>0)
			{
			    foreach my $region (@$dss_overlaps)
			    {
				my $r_chr=(split /:/,$region)[0];
				next if ($chrom ne $r_chr);  # ignore if overlap is on different chromosome
				push @dss_variant,$region;
			    }
			}
			##reconstruct splicing consense region to predict splicing-altering variants
			my $zero_one_state_of_dss=0;
			if(@dss_variant)
			{
			    my (@zero_one);
			    foreach my $var ( @dss_variant )
			    {
				push @zero_one,$var if $var =~/0\|1/ ||$var =~/1\|1/;
				#push @one_zero,$var if $var =~/1\|0/ ||$var =~/1\|1/;
			    }
			    if(@zero_one)
			    {
				my $dss_changed_bases=relative_position_altered($dss[$j],@zero_one_gene);
				my $wild_seq=substr($wild_gene_seq,$dss[$j]-3-$start-1,9);
				my $alter_seq=substr($alter_gene_seq,$dss[$j]-3+$dss_changed_bases-$start-1,9);
				my $wild_score=score5($wild_seq);
				my $alter_score=score5($alter_seq);
				$wild_score=$wild_score+0.001;
				my $final_score=($wild_score-$alter_score)/$wild_score;##to be continued;
				#my $zero_one_state;
				if($wild_score >0 && $final_score >= 0.217)
				{
				    $zero_one_state_of_dss=1;
				    push @splice_score,$final_score;
				}
				elsif($wild_score >0 && $final_score < 0.217)
				{
				    $zero_one_state_of_dss=0;
				}
				else
				{
				    if($wild_score>$alter_score)
				    {
					$zero_one_state_of_dss=1;##unknow if as 1, we don't add this constraint
					#push @splice_score,$;
				    }
				    else
				    {
					$zero_one_state_of_dss=0;
				    }
				}
			    }
			}
			if($zero_one_state_of_dss==0)
			{
			    my $alter_base=relative_position_altered($dss[$j],@zero_one_gene);
			    my $tmp=($dss[$j]+$alter_base-$start);
				#print $tmp;
				#print join(" ",@dss_zero_one);
			    push @dss_zero_one,$tmp;
			}
			else #search for cryptic dss site
			{
			    $splice_state=1;
			    my $alter_base=relative_position_altered($dss[$j],@zero_one_gene);
			    my $upper=0;
			    my $bottom=0;
			    if($j==0)
			    {
				$upper=$dss[$j]-$start-$context if $dss[$j]-$start-$context>0;
				$upper=1 if $dss[$j]-$start-$context <= 0;
			    }
			    else
			    {
				if($dss[$j]-$ass[$j-1]>$context)
				{
				    $upper=$dss[$j]-$context-$start;
				}
				else
				{
				    $upper=$ass[$j-1]-$start+1;
				}
			    }
			    if($ass[$j]-$dss[$j]>$context)
			    {
				$bottom=$dss[$j]+$context-$start;
			    }
			    else
			    {
				$bottom=$ass[$j]-1-$start;
			    }
			    my($site,$confident)=search_cryptic_dss($dss[$j]-$start+$alter_base,$alter_gene_seq,$upper+$alter_base,$bottom+$alter_base,$strand);
			    if($site==0) ## no cryptic splice site
			    {
				#$state_for_splice_donor=1;##means broken original structure
				## search for new splice dss
				my @site_of_dss;
				my $specific_start;
				my $specific_end;
				if($j==0)
				{
				    $specific_start=$start+1;
				}
				else
				{
				    $specific_start=$ass[$j-1];
				}
				$specific_end=$ass[$j];
				foreach my $variant (@zero_one_gene)
				{
				    my ($pos_first,$pos_second,$allele_alt)=(split /:/,$variant)[1,2,4];
				    if($specific_start<$pos_first && $pos_second<$specific_end) ## in specific region
				    {
					my $alter_base=relative_position_altered($pos_first,@zero_one_gene);
					my $real_position=$pos_first+$alter_base-$start;
					my $ref_pos=$pos_first-$start;
					my @new_dss=check_dss($alter_gene_seq,$wild_gene_seq,$real_position,$ref_pos);
					push @site_of_dss,@new_dss;
				    }
				}
				if(@site_of_dss)
				{
				    push @dss_zero_one,@site_of_dss;
				}
				else
				{
				    $state_for_splice_donor=1; ## ##means broken original structure
				}
			    }
			    else ## find cryptic splice site
			    {
				push @dss_zero_one,$site;
				push @confident_score,$confident;
			    }
			}
			####for ass
			my $ass_overlaps = $tree->fetch($ass[$j]-19,$ass[$j]+3+1);
			my @ass_variant;
			if(scalar(@$ass_overlaps)>0)
			{
			    foreach my $region (@$ass_overlaps)
			    {
				my $r_chr=(split /:/,$region)[0];
				next if ($chrom ne $r_chr);  # ignore if overlap is on different chromosome
				push @ass_variant,$region;
			    }
			}
			##reconstruct splicing consense region to predict splicing-altering variants
			my $zero_one_state_of_ass=0;
			if(@ass_variant )
			{
			    my (@zero_one,@one_zero);
			    foreach my $var ( @ass_variant )
			    {
				push @zero_one,$var if $var =~/0\|1/ ||$var =~/1\|1/;
			    }
			    #my $zero_one_state="unbroken";
			    if(@zero_one )
			    {
				my $ass_changed_bases=relative_position_altered($ass[$j],@zero_one_gene);
				my $wild_seq=substr($wild_gene_seq,$ass[$j]-19-$start-1,23);
				my $alter_seq=substr($alter_gene_seq,$ass[$j]-19+$ass_changed_bases-$start-1,23);
				my $wild_score=score3($wild_seq);
				my $alter_score=score3($alter_seq);
				$wild_score=$wild_score+0.001;
				my $final_score=($wild_score-$alter_score)/$wild_score;##to be continued;
				#my $zero_one_state;
				if($wild_score >0 && $final_score >= 0.217)
				{
				    $zero_one_state_of_ass=1;
				    push @splice_score,$final_score;
				}
				elsif($wild_score >0 && $final_score < 0.217)
				{
				    $zero_one_state_of_ass=0;
				}
				else
				{
				    if($wild_score>$alter_score)
				    {
					$zero_one_state_of_ass=1;##unknow if as 1, we don't add this constraint
				    }
				    else
				    {
					$zero_one_state_of_ass=0;
				    }
				}
			    }
			}
			if($zero_one_state_of_ass ==0)
			{
			    my $alter_base=relative_position_altered($ass[$j],@zero_one_gene);
			    my $tmp=($ass[$j]+$alter_base-$start);
				#print $tmp;
			    push @ass_zero_one,$tmp;
			}
			else
			{
			    $splice_state=1;
			    my $alter_base=relative_position_altered($ass[$j],@zero_one_gene);
			    my $upper=0;
			    my $bottom=0;
			    if($ass[$j]-$dss[$j]>$context)
			    {
				$upper=$ass[$j]-$context-$start;
			    }
			    else
			    {
				$upper=$dss[$j]+1-$start;
			    }
			    if(defined($dss[$j+1]))
			    {
				if($dss[$j+1]-$ass[$j]>$context)
				{
				    $bottom=$ass[$j]+$context-$start;
				}
				else
				{
				    $bottom=$dss[$j+1]-1-$start;
				}
			    }
			    else
			    {
				$bottom=$ass[$j]-$start+$context;
			    }
			    my($site,$confident)=search_cryptic_ass($ass[$j]-$start+$alter_base,$alter_gene_seq,$upper+$alter_base,$bottom+$alter_base,$strand);
			    if($site==0) ## without splice site
			    {
				#$state_for_splice_acceptor=1;
				my @site_of_ass;
				my $specific_start;
				my $specific_end;
				if(defined($dss[$j+1]))
				{
				    $specific_end=$dss[$j+1];
				}
				else
				{
				    $specific_end=$end;
				}
				$specific_start=$dss[$j];
				foreach my $variant (@zero_one_gene)
				{
				    my ($pos_first,$pos_second,$allele_alt)=(split /:/,$variant)[1,2,4];
				    if($specific_start<$pos_first && $pos_second<$specific_end) ## in specific region
				    {
					my $alter_base=relative_position_altered($pos_first,@zero_one_gene);
					my $real_position=$pos_first+$alter_base-$start;
					my $ref_pos=$pos_first-$start;
					my @new_ass=check_ass($alter_gene_seq,$wild_gene_seq,$real_position,$ref_pos);
					push @site_of_ass,@new_ass;
				    }
				}
				if(@site_of_ass)
				{
				    push @ass_zero_one,@site_of_ass;
				}
				else
				{
				    $state_for_splice_acceptor=1; ## means broken original structure
				}
			    }
			    else # with splice site
			    {
				push @ass_zero_one,$site;
				push @confident_score,$confident;
			    }
			}
		    }
		}
		else ##negative strand
		{
		    my $start_codon_altered=relative_position_altered($cds_end,@zero_one_gene);
		    $relative_start=$cds_end-$start+$start_codon_altered;
		    my @dss;
		    my @ass;
		    my $n=@length;
		    foreach my $i (0..$n-2)
		    {
			push @ass,$start+$start[$i]+$length[$i]+1; ##1-based ##need to change
			push @dss,$start+$start[$i+1];##1-based
		    }
		    foreach my $j (0..$n-2)
		    {
			####for ass
			my $ass_overlaps = $tree->fetch($ass[$j]-3,$ass[$j]+19+1);
			my @ass_variant;
			if(scalar(@$ass_overlaps)>0)
			{
			    foreach my $region (@$ass_overlaps)
			    {
				my $r_chr=(split /:/,$region)[0];
				next if ($chrom ne $r_chr);  # ignore if overlap is on different chromosome
				push @ass_variant,$region;
			    }
			}
			##reconstruct splicing consense region to predict splicing-altering variants
			my $zero_one_state_of_ass=0;
			if(@ass_variant)
			{
			    my @zero_one;
			    foreach my $var ( @ass_variant )
			    {
				push @zero_one,$var if $var =~/0\|1/ ||$var =~/1\|1/;
				#push @one_zero,$var if $var =~/1\|0/ ||$var =~/1\|1/;
			    }
			    #my $zero_one_state="unbroken";
			    if(@zero_one )
			    {
				my $ass_changed_bases=relative_position_altered($ass[$j],@zero_one_gene);
				my $wild=substr($wild_gene_seq,$ass[$j]-3-$start-1,23);
				my $alter=substr($alter_gene_seq,$ass[$j]-3-$start+$ass_changed_bases-1,23);
				my $wild_seq=reverse(to_complement($wild));
				my $alter_seq=reverse(to_complement($alter));
				my $wild_score=score3($wild_seq);
				my $alter_score=score3($alter_seq);
				$wild_score=$wild_score+0.001;
				my $final_score=($wild_score-$alter_score)/$wild_score;##to be continued;
				#my $zero_one_state;
				if($wild_score >0 && $final_score >= 0.217)
				{
				    $zero_one_state_of_ass=1;
				}
				elsif($wild_score >0 && $final_score < 0.217)
				{
				    $zero_one_state_of_ass=0;
				}
				else
				{
				    if($wild_score>$alter_score)
				    {
					$zero_one_state_of_ass=1;##unknow if as 1, we don't add this constraint
					push @splice_score,$final_score;
				    }
				    else
				    {
					$zero_one_state_of_ass=0;
				    }
				}
			    }
			}
			if($zero_one_state_of_ass==0)
			{
			    my $alter_base=relative_position_altered($ass[$j],@zero_one_gene);
			    my $tmp=$ass[$j]+$alter_base-$start;
			    push @ass_zero_one,$tmp;
			}
			else #search for cryptic ass
			{
			    $splice_state=1;
			    my $alter_base=relative_position_altered($ass[$j],@zero_one_gene);
			    my $upper=0;
			    my $bottom=0;
			    if($j==0)
			    {
				if($ass[$j]-$start-$context<1)
				{
				    $upper=1;
				}
				else
				{
				    $upper=$ass[$j]-$start-$context;
				}
			    }
			    else
			    {
				if($ass[$j]-$dss[$j-1]>$context)
				{
				    $upper=$ass[$j]-$context-$start;
				}
				else
				{
				    $upper=$dss[$j-1]+1-$start;
				}
			    }
			    if($dss[$j]-$ass[$j]>$context)
			    {
				$bottom=$ass[$j]+$context-$start;
			    }
			    else
			    {
				$bottom=$dss[$j]-1-$start;
			    }
			    my($site,$confident)=search_cryptic_ass($ass[$j]-$start+$alter_base,$alter_gene_seq,$upper+$alter_base,$bottom+$alter_base,$strand);
			    if($site==0) ## without cryptic splice site
			    {
				##$state_for_splice_acceptor=1;
				my @site_of_ass;
				my $specific_start;
				my $specific_end;
				if($j==0)
				{
				    $specific_start=$start+1;
				}
				else
				{
				    $specific_start=$dss[$j-1];
				}
				$specific_end=$dss[$j];
				foreach my $variant (@zero_one_gene)
				{
				    my ($pos_first,$pos_second,$allele_alt)=(split /:/,$variant)[1,2,4];
				    if($specific_start<$pos_first && $pos_second<$specific_end )
				    {
					my $alter_base=relative_position_altered($pos_first,@zero_one_gene);
					my $real_position=$pos_first+$alter_base-$start;
					my $ref_pos=$pos_first-$start;
					my @new_ass=check_negativestrand_ass($alter_gene_seq,$wild_gene_seq,$real_position,$ref_pos);
					push @site_of_ass,@new_ass;
				    }
				}
				if(@site_of_ass)
				{
				    push @ass_zero_one,@site_of_ass;
				}
				else
				{
				    $state_for_splice_acceptor=1;
				}
			    }
			    else #with cryptic splice site
			    {
				push @ass_zero_one,$site;
				push @confident_score,$confident;
			    }
			}
			####for dss
			my $dss_overlaps = $tree->fetch($dss[$j]-5,$dss[$j]+3+1);
			my @dss_variant;
			if(scalar(@$dss_overlaps)>0)
			{
			    foreach my $region (@$dss_overlaps)
			    {
				my $r_chr=(split /:/,$region)[0];
				next if ($chrom ne $r_chr);  # ignore if overlap is on different chromosome
				push @dss_variant,$region;
			    }
			}
			##reconstruct splicing consense region to predict splicing-altering variants
			my $zero_one_state_of_dss=0;
			if(@dss_variant )
			{
			    my @zero_one;
			    foreach my $var ( @dss_variant )
			    {
				push @zero_one,$var if $var =~/0\|1/ ||$var =~/1\|1/;
				#push @one_zero,$var if $var =~/1\|0/ ||$var =~/1\|1/;
			    }
			    #my $zero_one_state="unbroken";
			    if(@zero_one )
			    {
				my $dss_changed_bases=relative_position_altered($dss[$j],@zero_one_gene);
				my $wild=substr($wild_gene_seq,$dss[$j]-5-$start-1,9);
				my $alter=substr($alter_gene_seq,$dss[$j]-5-$start+$dss_changed_bases-1,9);
				my $wild_seq=reverse(to_complement($wild));
				my $alter_seq=reverse(to_complement($alter));
				my $wild_score=score5($wild_seq);
				my $alter_score=score5($alter_seq);
				$wild_score=$wild_score+0.001;
				my $final_score=($wild_score-$alter_score)/$wild_score;##to be continued;
				#my $zero_one_state;
				if($wild_score >0 && $final_score >= 0.217)
				{
				    $zero_one_state_of_dss=1;
				}
				elsif($wild_score >0 && $final_score < 0.217)
				{
				    $zero_one_state_of_dss=0;
				}
				else
				{
				    if($wild_score>$alter_score)
				    {
					$zero_one_state_of_dss=1;##unknow if as 1, we don't add this constraint
				    }
				    else
				    {
					$zero_one_state_of_dss=0;
				    }
				}
			    }
			}
			if($zero_one_state_of_dss==0)
			{
			    my $alter_base=relative_position_altered($dss[$j],@zero_one_gene);
			    my $tmp=$dss[$j]+$alter_base-$start;
			    push @dss_zero_one,$tmp;
			}
			else #search for cryptic dss
			{
			    $splice_state=1;
			    my $alter_base=relative_position_altered($dss[$j],@zero_one_gene);
			    my $upper=0;
			    my $bottom=0;
			    if($dss[$j]-$ass[$j]>$context)
			    {
				$upper=$dss[$j]-$context-$start;
			    }
			    else
			    {
				$upper=$ass[$j]+1-$start;
			    }
			    if(defined($ass[$j+1]))
			    {
				if($ass[$j+1]-$dss[$j]>$context)
				{
				    $bottom=$dss[$j]+$context-$start;
				}
				else
				{
				    $bottom=$ass[$j+1]-1-$start;
				}
			    }
			    else
			    {
				$bottom=$dss[$j]+$context-$start;
			    }
			    my($site,$confident)=search_cryptic_dss($dss[$j]-$start+$alter_base,$alter_gene_seq,$upper+$alter_base,$bottom+$alter_base,$strand);
			    if($site==0)## without splice site
			    {
				my @site_of_dss;
				my $specific_start;
				my $specific_end;
				if(defined($ass[$j+1]))
				{
				    $specific_end=$ass[$j+1];
				}
				else
				{
				    $specific_end=$end;
				}
				$specific_start=$ass[$j];
				foreach my $variant (@zero_one_gene)
				{
				    my ($pos_first,$pos_second,$allele_alt)=(split /:/,$variant)[1,2,4];
				    if($specific_start<$pos_first && $pos_second<$specific_end)
				    {
					my $alter_base=relative_position_altered($pos_first,@zero_one_gene);
					my $real_position=$pos_first+$alter_base-$start;
					my $ref_pos=$pos_first-$start;
					my @new_dss=check_negativestrand_dss($alter_gene_seq,$wild_gene_seq,$real_position,$ref_pos);
					push @site_of_dss,@new_dss;
				    }
				}
				if(@site_of_dss)
				{
				    push @dss_zero_one,@site_of_dss;
				}
				else
				{
				    $state_for_splice_donor=1;
				}
			    }
			    else
			    {
				push @dss_zero_one,$site;
				push @confident_score,$confident;
			    }
			}
		    }
		}

		## extract new transcript
		my $logo_for_splice_gain='N';
		if($splicegain ==1)## prediction of splice site gain
		{
		    my @dss_gain;
		    my @ass_gain;
		    if($strand eq "+")
		    {
			foreach my $variant (@zero_one_gene)
			{
			    my ($pos_first,$pos_second,$allele_alt)=(split /:/,$variant)[1,2,4];
			    my $alter_base=relative_position_altered($pos_first,@zero_one_gene);
			    my $real_position=$pos_first+$alter_base-$start;
			    my $ref_pos=$pos_first-$start;
			    my @site_of_dss=check_dss($alter_gene_seq,$wild_gene_seq,$real_position,$ref_pos);
			    my @site_of_ass=check_ass($alter_gene_seq,$wild_gene_seq,$real_position,$ref_pos);
			    push @dss_gain,@site_of_dss;
			    push @ass_gain,@site_of_ass;
			}
		    }
		    else
		    {
			foreach my $variant (@zero_one_gene)
			{
			    my ($pos_first,$pos_second,$allele_alt)=(split /:/,$variant)[1,2,4];
			    my $alter_base=relative_position_altered($pos_first,@zero_one_gene);
			    my $real_position=$pos_first+$alter_base-$start;
			    my $ref_pos=$pos_first-$start;
			    my @site_of_dss=check_negativestrand_dss($alter_gene_seq,$wild_gene_seq,$real_position,$ref_pos);
			    my @site_of_ass=check_negativestrand_ass($alter_gene_seq,$wild_gene_seq,$real_position,$ref_pos);
			    push @dss_gain,@site_of_dss;
			    push @ass_gain,@site_of_ass;
			}
		    }
		    #my $logo_for_splice_gain='N';
		    #my @novel_dss; ## the real splice gain site 
		    #my @novel_ass;
		    if(@dss_gain)
		    {
			foreach my $gain (@dss_gain)
			{
			    unless($gain ~~ @dss_zero_one)
			    {
				$logo_for_splice_gain='G';
				last;
				#push @novel_dss,$gain;
			    }
			}
		    }
		    if(@ass_gain)
		    {
			foreach my $gain (@ass_gain)
			{
			    unless($gain ~~ @ass_zero_one)
			    {
				$logo_for_splice_gain='G';
				last;
				#push @novel_ass,$gain;
			    }
			}
		    }
		}
		#push @dss_zero_one,@novel_dss;
		#push @ass_zero_one,@novel_ass;
		my $logo;
		if($splice_state==0)
		{
		    $logo="O";
		}
		elsif($splice_state==1 && $state_for_splice_donor==0 && $state_for_splice_acceptor==0)
		{
		    $logo="R"; ## strands for rescued
		}
		else
		{
		    $logo="B"; ## strands for broken
		}
		if(@dss_zero_one && @ass_zero_one)
		{
		    my $dss_string=join("_",@dss_zero_one);
		    my $ass_string=join("_",@ass_zero_one);
		    #open(OUT,">start_dss_ass.txt")or die "cannot open the file";
		    #print OUT $relative_start,"\t",join(" ",@dss_zero_one),"\t",join(" ",@ass_zero_one),"\n";
		    #close OUT;
		    my @new_structure=split /\n/,`python ./script/extract_all_possible_combination_version3.py $strand $relative_start $dss_string $ass_string`;
		    #open(REGION,"out.txt") or die "can not open the file";
		    foreach my $line (@new_structure)
		    {
			chomp $line;
			next if $line=~/dot/; ## a warning reported by python networkx module
			my @site=split / /,$line;
			my $n=@site;
			my @seq=split //,$alter_gene_seq;
			my $length=@seq;
			##final splicing site in transcript
			my $final_splicing_site;
			my $gap=0;
			my $cds_rna_seq;
			if ($strand eq "+")
			{
			    for( my $i=0;$i<=$n-2;$i=$i+2 )
			    {
				foreach my $j ( $site[$i]-1..$site[$i+1]-1 )
				{
				    $seq[$j]="";
				}
			    }
			    for( my $i=0;$i<=$n-4;$i=$i+2 )
			    {
				$gap=$gap+$site[$i+1]-$site[$i]+1;
			    }
			    pop @site;
			    my $tmp=pop @site;
			    $final_splicing_site=$tmp-$gap;
			    foreach my $a ($relative_start-1 .. $length-1)
			    {
				unless($seq[$a] eq "")
				{
				    $cds_rna_seq =$cds_rna_seq . $seq[$a];
				}
			    }
			}
			else
			{
			    for( my $i=0;$i<=$n-2;$i=$i+2 )
			    {
				foreach my $j ( $site[$i+1]-1..$site[$i]-1 )
				{
				    $seq[$j]="";
				}
			    }
			    $final_splicing_site=pop @site;
			    my $temp;
			    foreach my $a (0 .. $relative_start-1)
			    {
				unless($seq[$a] eq "")
				{
				    $temp = $temp. $seq[$a];
				}
			    }
			    $cds_rna_seq=to_complement($temp);
			}
			my $rna=$cds_rna_seq;
			#print $rna,"\n";
			my($pro_seq,$end)=translate_rna($rna,$strand);
			## sequence alignment and parse
			my $ref_pro=$pro_database->seq($gene_id);
			$pm->finish(0) unless $ref_pro;
			my $ref_length=length($ref_pro);
			my $align_result;
			if($pro_seq)
			{
			    open(ALT,">$gene_id.0.1.protein_alt.fa") or die;
			    print ALT ">transcript\n",to_fasta($pro_seq),"\n";
			    close ALT;
			    open(REF,">$gene_id.0.1.protein_ref.fa") or die;
			    print REF ">transcript\n",to_fasta($ref_pro),"\n";
			    close REF;
			    `$needle_path/needle -asequence $gene_id.0.1.protein_ref.fa -sprotein1 -bsequence $gene_id.0.1.protein_alt.fa -sprotein2 -gapopen 10 -gapextend 0.5 -outfile $gene_id.0.1.pro_needel.txt 2> $gene_id.0.1.needle.log`;
			    $align_result=parse_needle("$gene_id.0.1.pro_needel.txt");
			    unlink "$gene_id.0.1.protein_alt.fa";
			    unlink "$gene_id.0.1.protein_ref.fa";
			    unlink "$gene_id.0.1.pro_needel.txt";
			    unlink "$gene_id.0.1.needle.log";
			}
			else
			{
			    $align_result="full deletion";
			}
			unless(@splice_score)
			{
			    @splice_score="--";
			}
			unless(@confident_score)
			{
			    @confident_score="--";
			}
			if($align_result)
			{
			    my $output=$gene_id;
			    if($bed==1)
			    {
				if($splicegain==1)
				{
				    $output=$output.":$chrom:$t_start:$t_end" ."\t$GeneName\t0|1\t$logo+$logo_for_splice_gain\t$ref_length\t". $align_result."\t" . join(",",@zero_one_gene);
				}
				else
				{
				    $output=$output.":$chrom:$t_start:$t_end" ."\t$GeneName\t0|1\t$logo\t$ref_length\t". $align_result."\t" . join(",",@zero_one_gene);
				}
				if($score==1)
				{
				    $output=$output."\t".join(",",@splice_score);
				}
				if($confident==1)
				{
				    $output=$output."\t".join(",",@confident_score);
				}
			    }
			    else ## default output
			    {
				if($splicegain==1)
				{
				    $output=$output."\t$GeneName\t0|1\t$logo+$logo_for_splice_gain\t$ref_length\t". $align_result ."\t" . join(",",@zero_one_gene);
				}
				else
				{
				    $output=$output."\t$GeneName\t0|1\t$logo\t$ref_length\t". $align_result."\t" . join(",",@zero_one_gene);
				}
				if($score==1)
				{
				    $output=$output."\t".join(",",@splice_score);
				}
				if($confident==1)
				{
				    $output=$output."\t".join(",",@confident_score);
				}
			    }
			    push @novel_seq,$output;
			}
			if($end == 0)
			{
			    ##nonstop decay
			}
			else
			{
			    my $end_position;
			    $end_position=$end+$relative_start-1 if $strand eq "+";
			    $end_position=$end if $strand eq "-";
			    if(abs($final_splicing_site-$end_position) > 50)
			    {
				#Nonsense mediated decay
				#print "$gene_id 0|1\t",abs($final_splicing_site-$end_position),"\n";
			    }
			    else
			    {
				# normal transcript
				#print "$gene_id 0|1\t",$pro_seq,"\n";
			    }
			}
		    }
		}
		else ##for single exon
		{
		    my @seq=split //,$alter_gene_seq;
		    my $length=@seq;
		    my $cds_rna_seq;
		    if($strand eq "+")
		    {
			#my $cds_rna_seq;
			foreach my $i ($relative_start-1 .. $length-1)
			{
			    $cds_rna_seq .= $seq[$i];
			}
		    }
		    else
		    {
			my $temp;
			foreach my $a (0 .. $relative_start-1)
			{
			    $temp .= $seq[$a];
			}
			$cds_rna_seq=to_complement($temp);
		    }
		    my $rna=$cds_rna_seq;
		    my($pro_seq,$end)=translate_rna($rna,$strand);
		    #print "$gene_id 0|1 S\t", $pro_seq,"\n";
		    ## alignment
		    my $ref_pro=$pro_database->seq($gene_id);
		    $pm->finish(0) unless $ref_pro;
		    my $ref_length=length($ref_pro);
		    my $align_result;
		    if($pro_seq)
		    {
			open(ALT,">$gene_id.0.1.protein_alt.fa") or die;
			print ALT ">transcript\n",to_fasta($pro_seq),"\n";
			close ALT;
			open(REF,">$gene_id.0.1.protein_ref.fa") or die;
			print REF ">transcript\n",to_fasta($ref_pro),"\n";
			close REF;
			`$needle_path/needle -asequence $gene_id.0.1.protein_ref.fa -sprotein1 -bsequence $gene_id.0.1.protein_alt.fa -sprotein2 -gapopen 10 -gapextend 0.5 -outfile $gene_id.0.1.pro_needel.txt 2> $gene_id.0.1.needle.log`;
			$align_result=parse_needle("$gene_id.0.1.pro_needel.txt");
			unlink "$gene_id.0.1.protein_alt.fa";
			unlink "$gene_id.0.1.protein_ref.fa";
			unlink "$gene_id.0.1.pro_needel.txt";
			unlink "$gene_id.0.1.needle.log";
		    }
		    else
		    {
			$align_result="full deletion";
		    }
		    if($align_result)
		    {
			if($bed==1)
			{
			    push @novel_seq,$gene_id.":$chrom:$t_start:$t_end" ."\t$GeneName\t0|1\tS\t$ref_length\t". $align_result ."\t" . join(",",@zero_one_gene);
			}
			else
			{
			    push @novel_seq,$gene_id."\t$GeneName\t0|1\tS\t$ref_length\t". $align_result ."\t" . join(",",@zero_one_gene);
			}
		    }
		}
		#unlink "start_dss_ass.txt";
		#unlink "out.txt";
	    }
	    if(@one_zero_gene)
	    {
		my $splice_state=0;
		my $state_for_splice_donor=0;
		my $state_for_splice_acceptor=0;## 0 strands for unbroken
		my($wild_gene_seq,@alter_seq_array)=reconstruct_seq($database,$chrom,$start+1,$end,$strand,@one_zero_gene);
		$pm->finish(0) if $wild_gene_seq=~/N/;
		my $alter_gene_seq=join("",@alter_seq_array);
		my @dss_one_zero;
		my @ass_one_zero;
		my $relative_start;
		if($strand eq "+")
		{
		    my $start_codon_altered=relative_position_altered($cds_start+1,@one_zero_gene);
		    $relative_start=$cds_start+1-$start+$start_codon_altered;
		    my $start_codon_first=$cds_start+1;
		    my @dss;
		    my @ass;
		    my $n=@length;
		    foreach my $i (0..$n-2)
		    {
			push @dss,$start+$start[$i]+$length[$i]+1; ##1-based
			push @ass,$start+$start[$i+1];##1-based
		    }
		    foreach my $j (0..$n-2)
		    {
			##for dss
			my $dss_overlaps = $tree->fetch($dss[$j]-3,$dss[$j]+5+1);
			my @dss_variant;
			if(scalar(@$dss_overlaps)>0)
			{
			    foreach my $region (@$dss_overlaps)
			    {
				my $r_chr=(split /:/,$region)[0];
				next if ($chrom ne $r_chr);  # ignore if overlap is on different chromosome
				push @dss_variant,$region;
			    }
			}
			##reconstruct splicing consense region to predict splicing-altering variants
			my $one_zero_state_of_dss=0;
			if(@dss_variant)
			{
			    my (@one_zero);
			    foreach my $var ( @dss_variant )
			    {
				#push @zero_one,$var if $var =~/0\|1/ ||$var =~/1\|1/;
				push @one_zero,$var if $var =~/1\|0/ ||$var =~/1\|1/;
			    }
			    if(@one_zero)
			    {
				my $dss_changed_bases=relative_position_altered($dss[$j],@one_zero_gene);
				my $wild_seq=substr($wild_gene_seq,$dss[$j]-3-$start-1,9);
				my $alter_seq=substr($alter_gene_seq,$dss[$j]-3+$dss_changed_bases-$start-1,9);
				my $wild_score=score5($wild_seq);
				my $alter_score=score5($alter_seq);
				$wild_score=$wild_score+0.001;
				my $final_score=($wild_score-$alter_score)/$wild_score;##to be continued;
				#my $zero_one_state;
				if($wild_score >0 && $final_score >= 0.217)
				{
				    $one_zero_state_of_dss=1;
				    push @splice_score,$final_score;
				}
				elsif($wild_score >0 && $final_score < 0.217)
				{
				    $one_zero_state_of_dss=0;
				}
				else
				{
				    if($wild_score>$alter_score)
				    {
					$one_zero_state_of_dss=1;##unknow if as 1, we don't add this constraint
				    }
				    else
				    {
					$one_zero_state_of_dss=0;
				    }
				}
			    }
			}
			if($one_zero_state_of_dss==0)
			{
			    my $alter_base=relative_position_altered($dss[$j],@one_zero_gene);
			    my $tmp=($dss[$j]+$alter_base-$start);
			    push @dss_one_zero,$tmp;
			}
			else #search for cryptic dss site
			{
			    $splice_state=1;
			    my $alter_base=relative_position_altered($dss[$j],@one_zero_gene);
			    my $upper=0;
			    my $bottom=0;
			    if($j==0)
			    {
				$upper=$dss[$j]-$start-$context if $dss[$j]-$start-$context>0;
				$upper=1 if $dss[$j]-$start-$context <= 0;
			    }
			    else
			    {
				if($dss[$j]-$ass[$j-1]>$context)
				{
				    $upper=$dss[$j]-$context-$start;
				}
				else
				{
				    $upper=$ass[$j-1]-$start+1;
				}
			    }
			    if($ass[$j]-$dss[$j]>$context)
			    {
				$bottom=$dss[$j]+$context-$start;
			    }
			    else
			    {
				$bottom=$ass[$j]-1-$start;
			    }
			    my($site,$confident)=search_cryptic_dss($dss[$j]-$start+$alter_base,$alter_gene_seq,$upper+$alter_base,$bottom+$alter_base,$strand);
			    if($site==0) ## no cryptic splice site
			    {
				my @site_of_dss;
				my $specific_start;
				my $specific_end;
				if($j==0)
				{
				    $specific_start=$start+1;
				}
				else
				{
				    $specific_start=$ass[$j-1];
				}
				$specific_end=$ass[$j];
				foreach my $variant (@one_zero_gene)
				{
				    my ($pos_first,$pos_second,$allele_alt)=(split /:/,$variant)[1,2,4];
				    if($specific_start<$pos_first && $pos_second<$specific_end) ## in specific region
				    {
					my $alter_base=relative_position_altered($pos_first,@one_zero_gene);
					my $real_position=$pos_first+$alter_base-$start;
					my $ref_pos=$pos_first-$start;
					my @new_dss=check_dss($alter_gene_seq,$wild_gene_seq,$real_position,$ref_pos);
					push @site_of_dss,@new_dss;
				    }
				}
				if(@site_of_dss)
				{
				    push @dss_one_zero,@site_of_dss;
				}
				else
				{
				    $state_for_splice_donor=1; ## ##means broken original structure
				}
			    }
			    else
			    {
				push @dss_one_zero,$site;
				push @confident_score,$confident;
			    }
			}
			
			####for ass
			my $ass_overlaps = $tree->fetch($ass[$j]-19,$ass[$j]+3+1);
			my @ass_variant;
			if(scalar(@$ass_overlaps)>0)
			{
			    foreach my $region (@$ass_overlaps)
			    {
				my $r_chr=(split /:/,$region)[0];
				next if ($chrom ne $r_chr);  # ignore if overlap is on different chromosome
				push @ass_variant,$region;
			    }
			}
			##reconstruct splicing consense region to predict splicing-altering variants
			my $one_zero_state_of_ass=0;
			if(@ass_variant )
			{
			    my @one_zero;
			    foreach my $var ( @ass_variant )
			    {
				push @one_zero,$var if $var =~/1\|0/ ||$var =~/1\|1/;
			    }
			    #my $zero_one_state="unbroken";
			    if(@one_zero )
			    {
				my $ass_changed_bases=relative_position_altered($ass[$j],@one_zero_gene);
				my $wild_seq=substr($wild_gene_seq,$ass[$j]-19-$start-1,23);
				my $alter_seq=substr($alter_gene_seq,$ass[$j]-19+$ass_changed_bases-$start-1,23);
				my $wild_score=score3($wild_seq);
				my $alter_score=score3($alter_seq);
				$wild_score=$wild_score+0.001;
				my $final_score=($wild_score-$alter_score)/$wild_score;##to be continued;
				#my $zero_one_state;
				if($wild_score >0 && $final_score >= 0.217)
				{
				    $one_zero_state_of_ass=1;
				    push @splice_score,$final_score;
				}
				elsif($wild_score >0 && $final_score < 0.217)
				{
				    $one_zero_state_of_ass=0;
				}
				else
				{
				    if($wild_score>$alter_score)
				    {
					$one_zero_state_of_ass=1;##unknow if as 1, we don't add this constraint
				    }
				    else
				    {
					$one_zero_state_of_ass=0;
				    }
				}
			    }
			}
			if($one_zero_state_of_ass ==0)
			{
			    my $alter_base=relative_position_altered($ass[$j],@one_zero_gene);
			    my $tmp=($ass[$j]+$alter_base-$start);
			    push @ass_one_zero,$tmp;
			}
			else
			{
			    $splice_state=1;
			    my $alter_base=relative_position_altered($ass[$j],@one_zero_gene);
			    my $upper=0;
			    my $bottom=0;
			    if($ass[$j]-$dss[$j]>$context)
			    {
				$upper=$ass[$j]-$context-$start;
			    }
			    else
			    {
				$upper=$dss[$j]+1-$start;
			    }
			    if(defined($dss[$j+1]))
			    {
				if($dss[$j+1]-$ass[$j]>$context)
				{
				    $bottom=$ass[$j]+$context-$start;
				}
				else
				{
				    $bottom=$dss[$j+1]-1-$start;
				}
			    }
			    else
			    {
				$bottom=$ass[$j]-$start+$context;
			    }
			    my($site,$confident)=search_cryptic_ass($ass[$j]-$start+$alter_base,$alter_gene_seq,$upper+$alter_base,$bottom+$alter_base,$strand);
			    if($site==0)## without splice site
			    {
				#$state_for_splice_acceptor=1;
				my @site_of_ass;
				my $specific_start;
				my $specific_end;
				if(defined($dss[$j+1]))
				{
				    $specific_end=$dss[$j+1];
				}
				else
				{
				    $specific_end=$end;
				}
				$specific_start=$dss[$j];
				foreach my $variant (@one_zero_gene)
				{
				    my ($pos_first,$pos_second,$allele_alt)=(split /:/,$variant)[1,2,4];
				    if($specific_start<$pos_first && $pos_second<$specific_end) ## in specific region
				    {
					my $alter_base=relative_position_altered($pos_first,@one_zero_gene);
					my $real_position=$pos_first+$alter_base-$start;
					my $ref_pos=$pos_first-$start;
					my @new_ass=check_ass($alter_gene_seq,$wild_gene_seq,$real_position,$ref_pos);
					push @site_of_ass,@new_ass;
				    }
				}
				if(@site_of_ass)
				{
				    push @ass_one_zero,@site_of_ass;
				}
				else
				{
				    $state_for_splice_acceptor=1; ## means broken original structure
				}
			    }
			    else
			    {
				push @ass_one_zero,$site;
				push @confident_score,$confident;
			    }
			}
		    }
		}
		else ##negative strand
		{
		    my $start_codon_altered=relative_position_altered($cds_end,@one_zero_gene);
		    $relative_start=$cds_end-$start+$start_codon_altered;
		    my @dss;
		    my @ass;
		    my $n=@length;
		    foreach my $i (0..$n-2)
		    {
			push @ass,$start+$start[$i]+$length[$i]+1; ##1-based ##need to change
			push @dss,$start+$start[$i+1];##1-based
		    }
		    foreach my $j (0..$n-2)
		    {
			####for ass
			my $ass_overlaps = $tree->fetch($ass[$j]-3,$ass[$j]+19+1);
			my @ass_variant;
			if(scalar(@$ass_overlaps)>0)
			{
			    foreach my $region (@$ass_overlaps)
			    {
				my $r_chr=(split /:/,$region)[0];
				next if ($chrom ne $r_chr);  # ignore if overlap is on different chromosome
				push @ass_variant,$region;
			    }
			}
			##reconstruct splicing consense region to predict splicing-altering variants
			my $one_zero_state_of_ass=0;
			if(@ass_variant)
			{
			    my @one_zero;
			    foreach my $var ( @ass_variant )
			    {
				#push @zero_one,$var if $var =~/0\|1/ ||$var =~/1\|1/;
				push @one_zero,$var if $var =~/1\|0/ ||$var =~/1\|1/;
			    }
			    #my $zero_one_state="unbroken";
			    if(@one_zero )
			    {
				my $ass_changed_bases=relative_position_altered($ass[$j],@one_zero_gene);
				my $wild=substr($wild_gene_seq,$ass[$j]-3-$start-1,23);
				my $alter=substr($alter_gene_seq,$ass[$j]-3-$start+$ass_changed_bases-1,23);
				my $wild_seq=reverse(to_complement($wild));
				my $alter_seq=reverse(to_complement($alter));
				my $wild_score=score3($wild_seq);
				my $alter_score=score3($alter_seq);
				$wild_score=$wild_score+0.001;
				my $final_score=($wild_score-$alter_score)/$wild_score;##to be continued;
				#my $zero_one_state;
				if($wild_score >0 && $final_score >= 0.217)
				{
				    $one_zero_state_of_ass=1;
				    push @splice_score,$final_score;
				}
				elsif($wild_score >0 && $final_score < 0.217)
				{
				    $one_zero_state_of_ass=0;
				}
				else
				{
				    if($wild_score>$alter_score)
				    {
					$one_zero_state_of_ass=1;##unknow if as 1, we don't add this constraint
				    }
				    else
				    {
					$one_zero_state_of_ass=0;
				    }
				}
			    }
			}
			if($one_zero_state_of_ass==0)
			{
			    my $alter_base=relative_position_altered($ass[$j],@one_zero_gene);
			    my $tmp=$ass[$j]+$alter_base-$start;
			    push @ass_one_zero,$tmp;
			}
			else #search for cryptic ass
			{
			    $splice_state=1;
			    my $alter_base=relative_position_altered($ass[$j],@one_zero_gene);
			    my $upper=0;
			    my $bottom=0;
			    if($j==0)
			    {
				if($ass[$j]-$start-$context<1)
				{
				    $upper=1;
				}
				else
				{
				    $upper=$ass[$j]-$start-$context;
				}
			    }
			    else
			    {
				if($ass[$j]-$dss[$j-1]>$context)
				{
				    $upper=$ass[$j]-$context-$start;
				}
				else
				{
				    $upper=$dss[$j-1]+1-$start;
				}
			    }
			    if($dss[$j]-$ass[$j]>$context)
			    {
				$bottom=$ass[$j]+$context-$start;
			    }
			    else
			    {
				$bottom=$dss[$j]-1-$start;
			    }
			    my($site,$confident)=search_cryptic_ass($ass[$j]-$start+$alter_base,$alter_gene_seq,$upper+$alter_base,$bottom+$alter_base,$strand);
			    if($site==0)## without cryptic splice site
			    {
				my @site_of_ass;
				my $specific_start;
				my $specific_end;
				if($j==0)
				{
				    $specific_start=$start+1;
				}
				else
				{
				    $specific_start=$dss[$j-1];
				}
				$specific_end=$dss[$j];
				foreach my $variant (@one_zero_gene)
				{
				    my ($pos_first,$pos_second,$allele_alt)=(split /:/,$variant)[1,2,4];
				    if($specific_start<$pos_first && $pos_second<$specific_end )
				    {
					my $alter_base=relative_position_altered($pos_first,@one_zero_gene);
					my $real_position=$pos_first+$alter_base-$start;
					my $ref_pos=$pos_first-$start;
					my @new_ass=check_negativestrand_ass($alter_gene_seq,$wild_gene_seq,$real_position,$ref_pos);
					push @site_of_ass,@new_ass;
				    }
				}
				if(@site_of_ass)
				{
				    push @ass_one_zero,@site_of_ass;
				}
				else
				{
				    $state_for_splice_acceptor=1;
				}
			    }
			    else
			    {
				push @ass_one_zero,$site;
				push @confident_score,$confident;
			    }
			}
			
			####for dss
			my $dss_overlaps = $tree->fetch($dss[$j]-5,$dss[$j]+3+1);
			my @dss_variant;
			if(scalar(@$dss_overlaps)>0)
			{
			    foreach my $region (@$dss_overlaps)
			    {
				my $r_chr=(split /:/,$region)[0];
				next if ($chrom ne $r_chr);  # ignore if overlap is on different chromosome
				push @dss_variant,$region;
			    }
			}
			##reconstruct splicing consense region to predict splicing-altering variants
			my $one_zero_state_of_dss=0;
			if(@dss_variant )
			{
			    my @one_zero;
			    foreach my $var ( @dss_variant )
			    {
				#push @zero_one,$var if $var =~/0\|1/ ||$var =~/1\|1/;
				push @one_zero,$var if $var =~/1\|0/ ||$var =~/1\|1/;
			    }
			    #my $zero_one_state="unbroken";
			    if(@one_zero )
			    {
				my $dss_changed_bases=relative_position_altered($dss[$j],@one_zero_gene);
				my $wild=substr($wild_gene_seq,$dss[$j]-5-$start-1,9);
				my $alter=substr($alter_gene_seq,$dss[$j]-5-$start+$dss_changed_bases-1,9);
				my $wild_seq=reverse(to_complement($wild));
				my $alter_seq=reverse(to_complement($alter));
				my $wild_score=score5($wild_seq);
				my $alter_score=score5($alter_seq);
				$wild_score=$wild_score+0.001;
				my $final_score=($wild_score-$alter_score)/$wild_score;##to be continued;
				#my $zero_one_state;
				if($wild_score >0 && $final_score >= 0.217)
				{
				    $one_zero_state_of_dss=1;
				    push @splice_score,$final_score;
				}
				elsif($wild_score >0 && $final_score < 0.217)
				{
				    $one_zero_state_of_dss=0;
				}
				else
				{
				    if($wild_score>$alter_score)
				    {
					$one_zero_state_of_dss=1;##unknow if as 1, we don't add this constraint
				    }
				    else
				    {
					$one_zero_state_of_dss=0;
				    }
				}
			    }
			}
			if($one_zero_state_of_dss==0)
			{
			    my $alter_base=relative_position_altered($dss[$j],@one_zero_gene);
			    my $tmp=$dss[$j]+$alter_base-$start;
			    push @dss_one_zero,$tmp;
			}
			else #search for cryptic dss
			{
			    $splice_state=1;
			    my $alter_base=relative_position_altered($dss[$j],@one_zero_gene);
			    my $upper=0;
			    my $bottom=0;
			    if($dss[$j]-$ass[$j]>$context)
			    {
				$upper=$dss[$j]-$context-$start;
			    }
			    else
			    {
				$upper=$ass[$j]+1-$start;
			    }
			    if(defined($ass[$j+1]))
			    {
				if($ass[$j+1]-$dss[$j]>$context)
				{
				    $bottom=$dss[$j]+$context-$start;
				}
				else
				{
				    $bottom=$ass[$j+1]-1-$start;
				}
			    }
			    else
			    {
				$bottom=$dss[$j]+$context-$start;
			    }
			    my($site,$confident)=search_cryptic_dss($dss[$j]-$start+$alter_base,$alter_gene_seq,$upper+$alter_base,$bottom+$alter_base,$strand);
			    if($site==0)
			    {
				my @site_of_dss;
				my $specific_start;
				my $specific_end;
				if(defined($ass[$j+1]))
				{
				    $specific_end=$ass[$j+1];
				}
				else
				{
				    $specific_end=$end;
				}
				$specific_start=$ass[$j];
				foreach my $variant (@one_zero_gene)
				{
				    my ($pos_first,$pos_second,$allele_alt)=(split /:/,$variant)[1,2,4];
				    if($specific_start<$pos_first && $pos_second<$specific_end)
				    {
					my $alter_base=relative_position_altered($pos_first,@one_zero_gene);
					my $real_position=$pos_first+$alter_base-$start;
					my $ref_pos=$pos_first-$start;
					my @new_dss=check_negativestrand_dss($alter_gene_seq,$wild_gene_seq,$real_position,$ref_pos);
					push @site_of_dss,@new_dss;
				    }
				}
				if(@site_of_dss)
				{
				    push @dss_one_zero,@site_of_dss;
				}
				else
				{
				    $state_for_splice_donor=1;
				}
			    }
			    else
			    {
				push @dss_one_zero,$site;
				push @confident_score,$confident;
			    }
			}
		    }
		}
		## extract new transcript
		my $logo_for_splice_gain='N';
		if($splicegain==1)## prediction of splice site gain
		{
		    my @dss_gain;
		    my @ass_gain;
		    if($strand eq "+")
		    {
			foreach my $variant (@one_zero_gene)
			{
			    my ($pos_first,$pos_second,$allele_alt)=(split /:/,$variant)[1,2,4];
			    my $alter_base=relative_position_altered($pos_first,@one_zero_gene);
			    my $real_position=$pos_first+$alter_base-$start;
			    my $ref_pos=$pos_first-$start;
			    my @site_of_dss=check_dss($alter_gene_seq,$wild_gene_seq,$real_position,$ref_pos);
			    my @site_of_ass=check_ass($alter_gene_seq,$wild_gene_seq,$real_position,$ref_pos);
			    push @dss_gain,@site_of_dss;
			    push @ass_gain,@site_of_ass;
			}
		    }
		    else
		    {
			foreach my $variant (@one_zero_gene)
			{
			    my ($pos_first,$pos_second,$allele_alt)=(split /:/,$variant)[1,2,4];
			    my $alter_base=relative_position_altered($pos_first,@one_zero_gene);
			    my $real_position=$pos_first+$alter_base-$start;
			    my $ref_pos=$pos_first-$start;
			    my @site_of_dss=check_negativestrand_dss($alter_gene_seq,$wild_gene_seq,$real_position,$ref_pos);
			    my @site_of_ass=check_negativestrand_ass($alter_gene_seq,$wild_gene_seq,$real_position,$ref_pos);
			    push @dss_gain,@site_of_dss;
			    push @ass_gain,@site_of_ass;
			}
		    }
		    #my $logo_for_splice_gain='N';
		    #my @novel_dss; ## the real splice gain site
		    #my @novel_ass;
		    if(@dss_gain)
		    {
			foreach my $gain (@dss_gain)
			{
			    unless($gain ~~ @dss_one_zero)
			    {
				$logo_for_splice_gain='G';
				last;
				#push @novel_dss,$gain;
			    }
			}
		    }
		    if(@ass_gain)
		    {
			foreach my $gain (@ass_gain)
			{
			    unless($gain ~~ @ass_one_zero)
			    {
				$logo_for_splice_gain='G';
				last;
				#push @novel_ass,$gain;
			    }
			}
		    }
		}
		#push @dss_one_zero,@novel_dss;
		#push @ass_one_zero,@novel_ass;
		my $logo;
		if($splice_state==0)
		{
		    $logo="O";
		}
		elsif($splice_state==1 && $state_for_splice_donor==0 && $state_for_splice_acceptor==0)
		{
		    $logo="R"; ## strands for rescued
		}
		else
		{
		    $logo="B"; ## strands for broken
		}
		if(@dss_one_zero && @ass_one_zero)
		{
		    my $dss_string=join("_",@dss_one_zero);
		    my $ass_string=join("_",@ass_one_zero);
		    #open(OUT,">start_dss_ass.txt")or die "cannot open the file";
		    #print OUT $relative_start,"\t",join(" ",@dss_one_zero),"\t",join(" ",@ass_one_zero),"\n";
		    #close OUT;
		    my @new_structure=split /\n/,`python ./script/extract_all_possible_combination_version3.py $strand $relative_start $dss_string $ass_string`;
		    #open(REGION,"out.txt") or die "can not open the file";
		    foreach my $line (@new_structure) #while(my $line=<REGION>)
		    {
			next if $line=~/dot/;
			chomp $line;
			my @site=split / /,$line;
			my $n=@site;
			my @seq=split //,$alter_gene_seq;
			my $length=@seq;
			##final splicing site in transcript
			my $final_splicing_site;
			my $gap=0;
			my $cds_rna_seq;
			if ($strand eq "+")
			{
			    for( my $i=0;$i<=$n-2;$i=$i+2 )
			    {
				foreach my $j ( $site[$i]-1..$site[$i+1]-1 )
				{
				    $seq[$j]="";
				}
			    }
			    for( my $i=0;$i<=$n-4;$i=$i+2 )
			    {
				$gap=$gap+$site[$i+1]-$site[$i]+1;
			    }
			    pop @site;
			    my $tmp=pop @site;
			    $final_splicing_site=$tmp-$gap;
			    foreach my $a ($relative_start-1 .. $length-1)
			    {
				unless($seq[$a] eq "")
				{
				    $cds_rna_seq .= $seq[$a];
				}
			    }
			}
			else
			{
			    for( my $i=0;$i<=$n-2;$i=$i+2 )
			    {
				foreach my $j ( $site[$i+1]-1..$site[$i]-1 )
				{
				    $seq[$j]="";
				}
			    }
			    $final_splicing_site=pop @site;
			    my $temp;
			    foreach my $a (0 .. $relative_start-1)
			    {
				unless($seq[$a] eq "")
				{
				    $temp .= $seq[$a];
				}
			    }
			    $cds_rna_seq=to_complement($temp);
			}
			my $rna=$cds_rna_seq;
			my($pro_seq,$end)=translate_rna($rna,$strand);
			## sequence alignment and parse
			my $ref_pro=$pro_database->seq($gene_id);
			$pm->finish(0) unless $ref_pro;
			my $ref_length=length($ref_pro);
			my $align_result;
			if($pro_seq)
			{
			    open(ALT,">$gene_id.1.0.protein_alt.fa") or die;
			    print ALT ">transcript\n",to_fasta($pro_seq),"\n";
			    close ALT;
			    open(REF,">$gene_id.1.0.protein_ref.fa") or die;
			    print REF ">transcript\n",to_fasta($ref_pro),"\n";
			    close REF;
			    `$needle_path/needle -asequence $gene_id.1.0.protein_ref.fa -sprotein1 -bsequence $gene_id.1.0.protein_alt.fa -sprotein2 -gapopen 10 -gapextend 0.5 -outfile $gene_id.1.0.pro_needel.txt 2> $gene_id.1.0.needle.log`;
			    $align_result=parse_needle("$gene_id.1.0.pro_needel.txt");
			    unlink "$gene_id.1.0.protein_alt.fa";
			    unlink "$gene_id.1.0.protein_ref.fa";
			    unlink "$gene_id.1.0.pro_needel.txt";
			    unlink "$gene_id.1.0.needle.log";
			}
			else
			{
			    $align_result="full deletion";
			}
			unless(@splice_score)
			{
			    @splice_score="--";
			}
			unless(@confident_score)
			{
			    @confident_score="--";
			}
			if($align_result)
			{
			    my $output=$gene_id;
			    if($bed==1)
			    {
				if($splicegain==1)
				{
				    $output=$output.":$chrom:$t_start:$t_end" ."\t$GeneName\t1|0\t$logo+$logo_for_splice_gain\t$ref_length\t". $align_result."\t" . join(",",@one_zero_gene);
				}
				else
				{
				    $output=$output.":$chrom:$t_start:$t_end" ."\t$GeneName\t1|0\t$logo\t$ref_length\t". $align_result."\t" . join(",",@one_zero_gene);
				}
				if($score==1)
				{
				    $output=$output."\t".join(",",@splice_score);
				}
				if($confident==1)
				{
				    $output=$output."\t".join(",",@confident_score);
				}
			    }
			    else ## default output
			    {
				if($splicegain==1)
				{
				    $output=$output."\t$GeneName\t1|0\t$logo+$logo_for_splice_gain\t$ref_length\t". $align_result ."\t" . join(",",@one_zero_gene);
				}
				else
				{
				    $output=$output."\t$GeneName\t1|0\t$logo\t$ref_length\t". $align_result."\t" . join(",",@one_zero_gene);
				}
				if($score==1)
				{
				    $output=$output."\t".join(",",@splice_score);
				}
				if($confident==1)
				{
				    $output=$output."\t".join(",",@confident_score);
				}
			    }
			    push @novel_seq,$output;
			}
			if($end == 0)
			{
			    ##nonstop decay
			}
			else
			{
			    my $end_position;
			    $end_position=$end+$relative_start-1 if $strand eq "+";
			    $end_position=$end if $strand eq "-";
			    if(abs($final_splicing_site-$end_position) > 50)
			    {
				#Nonsense mediated decay
				#print "$gene_id 1|0\t",abs($final_splicing_site-$end_position),"\n";
			    }
			    else
			    {
				# normal transcript
				#print "$gene_id 1|0\t",$pro_seq,"\n";
			    }
			}
		    }
		}
		else ##single_exon 
		{
		    my @seq=split //,$alter_gene_seq;
		    my $length=@seq;
		    my $cds_rna_seq;
		    if($strand eq "+")
		    {
			#my $cds_rna_seq;
			foreach my $i ($relative_start-1 .. $length-1)
			{
			    $cds_rna_seq .= $seq[$i];
			}
		    }
		    else
		    {
			my $temp;
			foreach my $a (0 .. $relative_start-1)
			{
			    $temp .= $seq[$a];
			}
			$cds_rna_seq=to_complement($temp);
		    }
		    my $rna=$cds_rna_seq;
		    my($pro_seq,$end)=translate_rna($rna,$strand);
		    #print "$gene_id 1|0 S\t", $pro_seq,"\n";
		    ## alignment
		    my $ref_pro=$pro_database->seq($gene_id);
		    $pm->finish(0) unless $ref_pro;
		    my $ref_length=length($ref_pro);
		    my $align_result;
		    if($pro_seq)
		    {
			open(ALT,">$gene_id.1.0.protein_alt.fa") or die;
			print ALT ">transcript\n",to_fasta($pro_seq),"\n";
			close ALT;
			open(REF,">$gene_id.1.0.protein_ref.fa") or die;
			print REF ">transcript\n",to_fasta($ref_pro),"\n";
			close REF;
			`$needle_path/needle -asequence $gene_id.1.0.protein_ref.fa -sprotein1 -bsequence $gene_id.1.0.protein_alt.fa -sprotein2 -gapopen 10 -gapextend 0.5 -outfile $gene_id.1.0.pro_needel.txt 2> $gene_id.1.0.needle.log`;
			$align_result=parse_needle("$gene_id.1.0.pro_needel.txt");
			unlink "$gene_id.1.0.protein_alt.fa";
			unlink "$gene_id.1.0.protein_ref.fa";
			unlink "$gene_id.1.0.pro_needel.txt";
			unlink "$gene_id.1.0.needle.log";
		    }
		    else
		    {
			$align_result="full deletion";
		    }
		    if($align_result)
		    {
			if($bed ==1 )
			{
			    push @novel_seq,$gene_id.":$chrom:$t_start:$t_end" ."\t$GeneName\t1|0\tS\t$ref_length\t". $align_result."\t" . join(",",@one_zero_gene);
			}
			else
			{
			    push @novel_seq,$gene_id."\t$GeneName\t1|0\tS\t$ref_length\t". $align_result."\t" . join(",",@one_zero_gene);
			}
		    }
		}
	    }
	    #unlink "start_dss_ass.txt";
	    #unlink "out.txt";
	}
    }
    $pm->finish(0,\@novel_seq);
}
$pm->wait_all_children;

if($bed == 0)
{
    print join("\n",@total_seq),"\n";
}
else ## output bed format
{
    foreach my $item (@total_seq)
    {
	my ($info,$genesymbol,$hap,$code,$length,$change,$variant,$score,$confident)=split /\t/,$item;
	my($gene,$chr,$start,$end)=split /:/,$info;
	print "$chr\t$start\t$end\tTranscript=$gene;Gene=$genesymbol;Haplotype=$hap;Splicing_Code=$code;Protein_Length=$length;Amino_Acid=$change","Variant=$variant;";
	if($score)
	{
	    print "Score=$score;";
	}
	if($confident)
	{
	    print "splicing_confident_score=$confident;"
	}
	print "\n";
    }
}

sub to_complement
{
    my $DNA=shift;

    # The Perl translate/transliterate command is just what we need:
    $DNA =~ tr/ACGTacgt/TGCAtgca/;
    return $DNA;
}

sub reconstruct_seq
{
    my($db,$chr,$first,$second,$strand,@variation)=@_;
    my $ref_seq=$db->seq($chr,$first,$second);
    my @var_seq=split //,$ref_seq;
    foreach my $mutation ( @variation )
    {
	my ($pos_first,$pos_second,$allele_alt)=(split /:/,$mutation)[1,2,4];
	$pos_first = $pos_first-$first+1; #to change the start to relative position
	$pos_second = $pos_second-$first+1; # to change the end position
	#$info[3]="" if $info[3] eq "-";
	#$info[4]="" if $info[4] eq "-";
	for my $i ($pos_first..$pos_second)
	{
	        $var_seq[$i-1]="";
	    }
	$var_seq[$pos_first-1]=$allele_alt;
    }
    return($ref_seq,@var_seq);
}

sub relative_position_altered
{
    my($pos,@variant)=@_;
    my $alter_base=0;
    foreach my $variation (@variant)
    {
	my($chr,$start,$end,$ref,$alt)=(split /:/,$variation)[0,1,2,3,4];
	if($end<=$pos)
	{
	        $alter_base=$alter_base+length($alt)-length($ref);
	    }
    }
    return $alter_base;
}

sub overlap
{
    my($start1,$end1,$start2,$end2)=@_;
    my $state=0;
    unless($start1 > $end2 ||$end1 < $start2)
    {
	$state=1;
    }
    return $state; ## 1 strands for overlapped
}

sub to_fasta
{
    my $seq=shift;
    my $n=length($seq);
    my @seq_array=(split //,$seq);
    my @fasta="";
    for my $i(1..$n)
    {
	push @fasta,$seq_array[$i-1];
	push @fasta, "\n" if $i%60==0 && $i<$n;
    }
     return @fasta;
}

sub check_ass
{
    my($db,$refdb,$position,$ref_pos)=@_;
    ##sequence
    my $seq=substr($db,$position-22-1,45);
    my $ref_seq=substr($refdb,$ref_pos-22-1,45);
    my @sequence=split //,$seq;
    my $nn=@sequence-1;
    my $mm=length($ref_seq)-1;
    my $n;
    if($nn<=$mm)
    {
	$n=$nn;
    }
    else
    {
	$n=$mm;
    }
    my @ass;
    foreach my $i (0 .. $n-22)
    {
	my $part=substr($seq,$i,23);
	my $refpart=substr($ref_seq,$i,23);
	if($part=~/[ATGC]{18}AG[ATGC]{3}/ && (!($refpart=~/[ATGC]{18}AG[ATGC]{3}/)))
	{
	    #my $refpart=substr($ref_seq,$i,23);
	    my $state=score3fun($part,$refpart);
	    my $site=$i+19+$position-22;
	    if($state==1)
	    {
		my $branch=substr($db,$site-100-1,89);
		my $state_branch=check_branch($branch);
		if($state_branch==1)
		{
		    push @ass,$site;
		}
	    }
	}
    }
    return @ass;
}


sub check_negativestrand_ass
{
    my($db,$refdb,$position,$ref_pos)=@_;
    my $seq_tmp=substr($db,$position-22-1,45);
    my $ref_seq_tmp=substr($refdb,$ref_pos-22-1,45);
    my $seq=to_complement($seq_tmp);
    my $ref_seq=to_complement($ref_seq_tmp);
    my @sequence=split //,$seq;
    my $nn=@sequence-1;
    my @ass;
    my $n;
    my $mm=length($ref_seq)-1;
    if($nn<=$mm)
    {
	$n=$nn;
    }
    else
    {
	$n=$mm;
    }
    foreach my $i (0 .. $n-22)
    {
        my $part_tmp=substr($seq,$i,23);
	my $part=reverse $part_tmp;
	my $refpart_tmp=substr($ref_seq,$i,23);
	my $refpart=reverse $refpart_tmp;
	if($part=~/[ATGC]{18}AG[ATGC]{3}/ && (!($refpart=~/[ATGC]{18}AG[ATGC]{3}/)))
	{
	    #my $refpart_tmp=substr($ref_seq,$i,23);
	    #my $refpart=reverse $refpart_tmp;
	    my $state=score3fun($part,$refpart);
	    my $site=$i+19+$position-22;
	    if($state==1)
	    {
		my $seq_br=substr($db,$site+12-1,89);
		my $branch=to_complement($seq_br);
		my $state_branch=check_branch($branch);
		if($state_branch ==1)
		{
		    push @ass,$site;
		}
	    }
	}
    }
    return @ass;
}

sub check_branch
{
    my $branch=shift;
    my $le=length($branch);
    my $state=0;
    foreach my $i (0 .. $le-5)
    {
	my $part=substr($branch,$i,5);
	if($part=~/[ATGC]{3}A[ATGC]/)
	{
	    my @sequence=split //,$part;
	    my $score=0;
	    my $n=@sequence;
	    foreach my $i ( 0 .. $n-1 )
	    {
		$score=$score+$motif{'branch'}[$i]{$sequence[$i]};
	    }
	    $state=1 if $score >0;
	    last if $state ==1;
	}
    }
    return $state;
}

sub score3fun
{
    my($sequence,$ref)=@_;
    my $wild_score=score3($ref);
    my $alt_score=score3($sequence);
    $wild_score=$wild_score+0.001;
	#print $ref,"\t",$sequence,"\t",$wild_score,"\t",$alt_score,"\tchengsj\n";
    my $relative=($alt_score-$wild_score)/$wild_score;
    my $state=0;# means no ass
    if($wild_score >0 && $alt_score > 3 && $relative > 0.217)
    {
	$state=1;
    }
    elsif($alt_score > 3 && $wild_score<0)
    {
	$state=1;
    }
    return $state;
}

sub check_dss
{
    my($db,$refdb,$position,$ref_pos)=@_;
    ##sequence
    my $seq=substr($db,$position-8-1,17);
    my $ref_seq=substr($refdb,$ref_pos-8-1,17);
    my $nn=length($seq)-1;
    my $mm=length($ref_seq)-1;
    my $n;
    if($nn<=$mm)
    {
	$n=$nn;
    }
    else
    {
	$n=$mm;
    }
    my @dss;
    foreach my $i (0 .. $n-8)
    {
        my $part=substr($seq,$i,9);
	my $refpart=substr($ref_seq,$i,9);
	if($part=~/[ATGC]{3}GT[ATGC]{4}/ && (!($refpart=~/[ATGC]{3}GT[ATGC]{4}/)))
	{
	    #my $refpart=substr($ref_seq,$i,9);
	    my $state=score5fun($part,$refpart);
	    push @dss,$i+3+$position-8 if $state==1;
	}
    }
    return @dss;
}

sub check_negativestrand_dss
{
    my($db,$refdb,$position,$ref_pos)=@_;
    ##sequence
    my $seq_tmp=substr($db,$position-8-1,17);
    my $seq=to_complement($seq_tmp);
    my $ref_seq_tmp=substr($refdb,$ref_pos-8-1,17);
    my $ref_seq=to_complement($ref_seq_tmp);
    my $nn=length($seq)-1;
    my $mm=length($ref_seq)-1;
    my $n;
    if($nn<=$mm)
    {
	$n=$nn;
    }
    else
    {
	$n=$mm;
    }
    my @dss;
    foreach my $i (0 .. $n-8)
    {
        my $part=reverse(substr($seq,$i,9));
	my $refpart=reverse(substr($ref_seq,$i,9));
	if($part=~/[ATGC]{3}GT[ATGC]{4}/ && (!($refpart=~/[ATGC]{3}GT[ATGC]{4}/)))
	{
	    #my $refpart=reverse(substr($ref_seq,$i,9));
	    my $state=score5fun($part,$refpart);
	    push @dss,$i+3+$position-8 if $state==1;
	}
    }
    return @dss;
}

sub score5fun
{
    my($sequence,$ref)=@_;
    my $wild_score=score5($ref);
    my $alt_score=score5($sequence);
    $wild_score=$wild_score+0.001;
	#print $wild_score,"\t",$alt_score;
	#open(TT,">chengsj");
	#print TT $ref,"\n",$sequence unless defined($alt_score);
    my $relative=($alt_score-$wild_score)/$wild_score;
    my $state=0;# means no ass
    if($wild_score <3 && $alt_score > 3 && $relative > 0.217)
    {
        $state=1;
    }
    elsif($wild_score<0 && $alt_score>3)
    {
	$state=1;
    }
    return $state;
}

sub translate_rna
{
    my %aacode = (
		  TTT => "F", TTC => "F", TTA => "L", TTG => "L",
		  TCT => "S", TCC => "S", TCA => "S", TCG => "S",
		  TAT => "Y", TAC => "Y", TAA => "", TAG => "",
		  TGT => "C", TGC => "C", TGA => "", TGG => "W",
		  CTT => "L", CTC => "L", CTA => "L", CTG => "L",
		  CCT => "P", CCC => "P", CCA => "P", CCG => "P",
		  CAT => "H", CAC => "H", CAA => "Q", CAG => "Q",
		  CGT => "R", CGC => "R", CGA => "R", CGG => "R",
		  ATT => "I", ATC => "I", ATA => "I", ATG => "M",
		  ACT => "T", ACC => "T", ACA => "T", ACG => "T",
		  AAT => "N", AAC => "N", AAA => "K", AAG => "K",
		  AGT => "S", AGC => "S", AGA => "R", AGG => "R",
		  GTT => "V", GTC => "V", GTA => "V", GTG => "V",
		  GCT => "A", GCC => "A", GCA => "A", GCG => "A",
		  GAT => "D", GAC => "D", GAA => "E", GAG => "E",
		  GGT => "G", GGC => "G", GGA => "G", GGG => "G",
		  #TGG =>"W",
		 ); # this is the hash table for the amino acids
    my ($seq,$strand)=@_;
    my $length=length $seq;
    my @sequence=split //,$seq;
    my $aminoacid;
    my $end_position=0;
    if($strand eq "+")
    {
	for(my $i=0; $i< $length;$i=$i+3)
	{
	    if(defined $sequence[$i+1] && defined $sequence[$i+2])
	    {
		my $codon=$sequence[$i].$sequence[$i+1].$sequence[$i+2];
		#print $codon,"\n";
		#print $aacode{$codon},"\n";
		$aminoacid .= $aacode{$codon};
		$end_position += 3;
		if($codon eq 'TAA' || $codon eq 'TGA' ||$codon eq 'TAG')
		{
		    last;
		}
	    }
	    else
	    {
		$end_position=0;
	    }
	}
    }
    else
    {
	my $tmp_end=0;
	for(my $i=$length-1;$i>=0;$i=$i-3)
	{
	    if(defined $sequence[$i-1] && defined $sequence[$i-2])
	    {
		my $codon=$sequence[$i].$sequence[$i-1].$sequence[$i-2];
		$aminoacid .= $aacode{$codon};
		$tmp_end += 3;
		if($codon eq 'TAA' || $codon eq 'TGA' ||$codon eq 'TAG')
                {
                    last;
                }
	    }
	    else
            {
                $tmp_end=0;
            }
	}
	if($tmp_end ==0)
	{
	    $end_position=0;
	}
	else
	{
	    $end_position=$length-$tmp_end+1;
	}
    }
    return $aminoacid,$end_position;
}

sub search_cryptic_dss
{
    my($pos,$db,$upper,$bottom,$strand)=@_;
    my $l=$bottom-$upper+1;
    my $seq=substr($db,$upper-1,$l);
    my $length=length($seq);
    my $cryptic_site=0; ##means no cryptic splicing site
    my $confident=0; ##  means no cryptic splicing site
    my $begin=$pos-$upper;
    if($strand eq "+")
    {
	foreach my $i (1..$length)
	{
	    if($begin-$i-3>0 && $begin-$i-3+8< $length)
	    {
		my $part=substr($seq,$begin-$i-3,9);
		if($part=~/[ATGC]{3}GT[ATGC]{4}/)
		{
		    my $score=score5($part);
		    if($score > 3)
		    {
			$cryptic_site=$begin-$i+$upper;
			$confident=$score;
			last;
		    }
		}
	    }
	    if($begin+$i+5<$length)
	    {
		my $part=substr($seq,$begin+$i-3,9);
		if($part=~/[ATGC]{3}GT[ATGC]{4}/)
                {
                    my $score=score5($part);
                    if($score > 3)
                    {
                        $cryptic_site=$begin+$i+$upper;
			$confident=$score;
                        last;
                    }
                }
	    }
	}
    }
    else #negative strand
    {
	foreach my $i (1..$length)
	{
            if($begin-$i-5>0)
            {
                my $tmp=substr($seq,$begin-$i-5,9);
		my $part=reverse(to_complement($tmp));
                if($part=~/[ATGC]{3}GT[ATGC]{4}/)
                {
                    my $score=score5($part);
                    if($score > 3)
                    {
                        $cryptic_site=$begin-$i+$upper;
			$confident=$score;
                        last;
                    }
                }
            }
            if($begin+$i+3<$length)
            {
                my $tmp=substr($seq,$begin+$i-5,9);
		my $part=reverse(to_complement($tmp));
                if($part=~/[ATGC]{3}GT[ATGC]{4}/)
                {
                    my $score=score5($part);
                    if($score > 3)
                    {
                        $cryptic_site=$begin+$i+$upper;
			$confident=$score;
                        last;
                    }
                }
            }
        }
    }
    return $cryptic_site,$confident;
}

sub search_cryptic_ass
{
    my($pos,$db,$upper,$bottom,$strand)=@_;
    my $l=$bottom-$upper+1;
    my $seq=substr($db,$upper-1,$l);
    my $length=length($seq);
    my $cryptic_site=0; ##means no cryptic splicing site
    my $confident=0;
    my $begin=$pos-$upper;
    if($strand eq "+")
    {
	foreach my $i (1..$length)
	{
	    if($begin-$i-19>0)
	    {
		my $part=substr($seq,$begin-$i-19,23);
		if($part=~/[ATGC]{18}AG[ATGC]{3}/)
		{
		    my $score=score3($part);
		    if($score > 3)
		    {
			$cryptic_site=$begin-$i+$upper;
			$confident=$score;
			last;
		    }
		}
	    }
	    if($begin+$i+3<$length)
	    {
		my $part=substr($seq,$begin+$i-19,23);
		if($part=~/[ATGC]{18}AG[ATGC]{3}/)
		{
		    my $score=score3($part);
		    if($score > 3)
		    {
			$cryptic_site=$begin+$i+$upper;
			$confident=$score;
			last;
		    }
		}
	    }
	}
    }
    else #negative strand
    {
	foreach my $i (1..$length)
	{
            if($begin-$i-3>0)
            {
                my $tmp=substr($seq,$begin-$i-3,23);
		my $part=reverse(to_complement($tmp));
                if($part=~/[ATGC]{18}AG[ATGC]{3}/)
                {
                    my $score=score3($part);
                    if($score > 3)
                    {
                        $cryptic_site=$begin-$i+$upper;
			$confident=$score;
                        last;
                    }
                }
            }
            if($begin+$i+19<$length)
            {
                my $tmp=substr($seq,$begin+$i-3,23);
		my $part=reverse(to_complement($tmp));
                if($part=~/[ATGC]{18}AG[ATGC]{3}/)
                {
                    my $score=score3($part);
                    if($score > 3)
                    {
                        $cryptic_site=$begin+$i+$upper;
			$confident=$score;
                        last;
                    }
                }
            }
        }
    }
    return $cryptic_site;
}

## maxentscan 3ass scoring
sub asshashseq
{
    #returns hash of sequence in base 4
    # &hashseq('CAGAAGT') returns 4619
    my $seq = shift;
    $seq = uc($seq);
    $seq =~ tr/ACGT/0123/;
    my @seqa = split(//,$seq);
    my $sum = 0;
    my $len = length($seq);
    my @four = (1,4,16,64,256,1024,4096,16384);
    my $i=0;
    while ($i<$len)
    {
	$sum+= $seqa[$i] * $four[$len - $i -1] ;
	$i++;
    }
    return $sum;
}
sub assmaxentscore{
    my $seq = shift;
    my $table_ref = shift;
    my @metables = @$table_ref;
    my @sc;
    $sc[0] = $metables[0]{&asshashseq(substr($seq,0,7))};
    $sc[1] = $metables[1]{&asshashseq(substr($seq,7,7))};
    $sc[2] = $metables[2]{&asshashseq(substr($seq,14,7))};
    $sc[3] = $metables[3]{&asshashseq(substr($seq,4,7))};
    $sc[4] = $metables[4]{&asshashseq(substr($seq,11,7))};
    $sc[5] = $metables[5]{&asshashseq(substr($seq,4,3))};
    $sc[6] = $metables[6]{&asshashseq(substr($seq,7,4))};
    $sc[7] = $metables[7]{&asshashseq(substr($seq,11,3))};
    $sc[8] = $metables[8]{&asshashseq(substr($seq,14,4))};
    my $finalscore = $sc[0] * $sc[1] * $sc[2] * $sc[3] * $sc[4] / ($sc[5] * $sc[6] * $sc[7] * $sc[8]);
    return $finalscore;
}
sub assgetrest{
    my $seq = shift;
    my $seq_noconsensus = substr($seq,0,18).substr($seq,20,3);
    return $seq_noconsensus;
}
sub assscoreconsensus{
    my $seq = shift;
    my @seqa = split(//,uc($seq));
    my %bgd;
    $bgd{'A'} = 0.27;
    $bgd{'C'} = 0.23;
    $bgd{'G'} = 0.23;
    $bgd{'T'} = 0.27;
    my %cons1;
    $cons1{'A'} = 0.9903;
    $cons1{'C'} = 0.0032;
    $cons1{'G'} = 0.0034;
    $cons1{'T'} = 0.0030;
    my %cons2;
    $cons2{'A'} = 0.0027;
    $cons2{'C'} = 0.0037;
    $cons2{'G'} = 0.9905;
    $cons2{'T'} = 0.0030;
    my $addscore = $cons1{$seqa[18]} * $cons2{$seqa[19]}/ ($bgd{$seqa[18]} * $bgd{$seqa[19]});
    return $addscore;
}
sub score3
{
    my $str=shift;
    my $score=sprintf("%.2f",&log2(&assscoreconsensus($str)*&assmaxentscore(&assgetrest($str),\@metables)));
    return $score;
}

## maxentscan 5dss scoring
sub dssgetrest{
    my $seq = shift;
    my @seqa = split(//,uc($seq));
    return $seqa[0].$seqa[1].$seqa[2].$seqa[5].$seqa[6].$seqa[7].$seqa[8];
}
sub dssscoreconsensus{
    my $seq = shift;
    my @seqa = split(//,uc($seq));
    my %bgd;
    $bgd{'A'} = 0.27;
    $bgd{'C'} = 0.23;
    $bgd{'G'} = 0.23;
    $bgd{'T'} = 0.27;
    my %cons1;
    $cons1{'A'} = 0.004;
    $cons1{'C'} = 0.0032;
    $cons1{'G'} = 0.9896;
    $cons1{'T'} = 0.0032;
    my %cons2;
    $cons2{'A'} = 0.0034;
    $cons2{'C'} = 0.0039;
    $cons2{'G'} = 0.0042;
    $cons2{'T'} = 0.9884;
    my $addscore = $cons1{$seqa[3]}*$cons2{$seqa[4]}/($bgd{$seqa[3]}*$bgd{$seqa[4]});
    return $addscore;
}
sub score5
{
    my $str=shift;
    my $score=sprintf("%.2f",&log2(&dssscoreconsensus($str)*$me2x5{$seq{&dssgetrest($str)}}));
    return $score;
}

sub log2{
    my ($val) = @_;
    return log($val)/log(2);
}

sub parse_needle
{
    my $file=shift;
    my $ref="";
    my $query="";
    my $symbol="";
    my $n=0;
    open(NEEDLE,$file) or print "died $file\n";
    while( my $line=<NEEDLE>)
    {
	chomp $line;
	if($line =~ /^#/ || $line eq "")
	{
	    next;
	}
	$n++;
	if($line =~ /[0-9]+\s([\S-]+)\s+[0-9]+/)
	{
	    #$n++;
	    $ref .= $1 if $n%3==1;
	    #$symbol=$1 if $n%3==2;
	    $query .= $1 if $n%3==0;
	}
	elsif($line =~ /^\s{21}([\s\S]+)$/)                       ##if($line =~ /^\s{21}([\s]*[\S]+[\s]*)$/ || $line =~/^\s{21}(\s+)$/ )
	{
	    $symbol .= $1;
	}
	#print "ref:$ref\nsymbol:$symbol\nquery:$query\n";
    }
    #print "$ref\n$symbol\n$query\n";
    my @sym=split //,$symbol;
    my @refseq=split //,$ref;
    my @queryseq=split //,$query;
    my $ref_gap=0;
    my $le=@refseq;
    my $start=0;
    my $end=0;
    my $a=0;
    my $result;
    foreach my $i ( 1..$le )
    {
	$i--;
	if($sym[$i] ne "|")
	{
	    $a++;
	    $start=$i if $a==1;
	    $end=$i if $a==1;
	    next if $a==1;
	    if($i-$end == 1)
	    {
		#print $i+1 .":$refseq[$i]->$queryseq[$i];";
		$end=$i;
	    }
	    elsif($i-$end > 1)
	    {
		my $ref_tmp="";
		my $alt_tmp="";
		my $gap=0;
		foreach my $tmp ($start..$end)
		{
		    $gap++ if $refseq[$tmp] eq "-";
		    $refseq[$tmp]="" if $refseq[$tmp] eq "-";
		    $queryseq[$tmp]="" if $queryseq[$tmp] eq "-";
		    $ref_tmp .=$refseq[$tmp];
		    $alt_tmp .=$queryseq[$tmp];
		}
		$start++;
		$end++;
		$start=$start-$ref_gap;
		$end=$end-$ref_gap-$gap;
		if($ref_tmp eq "")
		{
		    #$ref_gap=$ref_gap+$end-$start+1;
		    $start--;
		    $ref_tmp=$refseq[$start-1];
		    $alt_tmp=$queryseq[$start-1].$alt_tmp;
		    #$end=$start;
		}
		$result .= $start ."-" .$end.":$ref_tmp->$alt_tmp;";
		$start=$i;
		$end=$i;
		$ref_gap=$ref_gap+$gap;
		# $a=1;
	    }
	}
    }
    my $gap=0;
    if($start || $end)
    {
	my $ref_tmp="";
	my $alt_tmp="";
	foreach my $temp ($start..$end)
	{
	    $gap++ if $refseq[$temp] eq "-";
	    $refseq[$temp]="" if $refseq[$temp] eq "-";
	    $queryseq[$temp]="" if $queryseq[$temp] eq "-";
	    $ref_tmp .=$refseq[$temp];
	    $alt_tmp .=$queryseq[$temp];
	}
	$start++;
	$end++;
	$start=$start-$ref_gap;
	$end=$end-$ref_gap-$gap;
	if($ref_tmp eq "")
	{
	    #$ref_gap=$ref_gap+$end-$start+1;
	    $start--;
	    $ref_tmp=$refseq[$start-1];
	    $alt_tmp=$queryseq[$start-1].$alt_tmp;
	    #$end=$start;
	}
	$result .= $start ."-" .$end.":$ref_tmp->$alt_tmp;";
    }
    close NEEDLE;
    return $result;
}

