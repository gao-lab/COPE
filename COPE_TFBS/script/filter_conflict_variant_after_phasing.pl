#!/usr/bin/perl
use warnings;

## input sorted VCF
my @block;
my $file=$ARGV[0];
open(FILTER,">$file.remove.Conflicts.vcf") or die;
print FILTER "##fileformat=VCFv4.1\n";
open(FH,$file) or die;
while(my $line=<FH>)
{
    chomp $line;
    next if $line=~/^#/;
    unless(@block)
    {
	push @block,$line;
	next;
    }
    if(@block)
    {
	if(overlap_with_block(\@block,$line)) ## overlap with previous block
	{
	    push @block,$line;
	}
	else ## no overlap, print specific lines in blocks
	{
	    my $n=@block;
	    if($n==1)
	    {
		my $tmp=shift @block;
		my $result=change_to_union($tmp);
		my @array=@$result;
		print FILTER join("\n",@array),"\n";
	    }
	    else
	    {
		foreach my $tmp (@block)
		{
		    unless(overlap_consider_phase(\@block,$tmp)) ## note the item itself
		    {
			my $result=change_to_union($tmp);
			my @array=@$result;
			print FILTER join("\n",@array),"\n";
			#print FILTER $tmp,"\n";
		    }
		}
	    }
	    @block=();
	    push @block,$line;
	}
    }
}

close FH;

foreach my $tmp (@block)
{
    unless(overlap_consider_phase(\@block,$tmp)) ## note the item itself
    {
	print FILTER $tmp,"\n";
    }
}
close FILTER;

sub overlap_with_block
{
    my($line1,$line2)=@_;
    my($chr_2,$pos_2,$ref_2)=(split /\t/,$line2)[0,1,3];
    my $end_2=$pos_2+length($ref_2)-1;
    my $result=0;
    foreach my $pre (@$line1)
    {
	my($chr,$pos,$ref)=(split /\t/,$pre)[0,1,3];
	my $end=$pos+length($ref)-1;
	if($chr eq $chr_2)
	{
	    if(overlap($pos,$end,$pos_2,$end_2))
	    {
		$result=1;
		last;
	    }
	}
    }
    return $result;
}

sub overlap_consider_phase
{
    my($line1,$line2)=@_;
    my($pos_2,$ref_2,$hap_2)=(split /\t/,$line2)[1,3,9];
    my $end_2=$pos_2+length($ref_2)-1;
    my $result=0;
    foreach my $pre (@$line1)
    {
	next if $pre eq $line2;
	my($pos,$ref,$hap)=(split /\t/,$pre)[1,3,9];
	my $end=$pos+length($ref)-1;
	if(same_hap($hap,$hap_2))
	{
	    if(overlap($pos,$end,$pos_2,$end_2))
	    {
		$result=1;
		last;
	    }
	}
    }
    return $result;
}

sub overlap ## two region overlapped or not
{
    my($pos,$end,$pos_2,$end_2)=@_;
    my $result=1;
    if($pos>$end_2 || $end<$pos_2)
    {
	$result=0;
    }
    return $result;
}

sub same_hap
{
    my($hap1,$hap2)=@_;
    my $result=0;
    if($hap1=~/\|[123456789]/ && $hap2=~/\|[123456789]/)
    {
	$result=1;
    }
    elsif($hap1=~/[123456789]\|/ && $hap2=~/[123456789]\|/)
    {
	$result=1;
    }
    return $result;
}

sub change_to_union ## change haplotype to 0|1 and 1|0
{
    my $line=shift;
    my($chr,$pos,$id,$ref,$alt,$score,$filter,$info,$gt,$hap)=split /\t/,$line;
    my @alts=split /,/,$alt;
    my @result;
    if($hap =~ /0\|([123456789])/)
    {
	push @result,"$chr\t$pos\t$id\t$ref\t$alts[$1-1]\t$score\t$filter\t$info\t$gt\t0|1";
    }
    elsif($hap=~ /([123456789])\|0/)
    {
	push @result,"$chr\t$pos\t$id\t$ref\t$alts[$1-1]\t$score\t$filter\t$info\t$gt\t1|0";
    }
    elsif($hap =~/1\|1/)
    {
	push @result,"$chr\t$pos\t$id\t$ref\t$alts[0]\t$score\t$filter\t$info\t$gt\t1|1";
    }
    elsif($hap =~/([123456789])\|([123456789])/)
    {
	push @result,"$chr\t$pos\t$id\t$ref\t$alts[$1-1]\t$score\t$filter\t$info\t$gt\t1|0";
	push @result,"$chr\t$pos\t$id\t$ref\t$alts[$2-1]\t$score\t$filter\t$info\t$gt\t0|1";
    }
    return \@result;
}
