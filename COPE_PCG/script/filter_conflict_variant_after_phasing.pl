#!/usr/bin/perl

# filter_conflict_variant_after_phasing.pl 
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
