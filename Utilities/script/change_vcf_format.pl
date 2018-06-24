#!/usr/bin/perl
use warnings;

open(VAR,$ARGV[0]) or die;
my $sampleid=$ARGV[1];
print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$sampleid\n";

while(<VAR>)
{
    chomp;
    next if $_=~/^#/;
    my @info=split /\t/,$_;
    #my @tmp=split / /,$info[7];
    my $type=(split /:/,$info[9])[0];
    print "$info[0]\t$info[1]\t$info[2]\t$info[3]\t$info[4]\t$info[5]\t$info[6]\t.\tGT\t$type\n";
}
