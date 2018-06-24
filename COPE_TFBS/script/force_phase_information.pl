#!/usr/bin/perl
use warnings;

## usage:
my $input=$ARGV[0];
my $output=$ARGV[1];

open(FH,$ARGV[0]) or die "cannot open the file!";
open(OUT,">$output") or die;
while(my $line=<FH>)
{
    chomp $line;
    if($line =~/^#/)
    {
	print OUT $line,"\n";
    }
    else
    {
	my @info=split /\t/,$line;
	print OUT "$info[0]\t$info[1]\t$info[2]\t$info[3]\t$info[4]\t$info[5]\t$info[6]\t$info[7]\t$info[8]\t0|1\n";
    }
}

close FH;
close OUT;
