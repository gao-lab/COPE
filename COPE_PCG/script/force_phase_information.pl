#!/usr/bin/perl

# force_phase_information.pl 
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
