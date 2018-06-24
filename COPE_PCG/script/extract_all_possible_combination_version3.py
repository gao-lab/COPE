#!/usr/bin/python

# extract_all_possible_combination_version3.py 
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




import networkx as nx
import sys

strand=sys.argv[1] ##strand information
start=sys.argv[2] ##start
dss=sys.argv[3].split("_")
ass=sys.argv[4].split("_")

#start="10928" #sys.argv[2]
#dss=["381","11095","12289","13586","15742","17052","36990"]#sys.argv[3].split(" ")
#ass=["10921","12153","15578","16938","36774","38916"]#sys.argv[4].split(" ")

G=nx.DiGraph()
G.add_node(start)
G.add_nodes_from(dss)
G.add_nodes_from(ass)
nodes=[]
if strand == "+":
    nodes=sorted(dss + ass,key=int)
else:
    nodes=sorted(dss + ass,key=int)
    nodes.reverse()
n=0
state=0; # 1 strands for dss; 2 strands for ass
pre_node=[]
late_node=[]
start_dss=[]
end_ass=[]
end_of_state=0; # 1 strands for in dss, 2 strands for ass
state_for_start_dss=0
for item in nodes:
    if item in dss:
        if state == 0 or state == 1:
            pre_node.append(item)
            state=1
        if state == 2:
            if state_for_start_dss == 0:
                start_dss=pre_node
            for pre in pre_node:
                for pre_two in late_node:
                    G.add_edge(pre,pre_two)
            pre_node=[]
            pre_node.append(item)
            state=1
            end_of_state=1
            state_for_start_dss=1
    if item in ass:
        if state == 0:
            next
        if state==2:
            late_node.append(item)
            state=2
        if state==1:
            for item_a in late_node:
                for item_b in pre_node:
                    G.add_edge(item_a,item_b)
            late_node=[]
            late_node.append(item)
            state=2
            end_of_state=2
if end_of_state == 2:
    for i in pre_node:
        for j in late_node:
            G.add_edge(i,j)
end_ass=late_node
if start_dss:
    pass
else:
    start_dss=pre_node
##find all combination
#out=open("out.txt","w")
for start in start_dss:
    for end in end_ass:
        paths = nx.all_simple_paths(G, source=start, target=end)
        #print(list(paths))
        for ii in list(paths):
            print  ' '.join(ii)
#out.close()

