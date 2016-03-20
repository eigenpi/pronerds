Author: Cristinel Ababei
Research collaborator: Rajesh Kavasseri
December 2009, Fargo ND
cristinel.ababei@ndsu.edu


Synopsis
========
This is the PRO-NERDS tool: efficient network reconfiguration of 
distribution systems for loss minimization. It is based on minimum-
cost maximum-flow algorithm. This tool is one result of my 
collaboration with Rajesh Kavasseri, also with NDSU.

The tool contains two network reconfiguration algorithms:
1.  Baran's method (used as reference)
2.  Proposed minimum-cost maximum-flow (mcmf) method: see our 
    paper for details [1]

Also, the tool contains two loss estimation techniques:
1.  DistFlow of Baran (used as reference)
2.  The proposed random-walks (rw) technique: see our paper 
    for details [1]

[1] C. Ababei and R. Kavasseri, Efficient network reconfiguration 
    using minimum cost maximum flow based branch exchanges and random 
    walks based loss estimations, IEEE Trans. on Power Systems, 
    vol. 26, no. 1, pp. 30-37, Feb. 2011.


Installation
============
The latest version of the tool can be downloaded from:
http://venus.ece.ndsu.nodak.edu/~cris/software.html
The tool was developed on a Linux machine running Fedora 8. 
It was last compiled on a Fedora 12. It should be compile-able 
on Windows too. On linux, first edit the Makefile to reflect 
the location where you want to compile and link. Then, just type:
> make


How to use the tool
===================
Type "nerds" at the command prompt to see the help menu.

Examples:
> nerds tests/bus_83_11.pos -reconfig_algo baran
> nerds tests/bus_83_11.pos -reconfig_algo mcmf 
> nerds tests/bus_83_11.pos -reconfig_algo mcmf -use_gui


Benchmarks
==========
The format is new. I proposed it in order to simplify the parsing
routines. See an example of benchmark file understood by the tool, 
at the end of this file as well as readme file from tests/.
All benchmarks (or testcases) are located tests/.


Things you might want to do
===========================
1. If you dig into it and find any bug, please let me know
2. If you use this code in a research project and want to 
include a reference to it, then please use:
[*] C. Ababei and R. Kavasseri, Efficient network reconfiguration 
    using minimum cost maximum flow based branch exchanges and random 
    walks based loss estimations, IEEE Trans. on Power Systems, 
    vol. 26, no. 1, pp. 30-37, Feb. 2011.
3. If you'll ever hit it big (to be read: make a lot of money :-) ),
and this code helped you in any way, then please consider 
donating some to support our research (we need it :-) )


Credits
=======
-- The initial version of CS2, which is here partially ported
to C++, was developed by Andrew Goldberg (goldberg@intertrust.com) and
Boris Cherkassky (cher@cher.msk.su)
-- Vaughn Betz (while at Univ. of Toronto) developed much of the GUI;


Copyright
=========
Copyright 2009-20xx by Cristinel Ababei (cristinel.ababei@ndsu.edu)
and Rajesh Kavasseri (rajesh.kavasseri@ndsu.edu).
This Copyright notice applies to all files, called hereafter 
"The Software".
Permission to use, copy, and modify this software and its 
documentation is hereby granted only under the following 
terms and conditions.  Both the above copyright notice and 
this permission notice must appear in all copies of the 
software, derivative works or modified versions, and any 
portions thereof, and both notices must appear in supporting 
documentation.  Permission is granted only for non-commercial 
use.  For commercial use, please contact the author.
This software may be distributed (but not offered for sale 
or transferred for compensation) to third parties, provided 
such third parties agree to abide by the terms and conditions
of this notice.
The Software is provided "as is", and the authors, the 
North Dakota State University (NDSU), as well as any and
all previous authors (of portions or modified portions of
the software) disclaim all warranties with regard to this 
software, including all implied warranties of merchantability
and fitness.  In no event shall the authors or NDSU or any and
all previous authors be liable for any special, direct, 
indirect, or consequential damages or any damages whatsoever
resulting from loss of use, data or profits, whether in an
action of contract, negligence or other tortious action,
arising out of or in connection with the use or performance
of this software.


Example of benchmark: bus_13_3.pos (Civanlar)
=============================================
.PU 0                                           <-- Testcase is given in per unit values or not. 0 No 1 Yes
.V_base                                         <-- in [ kV ]
.S_base                                         <-- in [ MVA ]
.Branch_count 16                                <-- Number of branches counted from 0..15
.Sectionalizing_count 13                        <-- How many branches are sectionalizing
.Tie_count 3                                    <-- How many are tie switches
.Nodes 16                                       <-- Nodes counted from 0..15
.Feeders 3                                      <-- Number of feeders in system
.Feeder_node_ids 0 1 2                          <-- Feeder IDs are from 0..Num_feeders - 1
.Branch Src_bus Rec_bus R X PL_kw QL_kvar S_NS  <-- Next lines characterize branches
0 0 3 0.075 0.1 2.0 1.6 S
1 3 4 0.08 0.11 3.0 1.5 S
2 3 5 0.09 0.18 2.0 0.8 S
3 5 6 0.04 0.04 1.5 1.2 S 
4 1 7 0.11 0.11 4.0 2.7 S
5 7 8 0.08 0.11 5.0 3.0 S
6 7 9 0.11 0.11 1.0 0.9 S 
7 8 10 0.11 0.11 0.6 0.1 S 
8 8 11 0.08 0.11 4.5 2.0 S
9 2 12 0.11 0.11 1.0 0.9 S
10 12 13 0.09 0.12 1.0 0.7 S 
11 12 14 0.08 0.11 1.0 0.9 S
12 14 15 0.04 0.04 2.1 1.0 S
.Tie switches          <-- Remaining branches are the tie switches
13 4 10 0.04 0.04 
14 9 13 0.04 0.04 
15 6 15 0.12 0.12 
.Node x_coord y_coord  <-- Next lines give the placement of each node. Used for GUI purposes only
.Grid_size 8 8         <-- GUI will display nodes on a grid of 8x8
0 1 5                  <-- Node "0" is placed at location 1,5
1 4 5
2 7 5
3 1 4
4 2 3
5 1 1 
6 3 1
7 4 4 
8 4 3
9 5 3
10 3 3
11 4 2
12 7 4
13 6 3
14 7 1
15 5 1
