cristinel.ababei@ndsu.edu
Feb. 2009, Fargo ND

New benchmark format
====================
The format is new and simple. Introduced for the sake of easy 
parsing. The main differences from old formats are:
- Branches are counted from 0 to Branch_count - 1
- Buses (i.e., nodes) are also counted from 0 to Bus_count - 1
It also includes placement info to be used for visual display by 
tools equipped with a GUI. Normally not required.


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
