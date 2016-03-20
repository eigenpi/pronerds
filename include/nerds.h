#ifndef _NERDS_H_
#define _NERDS_H_

#include <assert.h>
#include <stdio.h>
#include <vector>
#include <deque>
#include <utility>
#include <functional>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <algorithm>
#include <climits>
#include "nerds_mcmf.h"

using namespace std;

#define RECONFIG_LOSS_REDUCTION_ERROR 1E-12
// limit the number of iters for runtime purposes;
#define RECONFIG_ITERATION_BOUND 512 // 512
// when partial RW based PFS is used we perform a full DistFlow
// after every other RW_PFS_BRANCH_EXCHANGES_BOUND;
#define RW_PFS_BRANCH_EXCHANGES_BOUND 100 // 100
// POWER_FLOW_SOLUTION_ERROR represents the accuracy with which
// we want the power flow solution to be calculated; this is just
// a percentage, and we want the delta error between two consecutive
// iteration to be less than that many percentages, in order to stop;
#define POWER_FLOW_SOLUTION_ERROR 0.1
#define POWER_FLOW_ITERATION_BOUND 10

enum LODGING_TYPE { HOME_NODE, MOTEL_NODE };
enum NERDS_NODE_TYPE { REGULAR, FEEDER, SUPER_FEEDER };
enum NERDS_ARC_TYPE { NOT_SECTIONALIZING, SECTIONALIZING, TIE };
enum RECONFIG_ALGORITHM { RECONFIG_MCMF = 0, RECONFIG_BARAN = 1 };
enum POWER_FLOW_METHOD { BARAN_DIST_FLOW = 0, RANDOM_WALK_POWER_FLOW_FULL = 1,
     RANDOM_WALK_POWER_FLOW_PARTIAL = 2 };

class GUI_GRAPHICS;

////////////////////////////////////////////////////////////////////////////////
//
// PAIR_TWO
//
////////////////////////////////////////////////////////////////////////////////

class PAIR_TWO
{
 private:
    int _from_i, _to_i;
 public:
    PAIR_TWO() { _from_i = 0; _to_i = 0; }
    PAIR_TWO(int from_i, int to_i) : _from_i(from_i), _to_i(to_i) {}
    PAIR_TWO(const PAIR_TWO &pt) : 
        _from_i(pt._from_i), _to_i(pt._to_i) {}
    ~PAIR_TWO() {}

    int from_i() const { return _from_i; }
    int to_i() const { return _to_i; }
    void set_from_i(int from_i) { _from_i = from_i; }
    void set_to_i(int to_i) { _to_i = to_i; }
};

////////////////////////////////////////////////////////////////////////////////
//
// TUPLE_THREE
//
////////////////////////////////////////////////////////////////////////////////

class TUPLE_THREE
{
 private:
    int _tie_i;
    int _sec_i; // sectionalizing branch id;
    double _loss; // if tie is closed and sec is open;
 public:
    TUPLE_THREE() { _tie_i = 0; _sec_i = 0; _loss = 0.0; }
    TUPLE_THREE(int tie_i, int sec_i, double loss) : 
        _tie_i(tie_i), _sec_i(sec_i), _loss(loss){}
    TUPLE_THREE(const TUPLE_THREE &tt) : 
        _tie_i(tt._tie_i), _sec_i(tt._sec_i), _loss(tt._loss) {}
    ~TUPLE_THREE() {}

    int tie_i() const { return _tie_i; }
    int sec_i() const { return _sec_i; }
    double loss() const { return _loss; }
    void set_tie_i(int tie_i) { _tie_i = tie_i; }
    void set_sec_i(int sec_i) { _sec_i = sec_i; }
    void set_loss(double loss) { _loss = loss; }
    void operator=(const TUPLE_THREE &tt) {
        _tie_i = tt._tie_i; _sec_i =tt._sec_i; _loss = tt._loss;
    }
    bool operator==(const TUPLE_THREE &tt) const { return (_loss == tt._loss); }
    bool operator!=(const TUPLE_THREE &tt) const { return (_loss != tt._loss); }
    bool operator<(const TUPLE_THREE &tt)  const { return (_loss < tt._loss); }
};

////////////////////////////////////////////////////////////////////////////////
//
// RECT
//
////////////////////////////////////////////////////////////////////////////////

class RECT 
{
 private:
    int _xl, _xr, _yb, _yt;
    int _area;
    double _dist; // from this rect to any other rect in the other set;
 public:
    RECT() { _xl = 0; _xr = 0; _yb = 0; _yt = 0; _area = 0; _dist = 0.; }
    RECT(int xl, int xr, int yb, int yt) : 
        _xl(xl), _xr(xr), _yb(yb), _yt(yt), 
        _area((xr-xl)*(yt-yb)), _dist(0.) {}
    RECT(const RECT &r) : 
        _xl(r._xl), _xr(r._xr), _yb(r._yb), _yt(r._yt),
        _area(r._area), _dist(r._dist) {}
    ~RECT() {}

    int xl() const { return _xl; }
    int xr() const { return _xr; }
    int yb() const { return _yb; }
    int yt() const { return _yt; }
    int area() const { return _area; }
    double dist() const { return _dist; }
    void set_xl(int xl) { _xl = xl; }
    void set_xr(int xr) { _xr = xr; }
    void set_yb(int yb) { _yb = yb; }
    void set_yt(int yt) { _yt = yt; }
    void set_area(int area) { _area = area; }
    void set_dist(double dist) { _dist = dist; }

    void operator=(const RECT &r) {
        _xl = r._xl; _xr = r._xr; _yb = r._yb; _yt = r._yt; 
        _area = r._area;
        _dist = r._dist;
    }
    bool operator==(const RECT &r) const { return (_area == r._area); }
    bool operator!=(const RECT &r) const { return (_area != r._area); }
    bool operator<(const RECT &r)  const { return (_dist < r._dist); }
};

////////////////////////////////////////////////////////////////////////////////
//
// PRONERDS
//
////////////////////////////////////////////////////////////////////////////////

// this is the place where everything is driven from;

class NERDS_NODE;
class NERDS_GRAPH;
class RECT;

class PRONERDS {
 private:
    // _nerds_graph represents the power distribution network graph;
    // nodes are the feeder-nodes and bus-nodes; arcs are the branches which
    // can have switches or not - and switches can be sectionalizing
    // or tie switches;
    NERDS_GRAPH *_nerds_graph;
    bool _use_gui;
    int _rng_seed; // seed for random number generator; else is set to 1;
    // Edmunds seems to be faster than CS2 of Goldberg;
    MCMF_ALGORITHM _mcmf_algo;
    // reconfig can be traditional Baran or new MCMF based;
    RECONFIG_ALGORITHM _reconfig_algo;
    // power flow solution method can be DistFlow or RandomWalk based;
    POWER_FLOW_METHOD _pfs_algo;
    // runtime monitors; used for runtime profiling of various components;
    // cpu time spent on power flow solutions and updates of it during 
    // running reconfiguration;
    double _runtime_power_flow; // cumulated power flow cpu only;
    double _runtime_all; // whole reconfig algo;

 public:
    GUI_GRAPHICS *_gui;
    
 public:
    PRONERDS();
    ~PRONERDS() {}
    
    NERDS_GRAPH *nerds_graph() { return _nerds_graph; };
    GUI_GRAPHICS *gui() { return _gui; };
    void set_gui(GUI_GRAPHICS *gui) { _gui = gui; };
    bool use_gui() { return _use_gui; }
    POWER_FLOW_METHOD pfs_algo() const { return _pfs_algo; }
    RECONFIG_ALGORITHM reconfig_algo() const { return _reconfig_algo; }
    MCMF_ALGORITHM mcmf_algo() const { return _mcmf_algo; }
    bool parse_command_arguments( int argc, char *argv[]);
    void print_initial_stats( int argc, char *argv[]);
    bool create_network_graph( int argc, char *argv[]);

    // algos;
    bool run_DistFlow_solution_baran( bool simplified_distflow = false,
        bool use_super_feeder_node = false);
    bool run_PowerFlow_solution_rw( long version);
    bool run_reconfiguration_method1_baran( bool simplified_distflow,
        POWER_FLOW_METHOD power_flow_solution = BARAN_DIST_FLOW); // default;
    bool run_reconfiguration_method_new( bool simplified_distflow,
        POWER_FLOW_METHOD power_flow_solution = BARAN_DIST_FLOW); // default;

    // debug;
    void print_sketch_arrays();
    // utils;
    void add_to_runtime_power_flow( double delta) { _runtime_power_flow += delta; }
    void set_runtime_power_flow( double val) { _runtime_power_flow = val; }
    double runtime_power_flow() const { return _runtime_power_flow; }
    void add_to_runtime_all( double delta) { _runtime_all += delta; }
    void set_runtime_all( double val) { _runtime_all = val; }
    double runtime_all() const { return _runtime_all; }
    static bool compare_by_area(const RECT &r1, const RECT &r2) {
        return ( r1.area() > r2.area());
    }
    static bool compare_by_dist(const RECT &r1, const RECT &r2) {
        // want smallest distance, largest area;
        if ( r1.dist() > r2.dist()) {
            return true;
        }
        else if ( (r1.dist() == r2.dist())&&(r1.area() < r2.area()) ) {
            return true;
        }       
        return false;
    }
    static int my_rand_int(int max_val) 
    {
        // return a random integer between [0..max_val];
        // r is a random float value in the range [0,1);
        double r = (double)rand() / ((double)(RAND_MAX)+(double)(1));
        r = r * (max_val + 1);
        return int( r);
    }
};

////////////////////////////////////////////////////////////////////////////////
//
// NERDS_ARC
//
////////////////////////////////////////////////////////////////////////////////

// an arc is a branch;

class NERDS_ARC {
 private:
    int _id;
    NERDS_ARC_TYPE _type; // NOT_SECTIONALIZING, SECTIONALIZING, TIE;
    long _src_id; // id of node that is source of this arc;
    long _des_id; // id of node that is destination of this arc;
    double _r;
    double _x;
    double _g; // 1/_r;
    // _P and _Q for arc i-1 -> i (using terminology from Baran paper), 
    // are the Pi and Qi for thru the arc going downstream from upstream
    // to downstream node; Note: every node has _P,_Q but with different meaning
    // they are accumulated values;
    double _P, _Q;
    bool _processed; // flag to keep track of arcs that are foward/backward processed;
        
 public:
    NERDS_ARC() { 
        _id = -1;
        _type = NOT_SECTIONALIZING;
        _processed = false;
    }
    NERDS_ARC(int id) : _id(id) {
        _type = NOT_SECTIONALIZING;
        _processed = false;
    }
    NERDS_ARC(const NERDS_ARC &arc) : _id(arc._id) {}
    ~NERDS_ARC() {}

    void clear() {
        _P = 0.0; _Q = 0.0;
    }   
    int id() const { return _id; }
    void set_id( int id) { _id = id; }
    NERDS_ARC_TYPE type() { return _type; }
    void set_type(NERDS_ARC_TYPE value) { _type = value; }
    int src_id() const { return _src_id; }
    int des_id() const { return _des_id; }
    double r() const { return _r; }
    double x() const { return _x; }
    double g() const { return _g; }
    void set_g( double g) { _g = g; }
    bool processed() const { return _processed; }
    double P() { return _P; }
    double Q() { return _Q; }
    void get_P_Q( double &p, double &q) { p = _P; q = _Q; }
    void set_P( double p) { _P = p; }
    void set_Q( double q) { _Q = q; }
    void set_P_Q( double p, double q) { _P = p; _Q = q; }
    void set_processed( bool val) { _processed = val; }
    void set( long src_id, long des_id, NERDS_ARC_TYPE type,
        double r, double x) {
        _src_id = src_id;
        _des_id = des_id;
        _type = type;
        _r = r;
        _x = x;
        // set _g; I use an extra class variable for this in order to save
        // computation time during random walks; I may hurt the mem usage 
        // though;
        // for dc case, this should store only _g = 1/_r; but for the ac case
        // it becomes _g = sqrt( 1/(arc->r()*arc->r() + arc->x()*arc->x()) );
        if ( _r > 0) {
            //_g = 1/_r; // use only g and neglect x;
            _g = sqrt( 1 / (_r * _r + _x * _x) );
        } else {
            _g = 10E6;
        }
    }
};

////////////////////////////////////////////////////////////////////////////////
//
// NERDS_NODE
//
////////////////////////////////////////////////////////////////////////////////

// node of network distribution graph; represents basically a feeder or bus;

class NERDS_NODE {
 private:
    // id of itself; 0..._nodes_count-1;
    int _id;
    // id of arc that is parent of this one viewed from the side of the
    // root node of the feeder_tree of which this node is part of; used in
    // backward forward sweeps;
    int _parent_arc_id;
    int _parent_id; // id parent node;
    int _feeder_id; // id of feeder (tree) to which this node belongs currently;
    int _depth; // from the root in the tree;
    NERDS_NODE_TYPE _type; // SUPER_FEEDER, FEEDER, REGULAR (i.e., a bus);
    double _pl_kw, _ql_kvar;
    double _v_k_1; // voltage at iteration k-1, ie previous iteration;
    double _v_k; // voltage at iteration k, ie current iteration;
    // _fans is the list of neighbors/adjacent nodes; among them its
    // parent (as we come from the root feeder) as well; 
    // id's of arcs entering or exiting this node;
    vector<long> _fans;
    // _prob_intervals "models" the transition probabilities; see nerds.cpp for
    // details;
    vector<PAIR_TWO> _prob_intervals;
    double _lodging_cost; // reward for home, penalty for motel;
    LODGING_TYPE _lodging_type;
    int _x_coord, _y_coord; // placement for nice gui purposes only;
    // number of nodes as children; this is different from _fans.size
    // because we do not include the parent neighbor and neither the nodes
    // neighbors via tie switches;
    long _children_count;
    long _processed_children; // used for sweeps;
    // _P and _Q for a node i, are the Pi and Qi for thru the arc going downstream to
    // the only child (in the case there is only one child, one fanout) or the
    // sum of all Pi and Qi thru all the downstream arcs/branches of this node;
    // these Pi and Qi are the ones appearing in eq 1, 2 in Baran paper;
    double _P, _Q;
    double _Pk_1, _Qk_1; // previous iter values;

 public:
    NERDS_NODE() { 
        _id = -1; 
        _parent_arc_id = -1;
        _parent_id = -1;
        _feeder_id = -1;
        _depth = -1;
        _type = REGULAR;
        _pl_kw = 0; _ql_kvar = 0; 
        _v_k = 1; _v_k_1 = 1;
        _children_count = 0;
        _processed_children = 0;
        _P = 1.; _Q = 1.; _Pk_1 = 1.; _Qk_1 = 1.;
        _x_coord = 0; _y_coord = 0;
        _lodging_cost = 0.0; _lodging_type = MOTEL_NODE;
    }
    NERDS_NODE( int id) : _id(id) {
        _parent_arc_id = -1;
        _parent_id = -1;
        _feeder_id = -1;
        _depth = -1;
        _type = REGULAR;
        _pl_kw = 0; _ql_kvar = 0;
        _v_k = 1; _v_k_1 = 1;
        _children_count = 0;
        _processed_children = 0;
        _P = 1.; _Q = 1.; _Pk_1 = 1.; _Qk_1 = 1.;
        _x_coord = 0; _y_coord = 0;
        _lodging_cost = 0.0; _lodging_type = MOTEL_NODE;
    }
    ~NERDS_NODE() {}

    void clear() {
        _v_k = 1; _v_k_1 = 1;
        _P = 1.; _Q = 1.; _Pk_1 = 1.; _Qk_1 = 1.;
        _lodging_cost = 0.0; 
        _lodging_type = MOTEL_NODE;
    }
    void clear_partially() {
        _lodging_type = MOTEL_NODE;
    }
    int id() const { return _id; }
    void set_id( int id) { _id = id; }
    void set_parent_arc_id( int id) { _parent_arc_id = id; }
    int parent_arc_id() const { return _parent_arc_id; }
    void set_parent_id( int id) { _parent_id = id; }
    int parent_id() const { return _parent_id; }
    void set_feeder_id( int id) { _feeder_id = id; }
    int feeder_id() const { return _feeder_id; }
    void set_depth( int id) { _depth = id; }
    int depth() const { return _depth; }
    NERDS_NODE_TYPE type() { return _type; }
    double pl_kw() { return _pl_kw; }
    double ql_kvar() { return _ql_kvar; }
    void set_pl_ql( double pl_kw, double ql_kvar) { 
        _pl_kw = pl_kw;
        _ql_kvar = ql_kvar;
    }
    double v_k_1() { return _v_k_1; }
    double v_k() { return _v_k; }
    void set_v_k_1( double v) { _v_k_1 = v; }
    void set_v_k( double v) { _v_k = v; }
    void update_v_k( double v) { _v_k = max( _v_k, v); }
    double P() { return _P; }
    double Q() { return _Q; }
    double Pk_1() { return _Pk_1; }
    double Qk_1() { return _Qk_1; }
    void copy_P_Q_to_Pk_1_Qk_1() { _Pk_1 = _P; _Qk_1 = _Q; }
    void set_P( double p) { _P = p; }
    void set_Q( double q) { _Q = q; }
    void set_P_Q( double p, double q) { _P = p; _Q = q; }
    void accumulate_P( double p) { _P += p; }
    void accumulate_Q( double q) { _Q += q; }
    int x_coord() const { return _x_coord; }
    int y_coord() const { return _y_coord; }
    void set_x_coord( int x_coord) { _x_coord = x_coord; }
    void set_y_coord( int y_coord) { _y_coord = y_coord; }
    void set_coords( int x_coord, int y_coord) { 
        _x_coord = x_coord; _y_coord = y_coord;
    }
    void set_type( NERDS_NODE_TYPE type) { _type = type; }
    void set( long id, NERDS_NODE_TYPE type) {
        _id = id;
        _type = type;
    }
    double lodging_cost() { return _lodging_cost; }
    LODGING_TYPE lodging_type() { return _lodging_type; }
    void set_lodging_cost( double cost) { _lodging_cost = cost; }   
    void set_lodging_type( LODGING_TYPE type) { _lodging_type = type; }
    vector<long> &fans() { return _fans; }
    long fans( long id) { // return arc_id entering/exiting this node;
        assert(id >= 0 && id < _fans.size());   
        return _fans[ id]; 
    }
    void add_fan( long arc_id) { // id is index of arc entering/exiting;
        _fans.push_back( arc_id);
    }
    vector<PAIR_TWO> &prob_intervals() { return _prob_intervals; }
    PAIR_TWO prob_intervals( long id) {
        //assert(id >= 0 && id < _prob_intervals.size());   
        return _prob_intervals[ id]; 
    }
    void add_prob_intervals( PAIR_TWO pt) {
        _prob_intervals.push_back( pt);
    }
    void set_prob_intervals( long index, long prob) {
        //assert(index >= 0 && index < _prob_intervals.size()); 
        _prob_intervals[ index].set_to_i( prob); // to_i is just storage;
    }
    bool processed_children_all() { 
        return (_children_count == _processed_children);
    }
    void incr_processed_children() { _processed_children ++; }
    void decr_processed_children() { _processed_children --; }
    void reset_processed_children() { _processed_children = 0; }
    void incr_children_count() { _children_count ++; }
    void decr_children_count() { _children_count --; }
    long children_count() const { return _children_count; }
    void reset_children_count() { _children_count = 0; }
};

////////////////////////////////////////////////////////////////////////////////
//
// NERDS_GRAPH
//
////////////////////////////////////////////////////////////////////////////////

class NERDS_GRAPH {
 private:
    long _arcs_count;
    long _sectionalizing_count;
    long _ties_count;
    long _nodes_count;
    long _feeders_count;
    long _super_feeder_id;
    vector<NERDS_NODE> _nodes; // feeder-nodes, bus-nodes;
    vector<NERDS_ARC> _arcs; // branches;
    vector<int> _ties; // ids or arcs/branches that are ties;
    // for each feeder store the maximum depth among all its leaves;
    vector<int> _feeders_depth;
    // grid size: number of columns and rows of the grid on which all nodes
    // are placed for gui display; these values are read in from the 
    // input network_file;
    long _nx,_ny;
    // store all leaves of each feeder for speed up of calculations;
    vector<vector<int> > _leaves;
    // storage and sketch stuff used to implement Baran Method 1,2;
    vector<int> _sketch_left_arc_ids; // left part of loop;
    vector<int> _sketch_right_arc_ids;
    vector<TUPLE_THREE> _sketch_exchange_losses; // all 1st ord exchanges, curr iter;
    vector<PAIR_TWO> _sketch_feeder_pairs; // stores feeder id pairs from -> to;

    // next are arrays used to "implement" a generic flow graph "indexing"
    // of all possible vertices: s,t,donor,acceptor, and switch;
    // we need these to be able to identify nodes in the mcmf flow graph;
    // _sketch_donors_id has _feeders_count entries; for example
    // _sketch_donors_id[j] == "index of this feeder node in the generic flow graph"
    // whose only a subgraph will constitute the actual flow graph built during
    // each iteration;
    // _sketch_acceptors_id has _feeders_count entries; for example
    // _sketch_acceptors_id[j] == "index of this feeder node in the generic flow graph"
    // _sketch_switches has _ties_count entries;
    // Important: indexing within the generic flow graph is:
    // 0 --> s_index
    // 1 --> t_index
    // 2,3,...,_feeders_count+2 --> donor indices
    // _feeders_count+3,...,2*_feeders_count+2 --> acceptor indices
    // 2*_feeders_count+3,... --> switch nodes indices
    vector<int> _sketch_donors_id;
    vector<int> _sketch_acceptors_id;
    vector<int> _sketch_switches_id;
    // next one is a storage of id's of all affected feeder nodes whose 
    // trees were changed due to branch exchanges; PowerFlow solution
    // based on random walks willfocus only on these trees only for 
    // re-calculation; when partial rw in fact we'll focus only on 
    // subsets of nodes;
    vector<int> _sketch_affected_feeders;

 public:
    INTERFACE_MCMF _interface_mcmf; // constructed locally;
    PRONERDS *_pronerds; // its host;

 public:
    NERDS_GRAPH() {}
    NERDS_GRAPH(long branch_count, long sectionalizing_count, long ties_count,
        long nodes_count, long feeders_count) {
        _arcs_count = branch_count;
        _sectionalizing_count = sectionalizing_count;
        _ties_count = ties_count;
        _nodes_count = nodes_count;
        _feeders_count = feeders_count;
        _super_feeder_id = nodes_count - 1;
        // build nodes and arcs arrays;
        build_nodes( _nodes_count);
        build_arcs( _arcs_count);
        _ties.resize( ties_count);
    }
    ~NERDS_GRAPH() {}

    void clear_nodes() {
        for (int i = 0; i < _nodes_count; i++) {
            _nodes[i].clear();
            _nodes[i].prob_intervals().clear();
        }
        for (int i = 0; i < _arcs_count; i++) {
            _arcs[i].clear();
        }
    }
    void clear_nodes_partially( vector<int> &nodes_for_rw_update) {
        // for super-fast partial RW based PowerFlow this is commented;
        //for (int i = _feeders_count; i < _nodes_count; i++) {
        //  _nodes[i].clear_partially();
        //}
        // clear only prob intervals of nodes directly affected
        // by branch exchanges, because only these nodes have now
        // changed the number of their adjacent neighbors, and not even
        // all of them;
        long count = nodes_for_rw_update.size();
        for ( long j = 0; j < count; j ++) {
            _nodes[ nodes_for_rw_update[j] ].prob_intervals().clear();
        }
    }
    INTERFACE_MCMF *interface_mcmf() { return &_interface_mcmf; }
    PRONERDS *pronerds() { return _pronerds; }  
    void set_host_pronerds(PRONERDS *host) { _pronerds = host; }
    long arcs_count() { return _arcs_count; }
    long sectionalizing_count() { return _sectionalizing_count; }
    long ties_count() { return _ties_count; }
    long nodes_count() { return _nodes_count; }
    long feeders_count() { return _feeders_count; }
    long nx() { return _nx; }
    long ny() { return _ny; }
    void set_grid_size_nx_ny( long nx, long ny) { 
        _nx = nx; _ny = ny;
    }
    vector<NERDS_NODE> &nodes() { return _nodes; }
    NERDS_NODE *get_node(long id) {
        assert(id >= 0 && id < _nodes_count);   
        return &_nodes[ id]; 
    }
    NERDS_ARC *get_arc(long id) {
        return &_arcs[ id]; 
    }
    void build_nodes(long nodes_count) {
        _nodes_count = nodes_count;
        _nodes.resize( nodes_count);
        for (int i = 0; i < _nodes_count; i++) { _nodes[i].set_id( i); }
    }
    void build_arcs(long arcs_count){ 
        _arcs_count = arcs_count;
        _arcs.resize( arcs_count);
        for (int i = 0; i < _arcs_count; i++) { _arcs[i].set_id( i); }
    }
    void set_node( long node_id, NERDS_NODE_TYPE node_type) {
        assert( node_id >= 0 && node_id < _nodes_count);
        _nodes[ node_id].set( node_id, node_type);
    }
    void set_node_coords( long node_id, int x_coord, int y_coord) {
        assert( node_id >= 0 && node_id < _nodes_count);
        _nodes[ node_id].set_coords( x_coord, y_coord);
    }
    void set_arc( long arc_id, long src_id, long des_id, NERDS_ARC_TYPE arc_type,
        double r, double x) {
        assert( arc_id >= 0 && arc_id < _arcs_count);
        assert( src_id >= 0 && src_id < _nodes_count);
        assert( des_id >= 0 && des_id < _nodes_count);
        // arcs are bidirectional arcs; that is because after
        // a branch exchange the parent-child relationship may be reversed
        // for an arc, which may be transfered from one feeder-tree to another;
        // that is why we actualy store arc ids in the list of neighbours (adjacency
        // list) of a node;
        _nodes[ src_id].add_fan( arc_id);
        _nodes[ des_id].add_fan( arc_id);
        // set arc too; this could be avoided by storing this info together with 
        // every fanout of every node;
        // the arc is stored as an object only once; it does not matter
        // which node is src and which is des;
        _arcs[ arc_id].set( src_id, des_id, arc_type, r, x);
    }
    void add_tie( long i, long arc_id) {
        assert( i >= 0 && i < _ties_count);
        assert( arc_id >= 0 && arc_id < _arcs_count);
        _ties[ i] = arc_id;
    }
    
    // routines to calculate power flow solution using Baran's method;
    bool run_DistFlow_solution_baran( bool simplified_distflow = false,
        bool use_super_feeder_node = false,
        bool used_during_partial_rw = false);
    bool run_super_feeder_DistFlow_solution_baran( bool simplified_distflow);
    
    void backward_mixed_with_forward_sweep( long iter, NERDS_NODE *feeder_node,
        bool simplified_distflow);
    void backward_sweep( long iter, NERDS_NODE *feeder_node, 
        bool simplified_distflow);
    void forward_sweep( long iter, NERDS_NODE *feeder_node,
        bool simplified_distflow);
    void backward_branch_equations( NERDS_NODE *u_node, NERDS_NODE *d_node,
        bool simplified_distflow,
        bool is_for_lateral = false);
    void forward_branch_equations( NERDS_NODE *u_node, NERDS_NODE *d_node,
        bool is_fork_node, bool simplified_distflow,
        bool is_for_lateral = false);

    // routines to calculate power flow solution using random-walks method;
    bool run_PowerFlow_solution_rw_version1();
    bool run_PowerFlow_solution_rw_version2();
    bool run_PowerFlow_solution_rw_update_subset_of_nodes( 
        vector<int> &nodes_for_rw_update);
    void convert_to_home_nodes_at_half_depth();
    void convert_to_home_all_nodes_but_2nd_order_neighbours(
         vector<int> &nodes_for_rw_update); 
    void get_leaf_nodes( NERDS_NODE *feeder_node, vector<int> &nodes_for_rw);
    double estimate_voltage_by_random_walks( long node_id);
    bool perform_random_walk_from_node( NERDS_NODE *node, 
        double &v, long &steps_count);
    inline int arc_id_for_walk( NERDS_NODE *node, int rand_int);
    void compute_probability_intervals( NERDS_NODE *feeder_node);
    void forward_sweep_rw( NERDS_NODE *feeder_node);
    void compute_probability_intervals_one_node( NERDS_NODE *node);
    void forward_branch_equations_rw( NERDS_NODE *u_node,
        NERDS_NODE *d_node, bool is_fork_node);
    void add_feeder_id_of_this_node_to_affected_feeders( int feeder_id) {
        assert (feeder_id >= 0 && feeder_id < _feeders_count);
        long count = _sketch_affected_feeders.size();
        for ( long i = 0; i < count; i ++) {
            if ( _sketch_affected_feeders[ i] == feeder_id) {
                return; // already added;
            }
        }
        // if here, then, this feeder_id was not added yet to the "affected";
        _sketch_affected_feeders.push_back( feeder_id);
    }
    double compute_sum_nodes_voltages( NERDS_NODE *feeder_node, 
        long &tree_nodes_count);
    double compute_total_losses_of_feeder_tree( NERDS_NODE *feeder_node);
    double compute_total_losses();
    void initial_backward_P_Q_aggregation_feeder_tree( NERDS_NODE *feeder_node);
    // Baran Method 1,2;
    bool run_reconfiguration_method1_baran( bool simplified_distflow,
        POWER_FLOW_METHOD power_flow_solution = BARAN_DIST_FLOW); // default;
    bool explore_branch_exchanges_baran1();
    void find_loop_as_left_right_arc_sets(int tie_i);
    void estimate_loss_reduction_dueto_branch_exchanges_baran1_first_order(int tie_i);
    void estimate_loss_reduction_dueto_branch_exchanges_baran1_entire_loop(int tie_i);
    // new methods;
    bool run_reconfiguration_method_new( bool simplified_distflow,
        POWER_FLOW_METHOD power_flow_solution = BARAN_DIST_FLOW); // default;
    int initialize_indexing_of_generic_flow_graph();
    void initialize_interface_mcmf(
        MCMF_ALGORITHM engine_type, int max_outside_node_id, int bound_nn);
    bool explore_branch_exchanges_mcmf();
    void populate_interface_mcmf_with_flow_graph_info();
    void estimate_loss_reduction_dueto_branch_exchanges_new_first_order(int tie_i);
    void estimate_loss_reduction_dueto_branch_exchanges_new_entire_loop(int tie_i);

    void update_parents_of_all_nodes();
    void update_parents_of_all_nodes_in_feeder_tree( NERDS_NODE *feeder_node);
    void testing_distflow( bool simplified_distflow);
    void testing_hitting_times( bool simplified_distflow);
    void print_tie_switches();
    void print_graph();
    void print_voltages();
    void print_losses();
    void print_left_right_loop(int tie_i);
    void print_sketch_exchange_losses();
    void compute_children_counts() {
        // of every node; should be called after graph has been completely built;
        // first reset/clear all counts;
        for ( long i = 0; i < _nodes_count; i++) {
            _nodes[ i].reset_children_count();
        }
        for ( long i = 0; i < _nodes_count; i++) {
            NERDS_NODE *node = &_nodes[ i];
            long neighbors_count = node->fans().size();
            for ( long k = 0; k < neighbors_count; k++) {
                NERDS_ARC *arc = &_arcs[ node->fans( k)];
                if ( arc->type() != TIE) {
                    node->incr_children_count();
                }
            }
            // if this is not a feeder node then decrement w/ 1 final count
            // of children for this node i to account for its parent (because
            // all its neighbors have been counted including its parent);
            //if ( i >= _feeders_count) { node->decr_children_count(); }
            // this has changed a little with the super feeder introduction;
            // we do it now for all nodes oincluding normal feeder nodes, except
            // the super feeder node;
            if ( i != _super_feeder_id) { node->decr_children_count(); }
        }
    }
    void reset_feeders_depth() {
        // should be called only after feeders_count has been set;
        _feeders_depth.clear();
        for ( long i = 0; i < _feeders_count; i++) {
            _feeders_depth.push_back( 0);
        }
    }
    void compute_feeders_depth();
    void print_feeders_depth();
    static int my_rand_int(int max_val) 
    {
        // Note: use alternativ NRG from nerds_utils.cpp;
        // return a random integer between [0..max_val];
        // r is a random float value in the range [0,1);
        double r = (double)rand() / ((double)(RAND_MAX)+(double)(1));
        r = r * (max_val + 1);
        return int( r);
    }
};

#endif

