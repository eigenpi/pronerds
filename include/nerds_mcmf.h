#ifndef _NERDS_MCMF_H_
#define _NERDS_MCMF_H_

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <vector>
#include <deque>

using namespace std;

enum MCMF_ALGORITHM { MCMF_ALGO_EDMONDS = 0, MCMF_ALGO_CS2 = 1 }; // CS2 is Goldberg's;
enum FLOW_GRAPH_VERTEX_TYPE { NO_TYPE_VERTEX, SWITCH_VERTEX,
    DONOR_FEEDER_VERTEX, ACCEPTOR_FEEDER_VERTEX,
    SOURCE_S_VERTEX, SINK_T_VERTEX};

class TUPLE_THREE;

// here are a few different implementations of algorithms to solve minimum cost
// maximum flow (MCMF) problem;

////////////////////////////////////////////////////////////////////////////////
//
// MCMF_EDMONDS
//
////////////////////////////////////////////////////////////////////////////////

// 1
// This is an adapted (i.e., ported to C++) version of the Edmonds-Karp algo
// from: http://shygypsy.com/tools/mcmf4.cpp
// The original authors of that implementation: Frank Chu, Igor Naverniouk

// Note: NN is not used anymore; I allocate dynamically only as much as
// is needed via _NN;
#define NN 128 // 128 is maximum number of vertices+1 supported;
#define CLEAR_MEM_1D(p,v,n) { for(int i=0; i<n; i++) p[i]=v; }
#define CLEAR_MEM_2D(p,v,n) { for(int i=0; i<n; i++) for(int j=0; j<n; j++) p[i][j]=v; }
    
#define HALF_INFINITY (INT_MAX/2)

class MCMF_EDMONDS
{
 private:

#define BUBL { \
 t = _q[i]; _q[i] = _q[j]; _q[j] = t; \
 t = _inq[_q[i]]; _inq[_q[i]] = _inq[_q[j]]; _inq[_q[j]] = t; }

// Dijkstra's using non-negative edge weights (cost + potential)
#define POTENTIAL(u,v) ( _d[u] + _pi[u] - _pi[v])

 private:
    int _num_vertices;
    int _qs;
    int _NN;
    // adjacency matrix (fill this up)
    int **_cap;
    // cost per unit of flow matrix (fill this up)
    int **_cost;
    // flow network and adjacency list
    int **_fnet, **_adj;
    int *_deg;
    // labelling function
    int *_pi;
    // Dijkstra's predecessor, depth and priority queue
    int *_par, *_d, *_q, *_inq;

 public:
    MCMF_EDMONDS() { _NN = 0; }
    MCMF_EDMONDS(int num_vertices) {
        _NN = 0;
        _num_vertices = num_vertices;
    }
    ~MCMF_EDMONDS() {}

    void clear_all() {
        if ( _NN > 0) {
            CLEAR_MEM_2D( _cap, 0, _NN);
            CLEAR_MEM_2D( _cost, 0, _NN);
            CLEAR_MEM_2D( _fnet, 0, _NN);
            CLEAR_MEM_2D( _adj, 0, _NN);
            CLEAR_MEM_1D( _deg, 0, _NN);
            CLEAR_MEM_1D( _pi, 0, _NN);
            CLEAR_MEM_1D( _par, 0, _NN);
            CLEAR_MEM_1D( _d, 0, _NN);
            CLEAR_MEM_1D( _q, 0, _NN);
            CLEAR_MEM_1D( _inq, 0, _NN);
        }
    }
    void set_NN_and_allocate_arrays( int nn) {
        _NN = nn;
        _cap = (int**) malloc(_NN * sizeof(int *));
        _cost = (int**) malloc(_NN * sizeof(int *));
        _fnet = (int**) malloc(_NN * sizeof(int *));
        _adj = (int**) malloc(_NN * sizeof(int *));
        if ( _cap == NULL) { printf("\nError: malloc 1\n"); exit(1); }
        if ( _cost == NULL) { printf("\nError: malloc 2\n"); exit(1); }
        if ( _fnet == NULL) { printf("\nError: malloc 3\n"); exit(1); }
        if ( _adj == NULL) { printf("\nError: malloc 4\n"); exit(1); }
        for ( int i = 0; i < _NN; i++) {
            _cap[i] = (int *) malloc(_NN * sizeof(int));
            _cost[i] = (int *) malloc(_NN * sizeof(int));
            _fnet[i] = (int *) malloc(_NN * sizeof(int));
            _adj[i] = (int *) malloc(_NN * sizeof(int));
            if ( _cap[i] == NULL) { printf("\nError: malloc 5\n"); exit(1); }
            if ( _cost[i] == NULL) { printf("\nError: malloc 6\n"); exit(1); }
            if ( _fnet[i] == NULL) { printf("\nError: malloc 7\n"); exit(1); }
            if ( _adj[i] == NULL) { printf("\nError: malloc 8\n"); exit(1); }
        }
        _deg = (int*) malloc(_NN * sizeof(int));
        _pi = (int*) malloc(_NN * sizeof(int));
        _par = (int*) malloc(_NN * sizeof(int));
        _d = (int*) malloc(_NN * sizeof(int));
        _q = (int*) malloc(_NN * sizeof(int));
        _inq = (int*) malloc(_NN * sizeof(int));
        if ( _deg == NULL) { printf("\nError: malloc 9\n"); exit(1); }
        if ( _pi == NULL) { printf("\nError: malloc 10\n"); exit(1); }
        if ( _par == NULL) { printf("\nError: malloc 11\n"); exit(1); }
        if ( _d == NULL) { printf("\nError: malloc 12\n"); exit(1); }
        if ( _q == NULL) { printf("\nError: malloc 13\n"); exit(1); }
        if ( _inq == NULL) { printf("\nError: malloc 14\n"); exit(1); }
        //printf ("\nAllocated arrays of MCMF solver.");

        // reset all allocated mem;
        clear_all();
    }
    // Note: the CS2 engine has deallocate_arrays() counterpart;
    void free_arrays() {
        for ( int i = 0; i < _NN; i++) {
            free( _cap[i]); free( _cost[i]); free( _fnet[i]); free( _adj[i]);
        }
        free( _cap); free( _cost); free( _fnet); free( _adj);
        free( _deg); free( _pi); free( _par); free( _d); free( _q); free( _inq);
    }
    void set_num_vertices( int num_vertices) { _num_vertices = num_vertices; }
    void set_edge(int i, int j, int edge_capacity, int edge_cost) {
        assert( i >= 0 && i < _num_vertices);
        assert( j >= 0 && j < _num_vertices);
        _cap[i][j] = edge_capacity; //_cap[j][i] = edge_capacity;
        _cost[i][j] = edge_cost; //_cost[j][i] = edge_cost;
    }
    bool dijkstra( int n, int s, int t);
    int run_edmonds( int n, int s, int t, int &flow_cost);
    inline int potential(int u, int v) {
        return ( _d[u] + _pi[u] - _pi[v]);
    }
    void buble(int i, int j, int &t) {
        t = _q[i]; _q[i] = _q[j]; _q[j] = t;
        t = _inq[_q[i]]; _inq[_q[i]] = _inq[_q[j]]; _inq[_q[j]] = t;
    }
    int get_fnet(int u, int v) const { return _fnet[u][v]; }
    int get_cost(int u, int v) const { return _cost[u][v]; }
};


////////////////////////////////////////////////////////////////////////////////
//
// MCMF_CS2
//
////////////////////////////////////////////////////////////////////////////////

// 2
// This is an adapted (i.e., ported to C++) version of the Goldberg CS2 algo;
// CS2 is the second version of scaling algorithm for minimum-cost
// flow problems. For more detailed description, see "An
// Efficient Implementation of a Scaling Minimum-Cost Flow Algorithm" by
// A.V. Goldberg, J. Algorithms, Vol. 22 (1997), pp. 1--29.

#define MAX_64 (0x7fffffffffffffffLL)
#define MAX_32 (0x7fffffff)
#define PRICE_MAX MAX_64

#define MAXLINE       100 // max line length in the input file
#define ARC_FIELDS      5 // no of fields in arc line
#define NODE_FIELDS     2 // no of fields in node line
#define P_FIELDS        3 // /* no of fields in problem line
#define PROBLEM_TYPE "min" //  name of problem type

#define UNFEASIBLE          2
#define ALLOCATION_FAULT    5
#define PRICE_OFL           6

// parameters
#define UPDT_FREQ      0.4
#define UPDT_FREQ_S    30
#define SCALE_DEFAULT  12.0
// PRICE_OUT_START may not be less than 1
#define PRICE_OUT_START  1
#define CUT_OFF_POWER    0.44
#define CUT_OFF_COEF     1.5
#define CUT_OFF_POWER2   0.75
#define CUT_OFF_COEF2    1
#define CUT_OFF_GAP      0.8
#define CUT_OFF_MIN      12
#define CUT_OFF_INCREASE 4

#define TIME_FOR_PRICE_IN1    2
#define TIME_FOR_PRICE_IN2    4
#define TIME_FOR_PRICE_IN3    6

#define MAX_CYCLES_CANCELLED 0
#define START_CYCLE_CANCEL   100

#define MAX( x, y ) ( ( (x) > (y) ) ?  x : y )
#define MIN( x, y ) ( ( (x) < (y) ) ? x : y )
#define ABS( x ) ( (x) >= 0 ) ? (x) : -(x)
#define N_NODE( i ) ( ( (i) == NULL ) ? -1 : ( (i) - _nodes + _node_min ) )
#define N_ARC( a ) ( ( (a) == NULL )? -1 : (a) - _arcs )

class MCMF_CS2
{
 public:
    typedef long long int excess_t;
    typedef long long int price_t;

    class NODE;
    
    class ARC {
    public:
        long _rez_capacity; // residual capacity;
        price_t _cost; // cost of arc;
        NODE *_head; // head node;
        ARC *_sister; // opposite arc;
        // default values;
        long _rez_capacity_default;
        price_t _cost_default;
    public:
        ARC() {}
        ~ARC() {}

        void set_rez_capacity( long rez_capacity) { _rez_capacity = rez_capacity; }
        void dec_rez_capacity( long delta) { _rez_capacity -= delta; }
        void inc_rez_capacity( long delta) { _rez_capacity += delta; }
        void set_cost( price_t cost) { _cost = cost; }
        void multiply_cost( price_t mult) { _cost *= mult; }
        void set_head( NODE *head) { _head = head; }
        void set_sister( ARC *sister) { _sister = sister; }
        long rez_capacity() { return _rez_capacity; }
        price_t cost() { return _cost; }
        NODE *head() { return _head; }
        ARC *sister() { return _sister; }
        // defaults;
        void set_rez_capacity_default( long r_cap) { _rez_capacity_default = r_cap; }
        void set_cost_default( price_t c) { _cost_default = c; }
        long rez_capacity_default() { return _rez_capacity_default; }
        price_t cost_default() { return _cost_default; }
        void reset_default() {
            _rez_capacity = _rez_capacity_default;
            _cost = _cost_default;
        }
    };

    class NODE {
    public:
        excess_t _excess; // excess of the node;
        price_t _price; // distance from a sink;
        ARC *_first; // first outgoing arc;
        ARC *_current; // current outgoing arc;
        ARC *_suspended;
        NODE *_q_next; // next node in push queue
        NODE *_b_next; // next node in bucket-list;
        NODE *_b_prev; // previous node in bucket-list;
        long _rank; // bucket number;
        long _inp; // auxilary field;
    public:
        NODE() {}
        ~NODE() {}

        void set_excess( excess_t excess) { _excess = excess; }
        void dec_excess( long delta) { _excess -= delta; }
        void inc_excess( long delta) { _excess += delta; }
        void set_price( price_t price) { _price = price; }
        void dec_price( long delta) { _price -= delta; }
        void set_first( ARC *first) { _first = first; }
        void set_current( ARC *current) { _current = current; }
        void inc_current() { _current ++; }
        void set_suspended( ARC *suspended) { _suspended = suspended; }
        void set_q_next( NODE *q_next) { _q_next = q_next; }
        void set_b_next( NODE *b_next) { _b_next = b_next; }
        void set_b_prev( NODE *b_prev) { _b_prev = b_prev; }
        void set_rank( long rank) { _rank = rank; }
        void set_inp( long inp) { _inp = inp; }
        excess_t excess() { return _excess; }
        price_t price() { return _price; }
        ARC *first() { return _first; }
        void dec_first() { _first --; }
        void inc_first() { _first ++; }
        ARC *current() { return _current; }
        ARC *suspended() { return _suspended; }
        NODE *q_next() { return _q_next; }
        NODE *b_next() { return _b_next; }
        NODE *b_prev() { return _b_prev; }
        long rank() { return _rank; }
        long inp() { return _inp; }
        void reset() { _excess = 0;  _price = 0; _rank = 0; _inp = 0; }
    };
 
    class BUCKET {
    private:
        // 1st node with positive excess or simply 1st node in the buket;
        NODE *_p_first;
    public:
        BUCKET( NODE *p_first) : _p_first(p_first) {}
        BUCKET() {}
        ~BUCKET() {}

        void set_p_first( NODE *p_first) { _p_first = p_first; }
        NODE *p_first() { return _p_first; }
    };

 public:
    long _n; // number of nodes
    long _m; // number of arcs

    long *_cap; // array containig capacities
    NODE *_nodes; // array of nodes
    NODE *_sentinel_node; // next after last
    NODE *_excq_first; // first node in push-queue
    NODE *_excq_last; // last node in push-queue
    ARC *_arcs; // array of arcs
    ARC *_sentinel_arc; // next after last

    BUCKET *_buckets; // array of buckets
    BUCKET *_l_bucket; // last bucket
    long _linf; // number of l_bucket + 1
    int _time_for_price_in;

    price_t _epsilon; // quality measure
    price_t _dn; // cost multiplier = number of nodes  + 1
    price_t _price_min; // lowest bound for prices
    price_t _mmc; // multiplied maximal cost
    double _f_scale; // scale factor
    double _cut_off_factor; // multiplier to produce cut_on and cut_off from n and epsilon
    double _cut_on; // the bound for returning suspended arcs
    double _cut_off; // the bound for suspending arcs
    excess_t _total_excess; // total excess

    // if = 1 - signal to start price-in ASAP - 
    // maybe there is infeasibility because of suspended arcs 
    int _flag_price;
    // if = 1 - update failed some sources are unreachable: either the 
    // problem is unfeasible or you have to return suspended arcs
    int _flag_updt;
    // maximal number of cycles cancelled during price refine 
    int _snc_max;

    // dummy variables;
    ARC _d_arc; // dummy arc - for technical reasons
    NODE _d_node; // dummy node - for technical reasons
    NODE *_dummy_node; //  the address of d_node
    NODE *_dnode;

    long _n_rel; // number of relabels from last price update
    long _n_ref; // current number of refines
    long _n_src; // current number of nodes with excess
    long _n_push;
    long _n_relabel;
    long _n_discharge;
    long _n_refine;
    long _n_update;
    long _n_scan;
    long _n_prscan;
    long _n_prscan1;
    long _n_prscan2;
    long _n_bad_pricein;
    long _n_bad_relabel;
    long _n_prefine;

    bool _no_zero_cycles; // finds an optimal flow with no zero-cost cycles
    bool _check_solution; // check feasibility/optimality. HIGH OVERHEAD!
    bool _comp_duals; // compute prices?
    bool _cost_restart; // to be able to restart after a cost function change
    bool _print_ans;
    long long int *_node_balance;

    // sketch variables used during reading in arcs;
    long _node_min; // minimal no of nodes
    long _node_max; // maximal no of nodes
    long *_arc_first; // internal array for holding
                      //   - node degree
                      //   - position of the first outgoing arc
    long *_arc_tail; // internal array: tails of the arcs
    long _pos_current;
    ARC *_arc_current;
    ARC *_arc_new;
    ARC *_arc_tmp;
    price_t _max_cost; // maximum cost
    excess_t _total_p; // total supply
    excess_t _total_n; // total demand
    // pointers to the node structure
    NODE *_i_node;
    NODE *_j_node;
    // final solution, if any;
    double _flow_cost;
    int _max_flow;
    // next three variables are used only for iterative lowering of demand
    // until a solution is found, if any;
    long _source_id; // S;
    long _sink_id; // T;
    int _demand_ST;

 public:
    MCMF_CS2( long num_nodes, long num_arcs) {
        _n = num_nodes;
        _m = num_arcs;

        _flag_price = 0; _flag_updt = 0; 
        _n_push = 0; _n_relabel = 0; _n_discharge = 0; 
        _n_ref = 0; _n_refine = 0; _n_update = 0; 
        _n_scan = 0; _n_prscan = 0; _n_prscan1 = 0; _n_prscan2 = 0; 
        _n_bad_pricein = 0; _n_bad_relabel = 0; _n_prefine = 0;
        _no_zero_cycles = false; _check_solution = false;
        _comp_duals = false; _cost_restart = false;
        _print_ans = true;

        _cut_on = 0.0;  _cut_off = 0.0;
        _total_excess = 0; _n_src = 0; _n_rel = 0;
        _time_for_price_in = 0;
        _total_p = 0; _total_n = 0;
        _snc_max = 0;

        // allocate arrays and prepare for "receiving" arcs;
        // will also reset _pos_current, etc.;
        allocate_arrays();
    }
    MCMF_CS2() {
        // if this ctor is used, then immediately after, we'll have to set num
        // of nodes and arcs and call the arrays allocations routine; only after 
        // that, we'll be prepared for "receiving" arcs;
        clear_all();
    }
    ~MCMF_CS2() {}


    void clear_all() {
        _n = 0;
        _m = 0;

        _flag_price = 0; _flag_updt = 0; 
        _n_push = 0; _n_relabel = 0; _n_discharge = 0; 
        _n_ref = 0; _n_refine = 0; _n_update = 0; 
        _n_scan = 0; _n_prscan = 0; _n_prscan1 = 0; _n_prscan2 = 0; 
        _n_bad_pricein = 0; _n_bad_relabel = 0; _n_prefine = 0;
        _no_zero_cycles = false; _check_solution = false;
        _comp_duals = false; _cost_restart = false; 
        _print_ans = false;

        _cut_on = 0.0;  _cut_off = 0.0;
        _total_excess = 0; _n_src = 0; _n_rel = 0;
        _time_for_price_in = 0;
        _total_p = 0; _total_n = 0;
        _snc_max = 0;
    }
    void clear_all_light() {
        // should be used as the first call in the sequence of re-runing the 
        // cs2 algo; used only during iterative lowering of demand for ST 
        // in order to find a solution;
        _flag_updt = 0; 
        _n_push = 0; _n_relabel = 0; _n_discharge = 0; 
        _n_refine = 0; _n_update = 0; 
        _n_scan = 0; _n_prscan = 0; _n_prscan1 = 0; _n_prscan2 = 0; 
        _n_bad_pricein = 0; _n_bad_relabel = 0; _n_prefine = 0;

        _cut_on = 0.0;  _cut_off = 0.0;
        _total_excess = 0; _n_src = 0; _n_rel = 0;
        _time_for_price_in = 0;
        _total_p = 0; _total_n = 0;
        _snc_max = 0;

        // clear up nodes;
        for ( NODE *i = _nodes; i <= _nodes + _n; i ++ ) {
            i->reset();
        }
        for ( ARC *a = _arcs; a <= _arcs + _m; a ++ ) {
            a->reset_default();
        }
    }

    // Note: the Edmonds engine has free_arrays() counterpart;
    void deallocate_arrays() {
        if ( _arcs) free ( _arcs );
        if ( _dnode) delete _dnode;
        if ( _cap) free ( _cap );
        if ( _buckets) free ( _buckets );
        if ( _check_solution == true) free ( _node_balance );
        if ( _nodes) {
            _nodes = _nodes - _node_min;
            free ( _nodes );
        }

        // these two used to be inside pre_processing(); I brought them here
        // to have all "free's" clustered in one place;
        if ( _arc_first ) free ( _arc_first ); 
        if ( _arc_tail ) free ( _arc_tail );
    }
    // Note: inside MCMF_CS2 we work with nodes and arcs; inside MCMF_EDMONDS
    // we work with vertices and edges; inside the interface we use vertices 
    // and edges; this should eliminate any confusion;
    void set_num_nodes( int num_nodes) { _n = num_nodes; }
    void set_num_arcs( int num_arcs) { _m = num_arcs; }
    void err_end( int cc);
    void allocate_arrays();
    void set_arc( long tail_node_id, long head_node_id,
        long low_bound, long up_bound, price_t cost);
    void set_supply_demand_of_node( long id, long excess);
    void pre_processing();
    void cs2_initialize();
    void cs2_initialize_light();
    void up_node_scan( NODE *i);
    void price_update();
    int relabel( NODE *i);
    void discharge( NODE *i);
    int price_in();
    // I made refine() now to return an int; it's possible that I have source S
    // with say 6 outgoing arcs and sink T with 6 incomming arcs but the only
    // possible solution is only 5; in this case if I set demand of S and T as 6,-6
    // CS2 would return Error 2 (unfeasible solution);
    int refine();
    int price_refine();
    void compute_prices();
    void price_out();
    int update_epsilon();
    int check_feas();
    int check_cs();
    int check_eps_opt();
    void print_solution( bool print_flow_only = false);
    void print_graph();
    void print_stuff();
    int get_flow_of_arc( int src_id, int des_id);   
    void finishup( double *objective_cost);
    int cs2( double *objective_cost);
    int run_cs2(); // part 1 + part 2;
    bool decrement_demand_of_source_and_sink();
    
    void set_source_id( long id) { _source_id = id; }
    void set_sink_id( long id) { _sink_id = id; }
    void set_demand_ST( int d) { _demand_ST = d; }


    // shared utils;
    void increase_flow( NODE *i, NODE *j, ARC *a, long df) {
        i->dec_excess( df);
        j->inc_excess( df);
        a->dec_rez_capacity( df);
        a->sister()->inc_rez_capacity( df);
    }
    bool time_for_update() { 
        return ( _n_rel > _n * UPDT_FREQ + _n_src * UPDT_FREQ_S); 
    }
    // utils for excess queue;
    void reset_excess_q() {
        for ( ; _excq_first != NULL; _excq_first = _excq_last ) {
            _excq_last = _excq_first->q_next();
            _excq_first->set_q_next( _sentinel_node);
        }
    }
    bool out_of_excess_q( NODE *i) { return ( i->q_next() == _sentinel_node); }
    bool empty_excess_q() { return ( _excq_first == NULL); }
    bool nonempty_excess_q() { return ( _excq_first != NULL); }
    void insert_to_excess_q( NODE *i) {
        if ( nonempty_excess_q() ) {
            _excq_last->set_q_next( i);
        } else {
            _excq_first = i;
        }
        i->set_q_next( NULL);
        _excq_last = i;
    }
    void insert_to_front_excess_q( NODE *i) {
        if ( empty_excess_q() ) {
            _excq_last = i;
        }
        i->set_q_next( _excq_first);
        _excq_first = i;
    }
    void remove_from_excess_q( NODE *i) {
        i = _excq_first;
        _excq_first = i->q_next();
        i->set_q_next( _sentinel_node);
    }
    // utils for excess queue as a stack;
    bool empty_stackq() { return empty_excess_q(); }
    bool nonempty_stackq() { return nonempty_excess_q(); }
    void reset_stackq() { reset_excess_q(); }
    void stackq_push( NODE *i) {
        i->set_q_next( _excq_first);
        _excq_first = i;
    }
    void stackq_pop( NODE *i) {
        remove_from_excess_q( i);
    }
    // utils for buckets;
    void reset_bucket( BUCKET *b) { b->set_p_first( _dnode); }
    bool nonempty_bucket( BUCKET *b) { return ( (b->p_first()) != _dnode); }
    void insert_to_bucket( NODE *i, BUCKET *b) {
        i->set_b_next( b->p_first() );
        b->p_first()->set_b_prev( i);
        b->set_p_first( i);
    }
    void get_from_bucket( NODE *i, BUCKET *b) {
        i = b->p_first();
        b->set_p_first( i->b_next());
    }
    void remove_from_bucket( NODE *i, BUCKET *b) {
        if ( i == b->p_first() ) {
            b->set_p_first( i->b_next());
        } else {
            i->b_prev()->set_b_next( i->b_next());
            i->b_next()->set_b_prev( i->b_prev());
        }
    }
    // misc utils;
    void update_cut_off() {
        if ( _n_bad_pricein + _n_bad_relabel == 0) {
            _cut_off_factor = CUT_OFF_COEF2 * pow( (double)_n, CUT_OFF_POWER2 );
            _cut_off_factor = MAX ( _cut_off_factor, CUT_OFF_MIN );
            _cut_off = _cut_off_factor * _epsilon;
            _cut_on = _cut_off * CUT_OFF_GAP;
        } else {
            _cut_off_factor *= CUT_OFF_INCREASE;
            _cut_off = _cut_off_factor * _epsilon;
            _cut_on = _cut_off * CUT_OFF_GAP;
        }
    }
    void exchange( ARC *a, ARC *b) {
        if ( a != b ) {
            ARC *sa = a->sister();
            ARC *sb = b->sister();
            long d_cap;                     

            _d_arc.set_rez_capacity( a->rez_capacity());
            _d_arc.set_cost( a->cost());
            _d_arc.set_rez_capacity_default( a->rez_capacity_default());
            _d_arc.set_cost_default( a->cost_default());
            _d_arc.set_head( a->head());

            a->set_rez_capacity( b->rez_capacity());
            a->set_cost( b->cost());
            a->set_rez_capacity_default( b->rez_capacity_default());
            a->set_cost_default( b->cost_default());
            a->set_head( b->head());

            b->set_rez_capacity( _d_arc.rez_capacity());
            b->set_cost( _d_arc.cost());
            b->set_rez_capacity_default( _d_arc.rez_capacity_default());
            b->set_cost_default( _d_arc.cost_default());
            b->set_head( _d_arc.head());

            if ( a != sb ) {            
                b->set_sister( sa);
                a->set_sister( sb);
                sa->set_sister( b);
                sb->set_sister( a);
            }
                        
            d_cap = _cap[ a - _arcs];
            _cap[ a - _arcs] = _cap[ b - _arcs];
            _cap[ b - _arcs] = d_cap;   
        }
    }
};

////////////////////////////////////////////////////////////////////////////////
//
// INTERFACE_MCMF
//
////////////////////////////////////////////////////////////////////////////////

// create an interface/wrapper around any mcmf engine; need it so that all
// nodes in the graph that the mcmf engine works with are indexed continously;
// use mapping vectors to translate/interface;

class INTERFACE_MCMF
{
 private:
    // outside world works with nodes, internally inside this wrapper we work
    // with vertices; nodes outside may be 0...17 but only some will be added
    // as vertices to build up the flow graph, let's say 9 of them; these will
    // be indexed internally bteween 0...8; we need _max_outside_node_id as 18
    // in this example to allocate enough mem for the mapping between outside
    // and inside worlds;
    int _max_outside_node_id;
    // _num_vertices, _num_edges will be incremented as nodes are added (they become
    // vertices of the flow graph); at the end they will store the total number
    // of vertices, edges;
    int _num_vertices;
    int _num_edges;
    // number of donor,acceptor nodes; used only for setting up the cs2 mcmf engine;
    int _num_vertices_donors;
    int _num_vertices_acceptors;
    int _flow_cost;
    int _max_flow;
    bool _is_empty; // flags that there is no flow graph, no mcmf problem;
    // source s and sink t vertex indices of the flow graph; will be set from
    // the outside world;
    int _s_index, _t_index;
    // if _vertex_already_present[i] = j >= 0, then vertex with name/id "i" has
    // been already added; its index in the _vertices array is "j";
    vector<int> _vertex_already_present;
    // next storage arrays store info on vertices; using them avoids
    // introducing a new vertex class for now;
    vector<FLOW_GRAPH_VERTEX_TYPE> _vertex_type;
    vector<int> _vertex_brex_index; // indices of branch exchanges as tags for vertices;
    // we do not have a _edges_capacity too because every edge will have capacity 1;
    vector<int> _edges_src_vertex;
    vector<int> _edges_des_vertex;
    vector<int> _edges_cost;

 public:
    // these object engines will be reused; they need clear_all() functions;
    // only one engine is used; here we just have more implementations of 
    // various engines;
    MCMF_EDMONDS _engine_mcmf_edmonds;
    int _edmonds_bound_NN;
    MCMF_CS2 _engine_mcmf_cs2;

 public:
    INTERFACE_MCMF( int max_outside_node_id) {
        _edmonds_bound_NN = 0;
        _max_outside_node_id = max_outside_node_id;
        _vertex_already_present.resize( _max_outside_node_id);
        clear_all();
    }
    INTERFACE_MCMF() {}
    ~INTERFACE_MCMF() {}

    MCMF_EDMONDS &engine_mcmf_edmonds() { return _engine_mcmf_edmonds; }
    MCMF_CS2 &engine_mcmf_cs2() { return _engine_mcmf_cs2; }
    void set_edmonds_bound_NN( int bound_nn) {
        // set the size of arrays to be allocated for the Edmunds solver
        // and allocate mem;
        _edmonds_bound_NN = bound_nn; 
        _engine_mcmf_edmonds.set_NN_and_allocate_arrays( bound_nn);
    }
    // if default ctor is used to create objects of this type, then we'll have
    // to set its _max_outside_node_id and allocate memory;
    void set_max_outside_node_id( int max_outside_node_id) {
        _max_outside_node_id = max_outside_node_id;
        _vertex_already_present.resize( _max_outside_node_id);
        clear_all();
    }
    void inc_num_vertices_donors() { _num_vertices_donors ++; }
    int num_vertices_donors() const { return _num_vertices_donors; }
    void inc_num_vertices_acceptors() { _num_vertices_acceptors ++; }
    int num_vertices_acceptors() const { return _num_vertices_acceptors; }
    void set_s_index(int id) { _s_index = id; }
    void set_t_index(int id) { _t_index = id; }
    void set_is_empty(int b) { _is_empty = b; }
    bool is_empty() { return _is_empty; }
    // clear_all() is called during each main iteration, right after power flow
    // solutions is recalculated;
    void clear_all() {
        _is_empty = true;
        _s_index = -1; _t_index = -1;
        _num_vertices = 0;
        _num_edges = 0;
        _num_vertices_donors = 0;
        _num_vertices_acceptors = 0;
        for ( int i = 0; i < _max_outside_node_id; i++) {
            _vertex_already_present[ i] = -1; // i.e., vertex with id = i not present;
        }
        _vertex_type.clear();
        _vertex_brex_index.clear();
        _edges_src_vertex.clear();
        _edges_des_vertex.clear();      
        _edges_cost.clear();
        // clear the engines as well;
        _engine_mcmf_edmonds.clear_all();
        _engine_mcmf_cs2.clear_all();
    }
    void free_arrays_mcmf_edmonds() {
        _engine_mcmf_edmonds.free_arrays();
    }
    void deallocate_arrays_mcmf_cs2() {
        _engine_mcmf_cs2.deallocate_arrays();
    }
    void set_s_t() {
        // call only after interface has been populated with arcs,
        // created vertices also;
        _s_index = _vertex_already_present[ 0];
        _t_index = _vertex_already_present[ 1];
    }
    bool is_node_added( int node_id) {
        return ( _vertex_already_present[ node_id] >= 0);
    }
    bool add_node( int node_id,
        FLOW_GRAPH_VERTEX_TYPE type, int index_in_sketch) {
        // node_id in outside world becomes vertex index in inside world;
        // return the correspondix vertex index;
        if ( _vertex_already_present[ node_id] >= 0) {
            // been added already;
            return false;
        } else { // add it;
            _vertex_already_present[ node_id] = _num_vertices;
            _vertex_type.push_back( type);
            _vertex_brex_index.push_back( index_in_sketch);
            _num_vertices ++;
            return true;
        }
    }
    void add_arc( int src_node_id, int des_node_id, int cost,
        FLOW_GRAPH_VERTEX_TYPE src_type, int src_index_in_sketch,
        FLOW_GRAPH_VERTEX_TYPE des_type, int des_index_in_sketch) 
    {
        // src/des_index_in_sketch will store index of the br exchange 
        // inside _sketch_exchange_losses; this br ex is "attached" to this vertex;
        // src node;
        add_node( src_node_id, src_type, src_index_in_sketch);
        // des node;
        add_node( des_node_id, des_type, des_index_in_sketch);
        // new edge;
        int temp_vertex_id = _vertex_already_present[ src_node_id];
        _edges_src_vertex.push_back( temp_vertex_id); // edge index is _num_edges;
        temp_vertex_id = _vertex_already_present[ des_node_id];
        _edges_des_vertex.push_back( temp_vertex_id); // edge index is _num_edges;
        _edges_cost.push_back( cost);
        _num_edges ++;
    }
    // main launching point for the mcmf problem solution;
    void run_mcmf( MCMF_ALGORITHM engine_type) {
        if ( _is_empty) {
            printf("\n---> MCMF flow graph is empty! No MCMF problem!");
            return; // no mcmf problem;
        }
        if ( engine_type == MCMF_ALGO_EDMONDS) {

            // (a) construct the graph inside the engine;
            build_flow_graph_mcmf_edmonds();
            // (b) solve the mcmf problem;
            _max_flow = _engine_mcmf_edmonds.run_edmonds(
                _num_vertices, _s_index, _t_index, _flow_cost);
        }
        else if ( engine_type == MCMF_ALGO_CS2) {

            // (a) construct the graph inside the engine;
            build_flow_graph_mcmf_cs2();
            // (b) solve the mcmf problem;
            _engine_mcmf_cs2.run_cs2();
        }
    }

    // engine preparations;
    // () use MCMF_EDMONDS engine;
    void build_flow_graph_mcmf_edmonds() {
        // this has to be called only after all nodes and arcs have been
        // input to this interface;
        _engine_mcmf_edmonds.set_num_vertices( _num_vertices);
        // add all edges to the mcmf_edmonds engine; takes care of vertices
        // as well as it works with matrices;
        for ( int j = 0; j < _num_edges; j++) {
            _engine_mcmf_edmonds.set_edge(
                _edges_src_vertex[j], _edges_des_vertex[j],
                1, _edges_cost[j]); // capacity 1, cost;
        }
        // in the outside world S,T are indexed generically as 0,1;
        // inside the interface though they will have id's based on 
        // the order in which arcs were added to the interface;
        // so, we still have to now set _s_index and _t_index for
        // the inside world of interface;
        _s_index = _vertex_already_present[ 0];
        _t_index = _vertex_already_present[ 1];
    }
    // () use MCMF_CS2 engine;
    void build_flow_graph_mcmf_cs2() {
        // this has to be called only after all nodes and arcs have been
        // input to this interface;
        _engine_mcmf_cs2.set_num_nodes( _num_vertices);
        _engine_mcmf_cs2.set_num_arcs( _num_edges);
        _engine_mcmf_cs2.allocate_arrays();

        // add all edges to the mcmf_cs2 engine; takes care of vertices
        // as well as it works with matrices;
        for ( int j = 0; j < _num_edges; j++) {
            _engine_mcmf_cs2.set_arc(
                _edges_src_vertex[j], _edges_des_vertex[j],
                0, 1, _edges_cost[j]); // lower capacity 0, capacity 1, cost;
        }
        // in the outside world S,T are indexed generically as 0,1;
        // inside the interface though they will have id's based on 
        // the order in which arcs were added to the interface;
        // so, we still have to now set _s_index and _t_index for
        // the inside world of interface;
        _s_index = _vertex_already_present[ 0];
        _t_index = _vertex_already_present[ 1];
        _engine_mcmf_cs2.set_source_id( _s_index);
        _engine_mcmf_cs2.set_sink_id( _t_index);
        // add the supply and demands of S and T; take the min between
        // the arcs leaving S and the arcs arriving to T (otherwise CS2 algo
        // does not like it and errors out);
        int min_demand = min( _num_vertices_donors, _num_vertices_acceptors);
        _engine_mcmf_cs2.set_demand_ST( min_demand);
        _engine_mcmf_cs2.set_supply_demand_of_node( _s_index, min_demand);
        _engine_mcmf_cs2.set_supply_demand_of_node( _t_index, -min_demand);
        //printf("\n Demand at S and T: %d, %d", min_demand, -min_demand);
    }
    void get_branch_exhanges_to_be_realized( vector<int> &sketch_exchange_losses_ids,
        MCMF_ALGORITHM engine_type);
 
    // debug;
    void print_interface() {
        printf("\nInterface mcmf graph. Source S:%d Sink T:%d", _s_index, _t_index); 
        printf("\nVertices: %d Edges: %d", _num_vertices, _num_edges);
        printf("\nmax_outside_node_id: %d", _max_outside_node_id);
        for ( int i = 0; i < _num_edges; i ++) {
            printf("\n %d --> %d  cost:%d",
                _edges_src_vertex[i],
                _edges_des_vertex[i],
                _edges_cost[i]);
        }
    }
    void print_flow( MCMF_ALGORITHM engine_type) {
        if ( engine_type == MCMF_ALGO_EDMONDS) {
            printf("\n flow: %d", _max_flow);
            printf("\n cost: %d\n", _flow_cost);
            for ( int i = 0; i < _num_vertices; i++) {
                for ( int j = 0; j < _num_vertices; j++) {
                    int this_flow = _engine_mcmf_edmonds.get_fnet(i,j);         
                    if ( this_flow)
                        printf("%d -> %d flow:%d cost:%d\n",i,j,
                               _engine_mcmf_edmonds.get_fnet(i,j),
                               _engine_mcmf_edmonds.get_cost(i,j));
                }
            }
        } else if ( engine_type == MCMF_ALGO_CS2) {
            _engine_mcmf_cs2.print_solution( true); // print_flow_only;
        }       
    }
};

#endif
