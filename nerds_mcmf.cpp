#include "nerds.h"
#include "nerds_mcmf.h"

#include <assert.h>
#include <vector>
#include <deque>

using namespace std;


////////////////////////////////////////////////////////////////////////////////
//
// INTERFACE_MCMF
//
////////////////////////////////////////////////////////////////////////////////

void INTERFACE_MCMF::get_branch_exhanges_to_be_realized(
    vector<int> &sketch_exchange_losses_ids, MCMF_ALGORITHM engine_type)
{
    // this simply inspects the network flow solution
    // and returns inside sketch_exchange_losses all branch exchanges that
    // we will commit at the end of the current interation;

    sketch_exchange_losses_ids.clear();

    if ( _is_empty) return; // no mcmf problem;
    
    // go thru all arcs of type "donor feeder node" -> "switch node" and
    // check if they ar carrying flow as part of the max flow solution;
    // if so, then the close-tie-switch open-sec-switch will have to be 
    // realized;
    int src_vertex_id, des_vertex_id;
    if ( engine_type == MCMF_ALGO_EDMONDS) {
        for ( int j = 0; j < _num_edges; j++) {
            src_vertex_id = _edges_src_vertex[ j];
            des_vertex_id = _edges_des_vertex[ j];

            if ( _vertex_type[ des_vertex_id] == SWITCH_VERTEX) {
                int this_flow = 
                    _engine_mcmf_edmonds.get_fnet( src_vertex_id, des_vertex_id);       
                if ( this_flow > 0) {
                    // found arc in flow graph that contributes to max flow;
                    // record the iindex of the tie-sec switches inside the
                    // sketch array; these branch exchanges will be realized;
                    sketch_exchange_losses_ids.push_back( _vertex_brex_index[ des_vertex_id]);
                    //printf("\n Flow: %d -> %d", src_vertex_id, des_vertex_id);
                }
            }
        }
    }
    else if ( engine_type == MCMF_ALGO_CS2) {
        for ( int j = 0; j < _num_edges; j++) {
            src_vertex_id = _edges_src_vertex[ j];
            des_vertex_id = _edges_des_vertex[ j];

            if ( _vertex_type[ des_vertex_id] == SWITCH_VERTEX) {
                int this_flow = 
                    _engine_mcmf_cs2.get_flow_of_arc( src_vertex_id, des_vertex_id);        
                if ( this_flow > 0) {
                    // found arc in flow graph that contributes to max flow;
                    // record the iindex of the tie-sec switches inside the
                    // sketch array; these branch exchanges will be realized;
                    sketch_exchange_losses_ids.push_back( _vertex_brex_index[ des_vertex_id]);
                    //printf("\n Flow: %d -> %d", src_vertex_id, des_vertex_id);
                }
            }
        }
    }
}


////////////////////////////////////////////////////////////////////////////////
//
// MCMF_EDMONDS
//
////////////////////////////////////////////////////////////////////////////////

// Min cost max flow * (Edmonds-Karp relabelling + fast heap Dijkstra)
// Takes a directed graph where each edge has a capacity ('cap') and a 
// cost per unit of flow ('cost') and returns a maximum flow network
// of minimal cost ('fcost') from s to t. USE mcmf3.cpp FOR DENSE GRAPHS!
//
// PARAMETERS:
//      - cap (global): adjacency matrix where cap[u][v] is the capacity
//          of the edge u->v. cap[u][v] is 0 for non-existent edges.
//      - cost (global): a matrix where cost[u][v] is the cost per unit
//          of flow along the edge u->v. If cap[u][v] == 0, cost[u][v] is
//          ignored. ALL COSTS MUST BE NON-NEGATIVE!
//      - n: the number of vertices ([0, n-1] are considered as vertices).
//      - s: source vertex.
//      - t: sink.
// RETURNS:
//      - the flow
//      - the total cost through 'fcost'
//      - fnet contains the flow network. Careful: both fnet[u][v] and
//          fnet[v][u] could be positive. Take the difference.
// COMPLEXITY:
//      - Worst case: O(m*log(m)*flow  <?  n*m*log(m)*fcost)
// FIELD TESTING:
//      - Valladolid 10594: Data Flow
// REFERENCE:
//      Edmonds, J., Karp, R.  "Theoretical Improvements in Algorithmic
//          Efficieincy for Network Flow Problems".
//      This is a slight improvement of Frank Chu's implementation.

bool MCMF_EDMONDS::dijkstra( int n, int s, int t)
{
    CLEAR_MEM_1D( _d, HALF_INFINITY, _NN); // Note: HALF_INFINITY was 0x3F;
    CLEAR_MEM_1D( _par, -1, _NN);
    CLEAR_MEM_1D( _inq, -1, _NN);
    //for ( int i = 0; i < n; i++ ) _d[i] = HALF_INFINITY, _par[i] = -1;
    _d[s] = _qs = 0;
    _inq[_q[_qs++] = s] = 0;
    _par[s] = n;

    while ( _qs) {
        // get the minimum from _q and bubble down
        int u = _q[0];
        _inq[u] = -1;
        _q[0] = _q[--_qs];
        if ( _qs) { _inq[_q[0]] = 0; }

        for ( int i = 0, j = 2*i + 1, t; j < _qs; i = j, j = 2*i + 1) {
            if ( j + 1 < _qs && _d[_q[j + 1]] < _d[_q[j]]) j++;
            if ( _d[_q[j]] >= _d[_q[i]]) break;
            BUBL;
        }

        // relax edge (u,i) or (i,u) for all i;
        for ( int k = 0, v = _adj[u][k]; k < _deg[u]; v = _adj[u][++k]) {
            // try undoing edge v->u      
            if ( _fnet[v][u] && _d[v] > POTENTIAL(u,v) - _cost[v][u] ) 
                _d[v] = POTENTIAL(u,v) - _cost[v][_par[v] = u];

            // try using edge u->v
            if( _fnet[u][v] < _cap[u][v] && _d[v] > POTENTIAL(u,v) + _cost[u][v]) 
                _d[v] = POTENTIAL(u,v) + _cost[_par[v] = u][v];

            if ( _par[v] == u) {
                // bubble up or decrease key
                if ( _inq[v] < 0) { 
                    _inq[_q[_qs] = v] = _qs; _qs++;
                }
                for ( int i = _inq[v], j = ( i - 1 )/2, t; _d[_q[i]] < _d[_q[j]]; 
                      i = j, j = ( i - 1 )/2 ) {
                    BUBL;
                }
            }
        }
    }

    for ( int i = 0; i < n; i++) { 
        if ( _pi[i] < HALF_INFINITY ) _pi[i] += _d[i];
    }

    return ( _par[t] >= 0);
}

int MCMF_EDMONDS::run_edmonds( int n, int s, int t, int &flow_cost) 
{
    // build the adjacency list
    CLEAR_MEM_1D( _deg, 0, _NN);
    for ( int i = 0; i < n; i++ ) {
        for ( int j = 0; j < n; j++ ) {
            if ( _cap[i][j] || _cap[j][i] ) _adj[i][_deg[i]++] = j;
        }
    }

    CLEAR_MEM_2D( _fnet, 0, _NN);
    CLEAR_MEM_1D( _pi, 0, _NN);
    int flow = flow_cost = 0;

    // repeatedly, find a cheapest path from s to t
    while ( dijkstra( n, s, t) ) {
        // get the bottleneck capacity
        int bot = INT_MAX;
        int this_bot = INT_MAX;
        for ( int v = t, u = _par[v]; v != s; u = _par[v = u] ) {
            //bot <?= (_fnet[v][u] ? _fnet[v][u] : ( _cap[u][v] - _fnet[u][v]));
            this_bot = (_fnet[v][u] ? _fnet[v][u] : ( _cap[u][v] - _fnet[u][v]));
            bot = (bot < this_bot) ? bot : this_bot;
        }
        
        // update the flow network
        for ( int v = t, u = _par[v]; v != s; u = _par[v = u] ) {
            if ( _fnet[v][u] ) {
                _fnet[v][u] -= bot; 
                flow_cost -= bot * _cost[v][u]; }
            else { 
                _fnet[u][v] += bot; 
                flow_cost += bot * _cost[u][v]; 
            }
        }

        flow += bot;
    }

    return flow;
}


////////////////////////////////////////////////////////////////////////////////
//
// MCMF_CS2
//
////////////////////////////////////////////////////////////////////////////////

#define WHITE 0
#define GREY  1
#define BLACK 2
#define OPEN( a ) ( a->rez_capacity() > 0 )
#define CLOSED( a ) ( a->rez_capacity() <= 0 )
#define REDUCED_COST( i, j, a ) ( i->price() + a->cost() - j->price() )
#define FEASIBLE( i, j, a ) ( i->price() + a->cost() < j->price() )
#define ADMISSIBLE( i, j, a ) ( OPEN( a ) && FEASIBLE( i, j, a ) )
#define SUSPENDED( i, a ) ( a < i->first() )


#define REMOVE_FROM_EXCESS_Q( i )               \
    {                                           \
        i = _excq_first;                        \
        _excq_first = i -> q_next();            \
        i ->set_q_next( _sentinel_node );       \
    }

#define STACKQ_POP( i )                         \
    {                                           \
        i = _excq_first;                        \
        _excq_first = i -> q_next();            \
        i ->set_q_next( _sentinel_node );       \
    }

#define GET_FROM_BUCKET( i, b )                 \
    {                                           \
        i = ( b -> p_first() );                 \
        b ->set_p_first( i -> b_next() );       \
    }
#define REMOVE_FROM_BUCKET( i, b )                              \
    {                                                           \
        if ( i == ( b -> p_first() ) )                          \
            b ->set_p_first( i -> b_next() );                   \
        else                                                    \
            {                                                   \
                ( i -> b_prev() )->set_b_next( i -> b_next() ); \
                ( i -> b_next() )->set_b_prev( i -> b_prev() ); \
            }                                                   \
    }


void MCMF_CS2::err_end( int cc)
{
    // abnormal finish
    printf ("\nError %d\n", cc );
    // 2 - problem is unfeasible
    // 5 - allocation fault
    // 6 - price overflow
    exit( cc);
}

void MCMF_CS2::allocate_arrays()
{
    // (1) allocate memory for 'nodes', 'arcs' and internal arrays;

    _nodes    = (NODE*) calloc ( _n+2,   sizeof(NODE) );
    _arcs     = (ARC*)  calloc ( 2*_m+1, sizeof(ARC) );
    _cap      = (long*) calloc ( 2*_m,   sizeof(long) );

    _arc_tail = (long*) calloc ( 2*_m,   sizeof(long) );
    _arc_first= (long*) calloc ( _n+2,   sizeof(long) );
    // arc_first [ 0 .. n+1 ] = 0 - initialized by calloc;

    for ( NODE *in = _nodes; in <= _nodes + _n; in ++ ) {
        in->set_excess( 0);
    }
    if ( _nodes == NULL || _arcs == NULL || _arc_first == NULL || _arc_tail == NULL) {
        printf("Error:  Memory allocation problem inside CS2\n");
        exit( 1);
    }

    // (2) resets;
    _pos_current = 0;
    _arc_current = _arcs; // set "current" pointer to the first arc
    _node_max = 0;
    _node_min = _n;
    _max_cost = 0;
    _total_p = _total_n = 0;
    // at this moment we are ready to add arcs and build the network,
    // by using set_arc()...
}

void MCMF_CS2::set_arc( long tail_node_id, long head_node_id,
    long low_bound, long up_bound, // up_bound is basically capacity;
    price_t cost)
{
    // DIMACS format:
    // c arc has <tail> <head> <capacity l.b.> <capacity u.b> <cost>

    if ( tail_node_id < 0 || tail_node_id > _n || 
        head_node_id < 0 || head_node_id > _n ) {
        printf("Error:  Arc with head or tail out of bounds inside CS2\n");
        exit( 1);
    }
    if ( up_bound < 0 ) {
        up_bound = MAX_32;
        printf("Warning:  Infinite capacity replaced by BIGGEST_FLOW\n");
    }
    if ( low_bound < 0 || low_bound > up_bound ) {
        printf("Error:  Wrong capacity bounds inside CS2\n");
        exit( 1);
    }

    // no of arcs incident to node i is placed in _arc_first[i+1]
    _arc_first[tail_node_id + 1] ++; 
    _arc_first[head_node_id + 1] ++;
    _i_node = _nodes + tail_node_id;
    _j_node = _nodes + head_node_id;

    // store information about the arc
    _arc_tail[_pos_current]   = tail_node_id;
    _arc_tail[_pos_current+1] = head_node_id;

    _arc_current->set_head( _j_node );
    _arc_current->set_rez_capacity( up_bound - low_bound );
    _arc_current->set_rez_capacity_default( up_bound - low_bound ); // default;
    _cap[_pos_current] = up_bound;
    _arc_current->set_cost( cost );
    _arc_current->set_cost_default( cost ); // default;
    _arc_current->set_sister( _arc_current + 1 );


    ( _arc_current + 1 )->set_head( _nodes + tail_node_id );
    ( _arc_current + 1 )->set_rez_capacity( 0 );
    ( _arc_current + 1 )->set_rez_capacity_default( 0 ); // default;
    _cap[_pos_current+1] = 0;
    ( _arc_current + 1 )->set_cost( -cost );
    ( _arc_current + 1 )->set_cost_default( -cost ); // default;
    ( _arc_current + 1 )->set_sister( _arc_current );

    _i_node->dec_excess( low_bound );
    _j_node->inc_excess( low_bound );

    // searching for minimum and maximum node
    if ( head_node_id < _node_min ) _node_min = head_node_id;
    if ( tail_node_id < _node_min ) _node_min = tail_node_id;
    if ( head_node_id > _node_max ) _node_max = head_node_id;
    if ( tail_node_id > _node_max ) _node_max = tail_node_id;

    if ( cost < 0 ) cost = -cost;
    if ( cost > _max_cost && up_bound > 0 ) _max_cost = cost;

    // prepare for next arc to be added;
    _arc_current += 2;
    _pos_current += 2;
}

void MCMF_CS2::set_supply_demand_of_node( long id, long excess)
{
    // set supply and demand of nodes; not used for transhipment nodes;
    if ( id < 0 || id > _n ) {
        printf("Error:  Unbalanced problem inside CS2\n");
        exit( 1);
    }
    (_nodes + id)->set_excess( excess);
    if ( excess > 0) _total_p += excess;
    if ( excess < 0) _total_n -= excess;
}

void MCMF_CS2::pre_processing()
{
    // called after the arcs were just added and before run_cs2();
    // ordering arcs - linear time algorithm;
    long i;
    long last, arc_num, arc_new_num;;
    long tail_node_id;
    NODE *head_p;
    ARC *arc_new, *arc_tmp; 
    long up_bound;
    price_t cost;
    long up_bound_default;
    price_t cost_default;
    excess_t cap_out; // sum of outgoing capacities
    excess_t cap_in; // sum of incoming capacities

    if ( ABS( _total_p - _total_n ) > 0.5 ) {
        printf("Error:  Unbalanced problem inside CS2\n");
        exit( 1);
    }
    
    // first arc from the first node
    ( _nodes + _node_min )->set_first( _arcs );

    // before below loop arc_first[i+1] is the number of arcs outgoing from i;
    // after this loop arc_first[i] is the position of the first 
    // outgoing from node i arcs after they would be ordered;
    // this value is transformed to pointer and written to node.first[i]
    for ( i = _node_min + 1; i <= _node_max + 1; i ++ ) {
        _arc_first[i] += _arc_first[i-1];
        ( _nodes + i )->set_first( _arcs + _arc_first[i] );
    }

    // scanning all the nodes except the last
    for ( i = _node_min; i < _node_max; i ++ ) {

        last = ( ( _nodes + i + 1 )->first() ) - _arcs;
        // arcs outgoing from i must be cited    
        // from position arc_first[i] to the position
        // equal to initial value of arc_first[i+1]-1

        for ( arc_num = _arc_first[i]; arc_num < last; arc_num ++ ) {
            tail_node_id = _arc_tail[arc_num];

            while ( tail_node_id != i ) {
                // the arc no  arc_num  is not in place because arc cited here
                // must go out from i;
                // we'll put it to its place and continue this process
                // until an arc in this position would go out from i

                arc_new_num = _arc_first[tail_node_id];
                _arc_current = _arcs + arc_num;
                arc_new = _arcs + arc_new_num;
        
                // arc_current must be cited in the position arc_new    
                // swapping these arcs:

                head_p = arc_new->head();
                arc_new->set_head( _arc_current->head() );
                _arc_current->set_head( head_p );

                up_bound          = _cap[arc_new_num];
                _cap[arc_new_num] = _cap[arc_num];
                _cap[arc_num]     = up_bound;

                up_bound = arc_new->rez_capacity();
                arc_new->set_rez_capacity( _arc_current->rez_capacity() );
                _arc_current->set_rez_capacity( up_bound);
                cost = arc_new->cost();
                arc_new->set_cost( _arc_current->cost() );
                _arc_current->set_cost( cost );
                // defaults too;
                up_bound_default = arc_new->rez_capacity_default();
                arc_new->set_rez_capacity_default( _arc_current->rez_capacity_default() );
                _arc_current->set_rez_capacity_default( up_bound_default);
                cost_default = arc_new->cost_default();
                arc_new->set_cost_default( _arc_current->cost_default() );
                _arc_current->set_cost_default( cost_default );

                if ( arc_new != _arc_current->sister() ) {
                    arc_tmp = arc_new->sister();
                    arc_new->set_sister( _arc_current->sister() );
                    _arc_current->set_sister( arc_tmp );

                    _arc_current->sister()->set_sister( _arc_current );
                    arc_new->sister()->set_sister( arc_new );
                }

                _arc_tail[arc_num] = _arc_tail[arc_new_num];
                _arc_tail[arc_new_num] = tail_node_id;

                // we increase arc_first[tail_node_id]
                _arc_first[tail_node_id] ++ ;

                tail_node_id = _arc_tail[arc_num];
            }
        }
        // all arcs outgoing from  i  are in place
    }       
    // arcs are ordered by now!


    // testing network for possible excess overflow
    for ( NODE *ndp = _nodes + _node_min; ndp <= _nodes + _node_max; ndp ++ ) {
        cap_in  =   ( ndp->excess() );
        cap_out = - ( ndp->excess() );
        for ( _arc_current = ndp->first(); _arc_current != (ndp+1)->first(); 
              _arc_current ++ ) {
            arc_num = _arc_current - _arcs;
            if ( _cap[arc_num] > 0 ) cap_out += _cap[arc_num];
            if ( _cap[arc_num] == 0 ) 
                cap_in += _cap[ _arc_current->sister() - _arcs ];
        }
    }
    if ( _node_min < 0 || _node_min > 1 ) {
        printf("Error:  Node ids must start from 0 or 1 inside CS2\n");
        exit( 1);
    }

    // adjustments due to nodes' ids being between _node_min - _node_max;
    _n = _node_max - _node_min + 1;
    _nodes = _nodes + _node_min;

}

void MCMF_CS2::cs2_initialize()
{
    // initialization; 
    // called after allocate_arrays() and all nodes and arcs have been inputed;

    NODE *i; // current node
    ARC *a; // current arc
    ARC *a_stop;
    BUCKET *b; // current bucket
    long df;

    _f_scale = (long) SCALE_DEFAULT;
    _sentinel_node = _nodes + _n;
    _sentinel_arc  = _arcs + _m;

    for ( i = _nodes; i != _sentinel_node; i ++ ) {
        i->set_price( 0);
        i->set_suspended( i->first());
        i->set_q_next( _sentinel_node);
    }

    _sentinel_node->set_first( _sentinel_arc);
    _sentinel_node->set_suspended( _sentinel_arc);

    // saturate negative arcs, e.g. in the circulation problem case
    for ( i = _nodes; i != _sentinel_node; i ++ ) {
        for ( a = i->first(), a_stop = (i + 1)->suspended(); a != a_stop; a ++ ) {
            if ( a->cost() < 0) {
                if ( ( df = a->rez_capacity()) > 0) {
                    increase_flow( i, a->head(), a, df);
                }
            }
        }
    }

    _dn = _n + 1;
    if ( _no_zero_cycles == true) { // NO_ZERO_CYCLES
        _dn = 2 * _dn;
    }

    for ( a = _arcs; a != _sentinel_arc; a ++ ) {
        a->multiply_cost( _dn);
    }

    if ( _no_zero_cycles == true) { // NO_ZERO_CYCLES
        for ( a = _arcs; a != _sentinel_arc; a ++ ) {
            if ((a->cost() == 0) && (a->sister()->cost() == 0)) {
                a->set_cost( 1);
                a->sister()->set_cost( -1);
            }
        }
    }

    if ((double) _max_cost * (double) _dn > MAX_64) {
        printf("Warning:  Arc lengths too large, overflow possible\n");
    }
    _mmc = _max_cost * _dn;

    _linf = (long) (_dn * ceil(_f_scale) + 2);


    _buckets = (BUCKET*) calloc ( _linf, sizeof(BUCKET));
    if ( _buckets == NULL )
        err_end( ALLOCATION_FAULT);

    _l_bucket = _buckets + _linf;

    _dnode = new NODE; // used as reference;

    for ( b = _buckets; b != _l_bucket; b ++ ) {
        reset_bucket( b);
    }

    _epsilon = _mmc;
    if ( _epsilon < 1) {
        _epsilon = 1;
    }

    _price_min = -PRICE_MAX;

    _cut_off_factor = CUT_OFF_COEF * pow( (double)_n, CUT_OFF_POWER);

    _cut_off_factor = MAX( _cut_off_factor, CUT_OFF_MIN);

    _n_ref = 0;

    _flag_price = 0;

    _dummy_node = &_d_node;

    _excq_first = NULL;

    // debug;
    //print_graph();
    //print_stuff();
}

void MCMF_CS2::cs2_initialize_light()
{
    NODE *i; // current node
    ARC *a; // current arc
    ARC *a_stop;
    BUCKET *b; // current bucket
    long df;

    _f_scale = (long) SCALE_DEFAULT;
    _sentinel_node = _nodes + _n;
    _sentinel_arc  = _arcs + _m;

    for ( i = _nodes; i != _sentinel_node; i ++ ) {
        i->set_price( 0);
        i->set_suspended( i->first());
        i->set_q_next( _sentinel_node);
    }

    _sentinel_node->set_first( _sentinel_arc);
    _sentinel_node->set_suspended( _sentinel_arc);

    _dn = _n + 1;

    for ( a = _arcs; a != _sentinel_arc; a ++ ) {
        a->multiply_cost( _dn);
    }

    _mmc = _max_cost * _dn;

    _linf = (long) (_dn * ceil(_f_scale) + 2);

    _l_bucket = _buckets + _linf;

    for ( b = _buckets; b != _l_bucket; b ++ ) {
        reset_bucket( b);
    }

    _epsilon = _mmc;
    if ( _epsilon < 1) {
        _epsilon = 1;
    }

    _price_min = -PRICE_MAX;

    _cut_off_factor = CUT_OFF_COEF * pow( (double)_n, CUT_OFF_POWER);

    _cut_off_factor = MAX( _cut_off_factor, CUT_OFF_MIN);

    _n_ref = 0;

    _flag_price = 0;

    _dummy_node = &_d_node;

    _excq_first = NULL;

    // debug;
    //print_graph();
    //print_stuff();
}

void MCMF_CS2::up_node_scan( NODE *i)
{
    NODE *j; // opposite node
    ARC *a; // (i, j)
    ARC *a_stop; // first arc from the next node
    ARC *ra; // (j, i)
    BUCKET *b_old; // old bucket contained j
    BUCKET *b_new; // new bucket for j
    long i_rank;
    long j_rank; // ranks of nodes
    long j_new_rank;             
    price_t rc; // reduced cost of (j, i)
    price_t dr; // rank difference

    _n_scan ++;

    i_rank = i->rank();

    // scanning arcs;
    for ( a = i->first(), a_stop = (i + 1)->suspended(); a != a_stop; a ++ ) {

        ra = a->sister();

        if ( OPEN ( ra ) ) {
            j = a->head();
            j_rank = j->rank();

            if ( j_rank > i_rank ) {
                if ( ( rc = REDUCED_COST( j, i, ra ) ) < 0 ) {
                    j_new_rank = i_rank;
                } else {
                    dr = rc / _epsilon;
                    j_new_rank = ( dr < _linf ) ? i_rank + (long)dr + 1 : _linf;
                }

                if ( j_rank > j_new_rank ) {
                    j->set_rank( j_new_rank);
                    j->set_current( ra);

                    if ( j_rank < _linf ) {
                        b_old = _buckets + j_rank;
                        REMOVE_FROM_BUCKET( j, b_old );
                    }

                    b_new = _buckets + j_new_rank;
                    insert_to_bucket( j, b_new );
                }
            }
        }
    }

    i->dec_price( i_rank * _epsilon);
    i->set_rank( -1);
}

void MCMF_CS2::price_update()
{
    register NODE *i;
    excess_t remain;
    // total excess of unscanned nodes with positive excess;
    BUCKET *b; // current bucket;
    price_t dp; // amount to be subtracted from prices;

    _n_update ++;

    for ( i = _nodes; i != _sentinel_node; i ++ ) {
        if ( i->excess() < 0 ) {
            insert_to_bucket( i, _buckets );
            i->set_rank( 0);
        } else {
            i->set_rank( _linf);
        }
    }

    remain = _total_excess;
    if ( remain < 0.5 ) return;

    // scanning buckets, main loop;
    for ( b = _buckets; b != _l_bucket; b ++ ) {

        while ( nonempty_bucket( b) ) {

            GET_FROM_BUCKET( i, b );
            up_node_scan( i );

            if ( i ->excess() > 0 ) {
                remain -= ( i->excess());
                if ( remain <= 0 ) break; 
            }
        }
        if ( remain <= 0 ) break; 
    }

    if ( remain > 0.5 ) _flag_updt = 1;

    // finishup 
    // changing prices for nodes which were not scanned during main loop;
    dp = ( b - _buckets ) * _epsilon;

    for ( i = _nodes; i != _sentinel_node; i ++ ) {

        if ( i->rank() >= 0 ) {
            if ( i->rank() < _linf ) {
                REMOVE_FROM_BUCKET( i, ( _buckets + i->rank()) );
            }
            if ( i->price() > _price_min ) {
                i->dec_price( dp);
            }           
        }
    }
}

int MCMF_CS2::relabel( NODE *i)
{
    register ARC *a; // current arc from i
    register ARC *a_stop; // first arc from the next node
    register ARC *a_max; // arc which provides maximum price
    register price_t p_max; // current maximal price
    register price_t i_price; // price of node  i
    register price_t dp; // current arc partial residual cost

    p_max = _price_min;
    i_price = i->price();

    a_max = NULL;

    // 1/2 arcs are scanned;
    for ( a = i->current() + 1, a_stop = (i + 1)->suspended(); a != a_stop; a ++ ) {
        
        if ( OPEN(a) && ( (dp = (a->head()->price() - a->cost())) > p_max ) ) {
            if ( i_price < dp ) {
                i->set_current( a);
                return ( 1);
            }
            p_max = dp;
            a_max = a;
        }
    }

    // 2/2 arcs are scanned;
    for ( a = i->first(), a_stop = i->current() + 1; a != a_stop; a ++ ) {
        if ( OPEN( a) && ( (dp = (a->head()->price() - a->cost())) > p_max ) ) {
            if ( i_price < dp ) {
                i->set_current( a);
                return ( 1);
            }
            p_max = dp;
            a_max = a;
        }
    }

    // finishup
    if ( p_max != _price_min ) {
        i->set_price( p_max - _epsilon);
        i->set_current( a_max);
    } 
    else { // node can't be relabelled;
        if ( i->suspended() == i->first() ) {
            if ( i->excess() == 0 ) {
                i->set_price( _price_min);
            } else {
                if ( _n_ref == 1 ) {
                    printf("\nError: in relabel()");
                    err_end( UNFEASIBLE );
                } else {
                    err_end( PRICE_OFL );
                }
            }
        } else { // node can't be relabelled because of suspended arcs;
            _flag_price = 1;
        }
    }

    _n_relabel ++;
    _n_rel ++;
    return ( 0);
}

void MCMF_CS2::discharge( NODE *i)
{
    register ARC *a;// an arc from i
    register NODE *j; // head of a
    register long df; // amoumt of flow to be pushed through a
    excess_t j_exc; // former excess of j

    _n_discharge ++;
    //printf("\n-- Discharge node rank: %d  n_discharge: %d", i->rank(), _n_discharge);

    a = i->current();
    j = a->head();

    if ( !ADMISSIBLE( i, j, a ) ) { 
        relabel( i );
        a = i->current();
        j = a->head();
    }

    while ( 1 ) {

        j_exc = j->excess();
        if ( j_exc >= 0 ) {

            df = MIN( i->excess(), a->rez_capacity() );
            if ( j_exc == 0) _n_src++;
            increase_flow( i, j, a, df ); // INCREASE_FLOW 
            _n_push ++;

            if ( out_of_excess_q( j ) ) {
                insert_to_excess_q( j );
            }
        } 
        else { // j_exc < 0;

            df = MIN( i->excess(), a->rez_capacity() );
            increase_flow( i, j, a, df ); // INCREASE_FLOW 
            _n_push ++;

            if ( j->excess() >= 0 ) {
                if ( j->excess() > 0 ) {
                    _n_src ++;
                    relabel( j );
                    insert_to_excess_q( j );
                }
                _total_excess += j_exc;
            }
            else {
                _total_excess -= df;
            }
        }
  
        if ( i->excess() <= 0) _n_src --;
        if ( i->excess() <= 0 || _flag_price ) break;

        relabel( i );

        a = i->current();
        j = a->head();
    }

    i->set_current( a);
}

int MCMF_CS2::price_in()
{
    NODE *i; // current node
    NODE *j;
    ARC *a; // current arc from i
    ARC *a_stop; // first arc from the next node
    ARC *b; // arc to be exchanged with suspended
    ARC *ra; // opposite to a
    ARC *rb; // opposite to b
    price_t rc; // reduced cost
    int n_in_bad; // number of priced_in arcs with negative reduced cost
    int bad_found; // if 1 we are at the second scan if 0 we are at the first scan
    excess_t i_exc; // excess of i
    excess_t df; // an amount to increase flow


    bad_found = 0;
    n_in_bad = 0;

 restart:

    for ( i = _nodes; i != _sentinel_node; i ++ ) {

        for ( a = i->first() - 1, a_stop = i->suspended() - 1; a != a_stop; a -- ) {

            rc = REDUCED_COST( i, a->head(), a );
            if ( ( rc < 0) && ( a->rez_capacity() > 0) ) { // bad case;
                if ( bad_found == 0 ) {
                    bad_found = 1;
                    update_cut_off();
                    goto restart;
                }
                df = a->rez_capacity();
                increase_flow( i, a->head(), a, df );
        
                ra = a->sister();
                j  = a->head();
        
                i->dec_first();
                b = i->first();
                exchange( a, b );
        
                if ( SUSPENDED( j, ra ) ) {
                    j->dec_first();
                    rb = j->first();
                    exchange( ra, rb );
                }
        
                n_in_bad ++; 
            }
            else {
                if ( ( rc < _cut_on ) && ( rc > -_cut_on ) ) {
                    i->dec_first();
                    b = i->first();
                    exchange( a, b );
                }
            }
        }
    }


    if ( n_in_bad != 0 ) {

        _n_bad_pricein ++;

        // recalculating excess queue;
        _total_excess = 0;
        _n_src = 0;
        reset_excess_q();

        for ( i = _nodes; i != _sentinel_node; i ++ ) {
            i->set_current( i->first());
            i_exc = i->excess();
            if ( i_exc > 0 ) { // i is a source;
                _total_excess += i_exc;
                _n_src ++;
                insert_to_excess_q( i );
            }
        }

        insert_to_excess_q( _dummy_node );
    }

    if ( _time_for_price_in == TIME_FOR_PRICE_IN2)
        _time_for_price_in = TIME_FOR_PRICE_IN3;
    if ( _time_for_price_in == TIME_FOR_PRICE_IN1)
        _time_for_price_in = TIME_FOR_PRICE_IN2;

    return ( n_in_bad);
}

int MCMF_CS2::refine()
{
    NODE *i; // current node
    excess_t i_exc; // excess of i
    long np, nr, ns; // variables for additional print
    int pr_in_int; // current number of updates between price_in

    np = _n_push; 
    nr = _n_relabel; 
    ns = _n_scan;

    _n_refine ++;
    _n_ref ++;
    _n_rel = 0;
    pr_in_int = 0;

    // initialize;
    _total_excess = 0;
    _n_src = 0;
    reset_excess_q();

    _time_for_price_in = TIME_FOR_PRICE_IN1;

    for ( i = _nodes; i != _sentinel_node; i ++ ) {
        i->set_current( i->first());
        i_exc = i->excess();
        if ( i_exc > 0 ) { // i is a source 
            _total_excess += i_exc;
            _n_src++;
            insert_to_excess_q( i );
        }
    }

    if ( _total_excess <= 0 ) return 0;

    // (2) main loop

    while ( 1 ) {

        if ( empty_excess_q() ) {
            if ( _n_ref > PRICE_OUT_START ) {
                pr_in_int = 0;
                price_in();
            }
      
            if ( empty_excess_q() ) break;
        }
        
        REMOVE_FROM_EXCESS_Q( i );

        // push all excess out of i
        if ( i->excess() > 0 ) {
            discharge( i );

            if ( time_for_update() || _flag_price ) {
                if ( i->excess() > 0 ) {
                    insert_to_excess_q( i );
                }

                if ( _flag_price && ( _n_ref > PRICE_OUT_START ) ) {
                    pr_in_int = 0;
                    price_in();
                    _flag_price = 0;
                }

                price_update();

                while ( _flag_updt ) {
                    if ( _n_ref == 1 ) {
                        return 2;
                        printf("\nError: in refine()");
                        err_end( UNFEASIBLE );
                    } else {
                        _flag_updt = 0;
                        update_cut_off();
                        _n_bad_relabel ++;
                        pr_in_int = 0;
                        price_in();
                        price_update(); // will set _flag_updt = 1;
                    }
                }
                _n_rel = 0;

                if ( _n_ref > PRICE_OUT_START && (pr_in_int ++ > _time_for_price_in) ) {
                    pr_in_int = 0;
                    price_in();
                }
            }
        }
    }

    return 0;
}

int MCMF_CS2::price_refine()
{
    NODE *i; // current node
    NODE *j; // opposite node
    NODE *ir; // nodes for passing over the negative cycle
    NODE *is;
    ARC *a; // arc (i,j)
    ARC *a_stop; // first arc from the next node
    ARC *ar;
    long bmax;            // number of farest nonempty bucket
    long i_rank;          // rank of node i
    long j_rank;         // rank of node j
    long j_new_rank;      // new rank of node j
    BUCKET *b;              // current bucket
    BUCKET *b_old;          // old and new buckets of current node
    BUCKET *b_new;
    price_t rc = 0; // reduced cost of a
    price_t dr; // ranks difference
    price_t dp;
    int cc;              
    // return code: 1 - flow is epsilon optimal
    // 0 - refine is needed       
    long df; // cycle capacity
    int nnc; // number of negative cycles cancelled during one iteration
    int snc; // total number of negative cycle cancelled

    _n_prefine ++;

    cc = 1;
    snc = 0;

    _snc_max = ( _n_ref >= START_CYCLE_CANCEL) ? MAX_CYCLES_CANCELLED : 0;


    // (1) main loop
    // while negative cycle is found or eps-optimal solution is constructed
    while ( 1 ) { 

        nnc = 0;
        for ( i = _nodes; i != _sentinel_node; i ++ ) {
            i->set_rank( 0);
            i->set_inp( WHITE);
            i->set_current( i->first());
        }
        reset_stackq();

        for ( i = _nodes; i != _sentinel_node; i ++ ) {
            if ( i->inp() == BLACK ) continue;

            i->set_b_next( NULL);

            // deapth first search 
            while ( 1 ) {
                i->set_inp( GREY);

                // scanning arcs from node i starting from current
                for ( a = i->current(), a_stop = (i + 1)->suspended(); a != a_stop; a ++) {
                    if ( OPEN( a ) ) {
                        j = a->head();
                        if ( REDUCED_COST ( i, j, a ) < 0 ) {
                            if ( j->inp() == WHITE ) { // fresh node  - step forward 
                                i->set_current( a);
                                j->set_b_next( i);
                                i = j;
                                a = j->current();
                                a_stop = (j+1)->suspended();
                                break;
                            }

                            if ( j->inp() == GREY ) { // cycle detected 
                                cc = 0;
                                nnc ++;
                                i->set_current( a);
                                is = ir = i;
                                df = MAX_32;

                                while ( 1 ) {
                                    ar = ir->current();
                                    if ( ar->rez_capacity() <= df ) {
                                        df = ar->rez_capacity();
                                        is = ir;
                                    }
                                    if ( ir == j ) break;
                                    ir = ir->b_next();
                                }

                                ir = i;

                                while ( 1 ) {
                                    ar = ir->current();
                                    increase_flow( ir, ar->head(), ar, df);
                                    if ( ir == j ) break;
                                    ir = ir->b_next();
                                }

                                if ( is != i ) {
                                    for ( ir = i; ir != is; ir = ir->b_next() ) {
                                        ir->set_inp( WHITE);
                                    }
                                    i = is;
                                    a = is->current() + 1;
                                    a_stop = (is+1)->suspended();
                                    break;
                                }
                            }
                        }
                        // if j-color is BLACK - continue search from i
                    }
                } // all arcs from i are scanned 

                if ( a == a_stop ) {
                    // step back 
                    i->set_inp( BLACK);
                    _n_prscan1 ++;
                    j = i->b_next();
                    stackq_push( i );
                    if ( j == NULL ) break;
                    i = j;
                    i->inc_current();
                }

            } // end of deapth first search
        } // all nodes are scanned


        // () no negative cycle
        // computing longest paths with eps-precision

        snc += nnc;
        if ( snc < _snc_max ) cc = 1;
        if ( cc == 0 ) break;
        bmax = 0;

        while ( nonempty_stackq() ) {

            _n_prscan2 ++;
            STACKQ_POP( i );
            i_rank = i->rank();
            for ( a = i->first(), a_stop = (i + 1)->suspended(); a != a_stop; a ++) {

                if ( OPEN( a ) ) {
                    j  = a->head();
                    rc = REDUCED_COST( i, j, a );

                    if ( rc < 0 ) { // admissible arc;
                        dr = (price_t) (( - rc - 0.5 ) / _epsilon);
                        if (( j_rank = dr + i_rank ) < _linf ) {
                            if ( j_rank > j->rank() )
                                j->set_rank( j_rank);
                        }
                    }
                }
            } // all arcs from i are scanned 

            if ( i_rank > 0 ) {
                if ( i_rank > bmax ) bmax = i_rank;
                b = _buckets + i_rank;
                insert_to_bucket( i, b );
            }
        } // end of while-cycle: all nodes are scanned - longest distancess are computed;


        if ( bmax == 0 ) // preflow is eps-optimal;
            { break; }


        for ( b = _buckets + bmax; b != _buckets; b -- ) {
            i_rank = b - _buckets;
            dp = i_rank * _epsilon;

            while ( nonempty_bucket( b) ) {
                GET_FROM_BUCKET( i, b );
                _n_prscan ++;

                for ( a = i->first(), a_stop = (i + 1)->suspended(); a != a_stop; a ++) {
                    if ( OPEN( a ) ) {
                        j = a->head();
                        j_rank = j->rank();
                        if ( j_rank < i_rank ) {
                            rc = REDUCED_COST( i, j, a );
                            if ( rc < 0 ) {
                                j_new_rank = i_rank;
                            } else {
                                dr = rc / _epsilon;
                                j_new_rank = ( dr < _linf ) ? i_rank - ( (long)dr + 1 ) : 0;
                            }
                            if ( j_rank < j_new_rank ) {
                                if ( cc == 1 ) {
                                    j->set_rank( j_new_rank);
                                    if ( j_rank > 0 ) {
                                        b_old = _buckets + j_rank;
                                        REMOVE_FROM_BUCKET( j, b_old );
                                    }
                                    b_new = _buckets + j_new_rank;
                                    insert_to_bucket( j, b_new );
                                }
                                else {
                                    df = a->rez_capacity();
                                    increase_flow( i, j, a, df );
                                }
                            }
                        }
                    } // end if opened arc 
                } // all arcs are scanned

                i->dec_price( dp);

            } // end of while-cycle: the bucket is scanned 
        } // end of for-cycle: all buckets are scanned 

        if ( cc == 0 ) break;

    } // end of main loop



    // (2) finish
    // if refine needed - saturate non-epsilon-optimal arcs;

    if ( cc == 0 ) { 
        for ( i = _nodes; i != _sentinel_node; i ++) {
            for ( a = i->first(), a_stop = (i + 1)->suspended(); a != a_stop; a ++) {
                if ( REDUCED_COST( i, a->head(), a ) < - _epsilon ) {
                    if ( ( df = a->rez_capacity() ) > 0 ) {
                        increase_flow( i, a->head(), a, df );
                    }
                }
            }
        }
    }

    return ( cc );
}

void MCMF_CS2::compute_prices()
{
    NODE *i;          // current node
    NODE *j;             // opposite node
    ARC *a;             // arc (i,j)
    ARC *a_stop;        // first arc from the next node
    long bmax;           // number of farest nonempty bucket
    long i_rank;         // rank of node i
    long j_rank;         // rank of node j
    long j_new_rank;     // new rank of node j
    BUCKET *b; // current bucket
    BUCKET *b_old; // old and new buckets of current node
    BUCKET *b_new;
    price_t rc; // reduced cost of a
    price_t dr; // ranks difference
    price_t dp;
    int cc; // return code: 1 - flow is epsilon optimal 0 - refine is needed

    _n_prefine ++;
    cc = 1;

    // (1) main loop
    // while negative cycle is found or eps-optimal solution is constructed
    while ( 1 ) {

        for ( i = _nodes; i != _sentinel_node; i ++) {
            i->set_rank( 0);
            i->set_inp( WHITE);
            i->set_current( i->first());
        }
        reset_stackq();

        for ( i = _nodes; i != _sentinel_node; i ++ ) {
            if ( i->inp() == BLACK ) continue;

            i->set_b_next( NULL);
            // depth first search
            while ( 1 ) {
                i->set_inp( GREY);

                // scanning arcs from node i
                for ( a = i->suspended(), a_stop = (i + 1)->suspended(); a != a_stop; a ++) {
                    if ( OPEN( a ) ) {
                        j = a->head();
                        if ( REDUCED_COST( i, j, a ) < 0 ) {
                            if ( j->inp() == WHITE ) { // fresh node  - step forward
                                i->set_current( a);
                                j->set_b_next( i);
                                i = j;
                                a = j->current();
                                a_stop = (j+1)->suspended();
                                break;
                            }

                            if ( j->inp() == GREY ) { // cycle detected; should not happen
                                cc = 0;
                            }                     
                        }
                        // if j-color is BLACK - continue search from i
                    } 
                } // all arcs from i are scanned 

                if ( a == a_stop ) {
                    // step back
                    i->set_inp( BLACK);
                    _n_prscan1 ++;
                    j = i->b_next();
                    stackq_push( i );
                    if ( j == NULL ) break;
                    i = j;
                    i->inc_current();
                }

            } // end of deapth first search
        } // all nodes are scanned


        // no negative cycle
        // computing longest paths

        if ( cc == 0 ) break;
        bmax = 0;

        while ( nonempty_stackq() ) {
            _n_prscan2 ++;
            STACKQ_POP( i );
            i_rank = i->rank();
            for ( a = i->suspended(), a_stop = (i + 1)->suspended(); a != a_stop; a ++) {
                if ( OPEN( a ) ) {
                    j  = a->head();
                    rc = REDUCED_COST( i, j, a );


                    if ( rc < 0 ) {// admissible arc
                        dr = - rc;
                        if (( j_rank = dr + i_rank ) < _linf ) {
                            if ( j_rank > j->rank() )
                                j->set_rank( j_rank);
                        }
                    }
                }
            } // all arcs from i are scanned

            if ( i_rank > 0 ) {
                if ( i_rank > bmax ) bmax = i_rank;
                b = _buckets + i_rank;
                insert_to_bucket( i, b );
            }
        } // end of while-cycle: all nodes are scanned - longest distancess are computed;

        if ( bmax == 0 )
            { break; }

        for ( b = _buckets + bmax; b != _buckets; b -- ) {
            i_rank = b - _buckets;
            dp = i_rank;

            while ( nonempty_bucket( b) ) {
                GET_FROM_BUCKET( i, b );
                _n_prscan ++;

                for ( a = i->suspended(), a_stop = (i + 1)->suspended(); a != a_stop; a ++) {
                    if ( OPEN( a ) ) {
                        j = a->head();
                        j_rank = j->rank();
                        if ( j_rank < i_rank ) {
                            rc = REDUCED_COST( i, j, a );
 
                            if ( rc < 0 ) {
                                j_new_rank = i_rank;
                            } else {
                                dr = rc;
                                j_new_rank = ( dr < _linf ) ? i_rank - ( (long)dr + 1 ) : 0;
                            }
                            if ( j_rank < j_new_rank ) {
                                if ( cc == 1 ) {
                                    j->set_rank( j_new_rank);
                                    if ( j_rank > 0 ) {
                                        b_old = _buckets + j_rank;
                                        REMOVE_FROM_BUCKET( j, b_old );
                                    }
                                    b_new = _buckets + j_new_rank;
                                    insert_to_bucket( j, b_new );
                                }
                            }
                        }
                    } // end if opened arc
                } // all arcs are scanned

                i->dec_price( dp);

            } // end of while-cycle: the bucket is scanned
        } // end of for-cycle: all buckets are scanned

        if ( cc == 0 ) break;

    } // end of main loop
} 

void MCMF_CS2::price_out()
{
    NODE *i; // current node
    ARC *a; // current arc from i 
    ARC *a_stop; // first arc from the next node 
    ARC *b; // arc to be exchanged with suspended
    double n_cut_off; // -cut_off
    double rc; // reduced cost

    n_cut_off = - _cut_off;

    for ( i = _nodes; i != _sentinel_node; i ++) {
        for ( a = i->first(), a_stop = (i + 1)->suspended(); a != a_stop; a ++) {

            rc = REDUCED_COST( i, a->head(), a );       
            if ( ( rc > _cut_off && CLOSED(a->sister()) ) ||
                 ( rc < n_cut_off && CLOSED(a) ) ) { // suspend the arc

                b = i->first();
                i->inc_first();
                exchange( a, b );
            }
        }
    }
}

int MCMF_CS2::update_epsilon()
{
    // decrease epsilon after epsilon-optimal flow is constructed;
    if ( _epsilon <= 1 ) return ( 1 );

    _epsilon = (price_t) (ceil ( (double) _epsilon / _f_scale ));
    _cut_off = _cut_off_factor * _epsilon;
    _cut_on = _cut_off * CUT_OFF_GAP;

    return ( 0 );
}

int MCMF_CS2::check_feas()
{
    if ( _check_solution == false)
        return ( 0);
    
    NODE *i;
    ARC *a, *a_stop;
    long fa;
    int ans = 1;

    for ( i = _nodes; i != _sentinel_node; i ++) {
        for ( a = i->suspended(), a_stop = (i + 1)->suspended(); a != a_stop; a ++) {
            if ( _cap[ N_ARC(a) ] > 0) {
                fa = _cap[ N_ARC(a) ] - a->rez_capacity();
                if ( fa < 0) {
                    ans = 0;
                    break;
                }
                _node_balance[ i - _nodes ] -= fa;
                _node_balance[ a->head() - _nodes ] += fa;
            }
        }
    }

    for ( i = _nodes; i != _sentinel_node; i ++) {
        if ( _node_balance[ i - _nodes ] != 0) {
            ans = 0;
            break;
        }
    }

    return ( ans);
}

int MCMF_CS2::check_cs()
{
    // check complimentary slackness;
    NODE *i;
    ARC *a, *a_stop;

    for ( i = _nodes; i != _sentinel_node; i ++) {
        for ( a = i->suspended(), a_stop = (i + 1)->suspended(); a != a_stop; a ++) {

            if ( OPEN(a) && (REDUCED_COST(i, a->head(), a) < 0) ) {
                return ( 0);
            }
        }
    }
    return(1);
}

int MCMF_CS2::check_eps_opt()
{
    NODE *i;
    ARC *a, *a_stop;

    for ( i = _nodes; i != _sentinel_node; i ++) {
        for ( a = i->suspended(), a_stop = (i + 1)->suspended(); a != a_stop; a ++) {

            if ( OPEN(a) && (REDUCED_COST(i, a->head(), a) < - _epsilon) ) {
                return ( 0);
            }
        }
    }
    return(1);
}

void MCMF_CS2::print_solution( bool print_flow_only)
{
    NODE *i;
    ARC *a;
    long ni;
    price_t cost;
    printf("\n flow: %d", _max_flow);
    printf("\n cost: %.2f\n", _flow_cost );
    for ( i = _nodes; i < _nodes + _n; i ++ ) {
        ni = N_NODE( i );
        for ( a = i->suspended(); a != (i+1)->suspended(); a ++) {
            if ( _cap[ N_ARC (a) ]  > 0 ) {
                if ( print_flow_only) {
                    int flow = _cap[ N_ARC(a) ] - a->rez_capacity();
                    if ( flow > 0) {
                        printf("%d -> %d flow:%d cost:%d\n", ni, N_NODE(a->head()),
                               flow, a->cost());
                    }
                } else {
                    printf("f %7ld %7ld %10ld\n", 
                        ni, N_NODE(a->head()), _cap[ N_ARC(a) ] - a->rez_capacity());
                }
            }
        }
    }

    // COMP_DUALS?
    if ( _comp_duals == true) { // find minimum price;
        cost = MAX_32;
        for ( i = _nodes; i != _sentinel_node; i ++) {
            cost = MIN(cost, i->price());
        }
        for ( i = _nodes; i != _sentinel_node; i ++) {
            printf("p %7ld %7.2lld\n", N_NODE(i), i->price() - cost);
        }
    }
    printf("\n");
}

void MCMF_CS2::print_graph()
{
    NODE *i;
    ARC *a;
    long ni, na;
    printf ("\nGraph: %d\n", _n);
    for ( i = _nodes; i < _nodes + _n; i ++ ) {
        ni = N_NODE( i );
        printf("\nNode %d  demand: %d", ni, i->excess());
        for ( a = i->suspended(); a != (i+1)->suspended(); a ++) {
            na = N_ARC( a );
            printf("\n {%d} %d -> %d  cap: %d  cost: %d",
                   na, ni, N_NODE(a->head()), _cap[N_ARC(a)], a->cost());
        }
    }
}

void MCMF_CS2::print_stuff()
{
    printf ("\nStuff:");
    printf ("\n _n: %d  _m: %d", _n, _m);
    printf ("\n _linf: %d  _time_for_price_in: %d", _linf, _time_for_price_in);
    printf ("\n _epsilon: %d  _dn: %d", _epsilon, _dn);
    printf ("\n _price_min: %d  _mmc: %d", _price_min, _mmc);
    printf ("\n _f_scale: %.2f  _cut_off_factor: %.2f", _f_scale, _cut_off_factor);
    printf ("\n _cut_on: %.2f  _cut_off: %.2f", _cut_on, _cut_off);
    printf ("\n _total_excess: %d  _flag_price: %d", _total_excess, _flag_price);
    printf ("\n _flag_updt: %d  _snc_max: %d", _flag_updt, _snc_max);
    printf ("\n _n_rel: %d", _n_rel);
    printf ("\n _n_ref: %d", _n_ref);
    printf ("\n _n_src: %d", _n_src);
    printf ("\n _n_push: %d", _n_push);
    printf ("\n _n_relabel: %d", _n_relabel);
    printf ("\n _n_discharge: %d", _n_discharge);
    printf ("\n _n_refine: %d", _n_refine);
    printf ("\n _n_update: %d", _n_update);
    printf ("\n _n_scan: %d", _n_scan);
    printf ("\n _n_prscan: %d", _n_prscan);
    printf ("\n _n_prscan1: %d", _n_prscan1);
    printf ("\n _n_prscan2: %d", _n_prscan2);
    printf ("\n _n_bad_pricein: %d", _n_bad_pricein);
    printf ("\n _n_bad_relabel: %d", _n_bad_relabel);
    printf ("\n _n_prefine: %d", _n_prefine);
    printf ("\n _max_cost: %d", _max_cost);
    printf ("\n _total_p: %d  _total_n: %d", _total_p, _total_n);
}

int MCMF_CS2::get_flow_of_arc( int src_id, int des_id)
{
    int this_flow = 0;
    NODE *i;
    ARC *a;
    int ni;
    for ( i = _nodes; i < _nodes + _n; i ++ ) {
        ni = N_NODE( i );
        if ( ni != src_id) continue;
        for ( a = i->suspended(); a != (i+1)->suspended(); a ++) {
            if ( N_NODE(a->head()) != des_id) continue;
            this_flow = _cap[ N_ARC(a) ] - a->rez_capacity();
        }
    }
    return this_flow;
}

void MCMF_CS2::finishup( double *objective_cost)
{
    ARC *a; // current arc
    long na; // corresponding position in capacity array
    double obj_internal = 0; // objective
    price_t cs; // actual arc cost
    long flow; // flow through an arc
    NODE *i;

    // (1) NO_ZERO_CYCLES?
    if ( _no_zero_cycles == true) {
        for ( a = _arcs; a != _sentinel_arc; a ++ ) {
            if ( a->cost() == 1) {
                assert( a->sister()->cost() == -1);
                a->set_cost( 0);
                a->sister()->set_cost( 0);
            }
        }
    }

    // (2)
    for ( a = _arcs, na = 0; a != _sentinel_arc ; a ++, na ++ ) {
        cs = a->cost() / _dn;
        if ( _cap[na]  > 0 && (flow = _cap[na] - a->rez_capacity()) != 0 )
            obj_internal += (double) cs * (double) flow;
        a->set_cost( cs);
    }

    for ( i = _nodes; i != _sentinel_node; i ++) {
        i->set_price( (i->price() / _dn));
    }

    // (3) COMP_DUALS?
    if ( _comp_duals == true) {
        compute_prices();
    }

    *objective_cost = obj_internal;
}

int MCMF_CS2::cs2( double *objective_cost)
{
    // the main calling function;
    int cc = 0; // for storing return code;

  
    // (1) update epsilon first;
    update_epsilon();


    // (2) scaling loop;
    int iter_counter = 0;

    do {
        iter_counter ++;

        int unfeasible_sol = refine();

        if ( unfeasible_sol > 0) {
            // in this case, the demand of S,T is too large; we'll decrement
            // it and re-run cs2 until a feasible solution is found;
            return 2; // unfeasible
        }

        if ( _n_ref >= PRICE_OUT_START ) // >= 1;
            price_out();

        if ( update_epsilon() ) 
            break;

        while (1) {
            if ( ! price_refine() ) 
                break;

            if ( _n_ref >= PRICE_OUT_START ) { // >= 1;
                if ( price_in() ) break; 
                if ( (cc = update_epsilon()) ) break;
            }
        }
    } while ( cc == 0 );


    // (3) finishup;
    finishup( objective_cost );

    return 0;
}

int MCMF_CS2::run_cs2()
{
    // example of flow network in DIMACS format:
    //
    //"p min 6 8
    //c min-cost flow problem with 6 nodes and 8 arcs
    //n 1 10
    //c supply of 10 at node 1
    //n 6 -10
    //c demand of 10 at node 6
    //c arc list follows
    //c arc has <tail> <head> <capacity l.b.> <capacity u.b> <cost>
    //a 1 2 0 4 1
    //a 1 3 0 8 5
    //a 2 3 0 5 0
    //a 3 5 0 10 1
    //a 5 4 0 8 0
    //a 5 6 0 8 9
    //a 4 2 0 8 1
    //a 4 6 0 8 1"
    //
    // in order to solve this flow problem we have to follow these steps:
    // 1. ctor of MCMF_CS2 // sets num of nodes and arcs
    //                     // it also calls allocate_arrays()
    // 2. call set_arc() for each arc
    // 3. call set_supply_demand_of_node() for non-transhipment nodes
    //    Part 1:
    // 4. pre_processing()
    // 5. cs2_initialize()
    //    Part 2:
    // 6. cs2()
    // 7. retreive results
    //
    // this function is basically a wrapper to implement Steps 4, 5, 6;



    // Part 1:

    // (4) ordering, etc.;
    pre_processing();


    // () CHECK_SOLUTION?
    if ( _check_solution == true) {
        _node_balance = (long long int *) calloc (_n+1, sizeof(long long int));
        for ( NODE *i = _nodes; i < _nodes + _n; i ++ ) {
            _node_balance[i - _nodes] = i->excess();
        }
    }


    // (5) initializations;
    _m = 2 * _m;
    cs2_initialize(); // works already with 2*m;


    if ( _print_ans == true ) {
        printf("c CS 4.3\n");
        printf("c nodes: %ld  arcs: %ld\n", _n, _m/2 );
        printf("c scale-factor: %f  cut-off-factor: %f\nc\n",
               _f_scale, _cut_off_factor);
    }



    // Part 2:

    // (6) run CS2;
    double objective_cost;


    int unfeasible_sol = 0;
    unfeasible_sol = cs2( &objective_cost );

    bool lower_demand_until_solution_is_found = true;
    if ( lower_demand_until_solution_is_found == true) {
        while ( unfeasible_sol > 0) {
            // () re-clear;
            clear_all_light();

            // () lower demand of S,T;
            bool solution_may_still_exist = decrement_demand_of_source_and_sink();
            if ( solution_may_still_exist == false) {
                printf("\nError: in run_cs2()");
                err_end( UNFEASIBLE ); // terminate;
            }
            // () re-initialize;
            cs2_initialize_light();

            // () re-run cs2;
            unfeasible_sol = cs2( &objective_cost );
        }
    }

    


    double t = 0.0;  
    _flow_cost = objective_cost;
    if ( _print_ans == true ) {
        printf("c time:         %15.2f    cost:       %15.0f\n", t, objective_cost);
        printf("c refines:      %10ld     discharges: %10ld\n", _n_refine, _n_discharge);
        printf("c pushes:       %10ld     relabels:   %10ld\n", _n_push, _n_relabel);
        printf("c updates:      %10ld     u-scans:    %10ld\n", _n_update, _n_scan);
        printf("c p-refines:    %10ld     r-scans:    %10ld\n", _n_prefine, _n_prscan);
        printf("c dfs-scans:    %10ld     bad-in:     %4ld  + %2ld\n",
               _n_prscan1, _n_bad_pricein, _n_bad_relabel);
    }
    

    // () CHECK_SOLUTION?
    if ( _check_solution == true ) {
        printf("c checking feasibility...\n"); 
        if ( check_feas() )
            printf("c ...OK\n");
        else
            printf("c ERROR: solution infeasible\n");
        printf("c computing prices and checking CS...\n");
        compute_prices();
        if ( check_cs() )
            printf("c ...OK\n");
        else
            printf("ERROR: CS violation\n");
    }

    // () PRINT_ANS?
    if ( _print_ans == true ) {
        print_solution();
    }




    // () cleanup; not done here because after the mcmf is solved, we still need it
    // for printing and retrieval later on; will be done instead separately
    // during every iteration, inside clear_all();
    // deallocate_arrays();

    return 0;
}

bool MCMF_CS2::decrement_demand_of_source_and_sink()
{
    // this function is introduced to facillitate the iterative
    // call of cs2 until a solution is found;
    // it's possible that I have source S
    // with say 6 outgoing arcs and sink T with 6 incomming arcs but the only
    // possible solution is only 5; in this case if I set demand of S and T as 6,-6
    // CS2 would return Error 2 (unfeasible solution);

    NODE *i;
    bool solution_may_still_exist = false;


    if ( _demand_ST > 1) {
        solution_may_still_exist = true;
        _demand_ST --;
        //printf("\n New _demand_ST: %d", _demand_ST);
        set_supply_demand_of_node( _source_id, _demand_ST);
        set_supply_demand_of_node( _sink_id, -_demand_ST);
    }


    //for ( i = _nodes; i < _nodes + _n; i ++ ) {
    //  excess_t this_excess = i->excess(); // long long int;
    //  if ( this_excess > 0) {
    //      if ( this_excess > 1) {
    //          solution_may_still_exist = true;
    //      }
    //      printf("\n new demand: %d", this_excess - 1);
    //      i->set_excess( this_excess - 1);
    //  } else if ( this_excess < 0) {
    //      i->set_excess( this_excess + 1);
    //      printf("\n new demand: %d", this_excess + 1);
    //  }
    //}

    return solution_may_still_exist;
}

