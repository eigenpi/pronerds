#include "nerds_gui.h"
#include "nerds.h"
#include <string>
#include <math.h>
#include <time.h>

using namespace std;

////////////////////////////////////////////////////////////////////////////////
//
// PRONERDS
//
////////////////////////////////////////////////////////////////////////////////

PRONERDS::PRONERDS() :
    _use_gui(false),
    _rng_seed(1)
{
    // _nx and _ny will be read from the network_file;
    _nerds_graph = 0;
    _mcmf_algo = MCMF_ALGO_EDMONDS;
    _reconfig_algo = RECONFIG_MCMF;
    _pfs_algo = BARAN_DIST_FLOW;
    _runtime_power_flow = 0.0;
    _runtime_all = 0.0;
}

bool PRONERDS::parse_command_arguments( int argc, char *argv[])
{
    bool result = true;
    int i = 0;

    // (1) network file is mandatory;
    if (argc < 2) {
        printf("\nUsage:   nerds network_file [Options...]\n");
        printf("Example: nerds tests/bus_83_11.pos\n");
        printf("Options: [-option val|...] Default\n");
        printf("\t[-seed Int] 1\n");
        printf("\t[-use_gui]\n");
        printf("\t[-pfs_algo rw | distflow] distflow\n");
        printf("\t[-reconfig_algo mcmf | baran] baran\n");
        //printf("\t[-mcmf_algo edmonds | cs2] edmonds\n");
        exit(1);
    }

    // (2) here we read in the other Options;
    _mcmf_algo = MCMF_ALGO_EDMONDS; // by default we run Edmonds as slightly faster one;
    _reconfig_algo = RECONFIG_BARAN; // by default we run the oldy Baran heuristic;
    _pfs_algo = BARAN_DIST_FLOW;
    _use_gui = false; // by default we do not show any graphics;
    i = 2;
    while ( i < argc) {

        if (strcmp (argv[i],"-mcmf_algo") == 0) {
            if (argc <= i+1) {
                printf ("Error:  -mcmf_algo option requires a string parameter.\n");
                exit (1);
            } 
            if (strcmp(argv[i+1], "edmonds") == 0) {
                _mcmf_algo = MCMF_ALGO_EDMONDS;
            } 
            else if (strcmp(argv[i+1], "cs2") == 0) {
                _mcmf_algo = MCMF_ALGO_CS2;
            } else {
                printf("Error:  -mcmf_algo must be edmonds or cs2.\n");
                exit (1);
            }
            i += 2;
            continue;
        }

        if (strcmp (argv[i],"-reconfig_algo") == 0) {
            if (argc <= i+1) {
                printf ("Error:  -reconfig_algo option requires a string parameter.\n");
                exit (1);
            } 
            if (strcmp(argv[i+1], "mcmf") == 0) {
                _reconfig_algo = RECONFIG_MCMF;
            } 
            else if (strcmp(argv[i+1], "baran") == 0) {
                _reconfig_algo = RECONFIG_BARAN;
            } else {
                printf("Error:  -reconfig_algo must be baran or mcmf.\n");
                exit (1);
            }
            i += 2;
            continue;
        }

        if (strcmp (argv[i],"-pfs_algo") == 0) {
            if (argc <= i+1) {
                printf ("Error:  -pfs_algo option requires a string parameter.\n");
                exit (1);
            } 
            if (strcmp(argv[i+1], "rw") == 0) {
                _pfs_algo = RANDOM_WALK_POWER_FLOW_PARTIAL;
            }
            else if (strcmp(argv[i+1], "rw_full") == 0) {
                _pfs_algo = RANDOM_WALK_POWER_FLOW_FULL;
            }
            else if (strcmp(argv[i+1], "distflow") == 0) {
                _pfs_algo = BARAN_DIST_FLOW;
            } else {
                printf("Error:  -pfs_algo must be distflow or rw.\n");
                exit (1);
            }
            i += 2;
            continue;
        }

        if (strcmp(argv[i],"-use_gui") == 0) {
            _use_gui = true;
            i += 1; // do not take any parameter value;
            continue;
        }

        if (strcmp(argv[i],"-seed") == 0) {
            _rng_seed = atoi( argv[i+1]);
            // just set it right away :);
            srand ( _rng_seed);
            if ((_rng_seed < 1) || (_rng_seed > INT_MAX)) { 
                printf("Error:  -seed value must be between [1 %d].\n", INT_MAX);
                exit(1); 
            }
            i += 2; 
            continue;
        }
    }
    return result;      
}

void PRONERDS::print_initial_stats( int argc, char *argv[])
{
    printf("network_file:           %s \n", argv[1]);
    printf("Nodes count:            %d \n", _nerds_graph->nodes_count());
    printf("  Feeders:              %d \n", _nerds_graph->feeders_count());
    printf("  Buses:                %d \n", _nerds_graph->nodes_count() - _nerds_graph->feeders_count());
    printf("Arcs count:             %d \n", _nerds_graph->arcs_count());
    printf("  Sectionalizing:       %d \n", _nerds_graph->sectionalizing_count());
    printf("  Ties:                 %d \n", _nerds_graph->ties_count());
    printf("RNG seed:               %d \n", _rng_seed);
    printf("Power Flow Solution:    %s \n",
        (_pfs_algo == RANDOM_WALK_POWER_FLOW_FULL ||
        _pfs_algo == RANDOM_WALK_POWER_FLOW_PARTIAL) ? "rw" : "distflow");
    printf("Reconfiguration algo:   %s \n",
        (_reconfig_algo == RECONFIG_MCMF) ? "mcmf" : "baran");
    if ( _reconfig_algo == RECONFIG_MCMF) {
        printf("MCMF algo:              %s \n",
            (_mcmf_algo == MCMF_ALGO_EDMONDS) ? "edmonds" : "cs2");
    }
}

bool PRONERDS::create_network_graph( int argc, char *argv[])
{
    // create power distribution network graph;
    // Note: this has to be called after a call to parse_command_arguments();
    bool result = true;
    int i = 0;

    // (1) network file is mandatory;
    string network_graph_file_name = argv[1];
    ifstream network_graph_ifstream;
    network_graph_ifstream.open( network_graph_file_name.c_str());
    if ( !network_graph_ifstream) {
        printf("\nError: Cannot open network file: %s\n",
            network_graph_file_name.c_str());
        exit(1);
    }

    // the power distribution network file contains the netlist of the whole
    // thing; the format should be:
    //
    //.PU 1 // meaning that data are in p.u. values
    //.V_base 11.40
    //.S_base 5.68
    //.Branch_count 16
    //.Sectionalizing_count 13
    //.Tie_count 3
    //.Nodes 16
    //.Feeders 3
    //.Feeder_node_ids 0 1 2
    //.Branch Src_bus Rec_bus R X PL_kw QL_kvar S_NS
    //0 0 3 0.075 0.1 2.0 1.6 S
    //1 3 4 0.08 0.11 3.0 1.5 S
    //2 3 5 0.09 0.18 2.0 0.8 S
    //3 5 6 0.04 0.04 1.5 1.2 S
    //...
    //.Tie switches
    //13 4 10 0.04 0.04 
    //14 9 13 0.04 0.04 
    //15 6 15 0.12 0.12 
    //.Node x_coord y_coord
    //.Grid_size 8 8
    //0 1 5
    //1 4 5
    //2 7 5
    //...

    // (2) read from the netlist file; 
    // I assume that the input files are correct for now; I do not do any
    // checking and printouts to educate the user;


    // () read in various counts;
    string my_string;


    // () find out if netlist is in p.u. values;
    bool is_data_in_pu = true; // preferred;
    long pu_value;
    network_graph_ifstream >> my_string; // skip ".PU";
    network_graph_ifstream >> pu_value;
    if ( pu_value == 0) is_data_in_pu = false;


    // () read in V_base and S_base from user;
    // initially I used same V_base and Q_base for all testcases;
    double P_base = 5084.26; // kW;
    double Q_base = 2547.32; // kvar;
    double V_base = 12.66; // kV; it was 12.66 initially for su testcase;
    double S_base = 0.001 * sqrt(P_base*P_base + Q_base*Q_base); // MVA megavoltamper;
    // OR a second set of default values:
    //double V_base = 1000; // 1000 kV;
    //double S_base = 10000; // 10000 MVA megavoltamper;

    // now I read them from file;
    double V_base_from_user = 0.0; // kW;
    double S_base_from_user = 0.0; // kvar;
    network_graph_ifstream >> my_string; // skip ".V_base";
    network_graph_ifstream >> V_base_from_user;
    network_graph_ifstream >> my_string; // skip ".S_base";
    network_graph_ifstream >> S_base_from_user;
    // overwrite default values;
    V_base = V_base_from_user;
    S_base = S_base_from_user;


    // () calculate z_pu;
    double z_pu = (S_base / (V_base * V_base));


    // () read in value of Branch_count;
    long branch_count;
    network_graph_ifstream >> my_string; // skip ".Branch_count";
    network_graph_ifstream >> branch_count;


    // () read in value of Sectionalizing_count;
    long sectionalizing_count;
    network_graph_ifstream >> my_string; // skip ".Sectionalizing_count";
    network_graph_ifstream >> sectionalizing_count;

    // () read in value of Tie_count;
    long ties_count;
    network_graph_ifstream >> my_string; // skip ".Tie_count";
    network_graph_ifstream >> ties_count;

    // read in value of Nodes;
    long nodes_count;
    network_graph_ifstream >> my_string; // skip ".Nodes";
    network_graph_ifstream >> nodes_count;
    // () increase with one nodes_count to account for the super feeder node 
    // as well; super feeder node will have the largest index; regular feeders
    // are 0,1,2,...;
    long nodes_count_extended = nodes_count + 1;

    // read in value of Feeders;
    long feeders_count;
    network_graph_ifstream >> my_string; // skip ".Feeders";
    network_graph_ifstream >> feeders_count;

    // () increase with feeders_count the number of branches; these branches
    // will connect super feeder node to all feeder nodes; each such super branch
    // will have r,x set to zero; load P,Q of each normal feeder node will be 
    // set to 0;
    long branch_count_extended = branch_count + feeders_count;

    // () create the new nerds graph object; Note: nr. of branches and nodes
    // account for the super feeder node and arcs;
    // this ctor calls also:
    // _nerds_graph->build_nodes( nodes_count)
    // _nerds_graph->build_arcs( arcs_count)
    _nerds_graph = new NERDS_GRAPH(
        branch_count_extended, sectionalizing_count, ties_count,
        nodes_count_extended, feeders_count);


    // () read in the ids of nodes that are feeder type and record that;
    // ".Feeder_node_ids 0 1 2"
    network_graph_ifstream >> my_string; // skip ".Feeder_node_ids";
    for ( long j = 0; j < feeders_count; j ++) {
        long feeder_id;
        network_graph_ifstream >> feeder_id;
        _nerds_graph->get_node( feeder_id)->set_type( FEEDER);
    }

    // () read in arcs;
    // skip this entire line: ".Branch Src_bus Rec_bus R X PL_kw QL_kvar S_NS";
    // read in arcs with or without sectionalizing switches;
    for ( long i = 0; i < 8; i ++) {
        network_graph_ifstream >> my_string;
    }
    long arcs_count = branch_count - ties_count;
    //printf("\n arcs_count: %d", arcs_count);
    for ( int j = 0; j < arcs_count; j ++) {
        long arc_id, src_id, des_id;
        double r, x, pl_kw, ql_kvar;
        network_graph_ifstream >> arc_id;
        network_graph_ifstream >> src_id;
        network_graph_ifstream >> des_id;
        network_graph_ifstream >> r;
        network_graph_ifstream >> x;
        network_graph_ifstream >> pl_kw;
        network_graph_ifstream >> ql_kvar;
        network_graph_ifstream >> my_string;
        assert( j == arc_id);
        NERDS_ARC_TYPE arc_type = NOT_SECTIONALIZING; // doesn't have sect. switch;
        if ( my_string.compare("S") == 0) {
            arc_type = SECTIONALIZING;
        }
        if ( is_data_in_pu == false) {
            r *= z_pu;
            x *= z_pu;
            pl_kw /= (1000 * S_base);
            ql_kvar /= (1000 * S_base);
        }
        _nerds_graph->set_arc( arc_id, src_id, des_id, arc_type, r, x);
        _nerds_graph->get_node( des_id)->set_pl_ql( pl_kw, ql_kvar);
        //printf("\n%d %d->%d", arc_id, src_id, des_id);
    }
    
    // read in arcs with tie switches;
    // skip" .Tie switches";
    network_graph_ifstream >> my_string;
    network_graph_ifstream >> my_string;
    for ( int j = 0; j < ties_count; j ++) {
        long arc_id, src_id, des_id;
        double r, x;
        network_graph_ifstream >> arc_id;
        network_graph_ifstream >> src_id;
        network_graph_ifstream >> des_id;
        network_graph_ifstream >> r;
        network_graph_ifstream >> x;
        if ( is_data_in_pu == false) {
            r *= z_pu;
            x *= z_pu;
        }
        _nerds_graph->set_arc( arc_id, src_id, des_id, TIE, r, x);
        _nerds_graph->add_tie(j, arc_id); // record initial tie switches;
        //printf("\n%d %d->%d", arc_id, src_id, des_id);
    }
    

    // () read coordinates of all nodes; their placement for gui display;
    // skip ".Node x_coord y_coord";
    network_graph_ifstream >> my_string;
    network_graph_ifstream >> my_string;
    network_graph_ifstream >> my_string;
    // skip ".Grid_size" and read in grid size;
    network_graph_ifstream >> my_string;
    long nx, ny;
    network_graph_ifstream >> nx;
    network_graph_ifstream >> ny;
    _nerds_graph->set_grid_size_nx_ny(nx, ny);
    // now read placement coords;
    for ( int i = 0; i < nodes_count ; i ++) {
        long node_id;
        int x_coord, y_coord;
        network_graph_ifstream >> node_id;
        network_graph_ifstream >> x_coord;
        network_graph_ifstream >> y_coord;
        _nerds_graph->set_node_coords( node_id, x_coord, y_coord);
        //printf("\n%d %d,%d", node_id, x_coord, y_coord);
    }
    

    // (3) clean up;
    network_graph_ifstream.close();


    // (4) now graph is almost built;
    
    // (a) we create a super feeder node as well and its arcs towards its children
    // that are all normal feeder nodes; this way the whole system becomes
    // a giant feder tree; DistFlow will be called for this;
    //
    long super_feeder_id = nodes_count; // i.e., nodes_count_extended - 1;
    _nerds_graph->get_node( super_feeder_id)->set_type( SUPER_FEEDER);
    long super_arc_id = branch_count;
    for ( long f_i = 0; f_i < feeders_count; f_i++) {
        _nerds_graph->set_arc( super_arc_id, super_feeder_id, f_i,
            NOT_SECTIONALIZING, 0.0, 0.0); // r = 0, x = 0;
        _nerds_graph->get_node( f_i)->set_pl_ql( 0.0, 0.0);
        super_arc_id ++;
    }

    // (b) stuff needed during sweeps for power flow solution calculation;
    _nerds_graph->compute_children_counts();
    _nerds_graph->update_parents_of_all_nodes(); // arc parent ids;
    _nerds_graph->set_host_pronerds( this); // needed for passing the gui object;
    _nerds_graph->reset_feeders_depth(); // one time only;
    _nerds_graph->compute_feeders_depth(); // initial depths;
    
    
    // (5) sanity check;
    //_nerds_graph->print_graph();
    //_nerds_graph->print_voltages();

    return result;
}

bool PRONERDS::run_DistFlow_solution_baran( bool simplified_distflow,
    bool use_super_feeder_node)
{
    // use Baran's method;
    _nerds_graph->run_DistFlow_solution_baran( simplified_distflow, use_super_feeder_node);
}

bool PRONERDS::run_PowerFlow_solution_rw( long version) // 1 or 2;
{
    // use new method based on random walks;
    if ( version == 1) {
        _nerds_graph->run_PowerFlow_solution_rw_version1();
    } else if ( version == 2) {
        _nerds_graph->run_PowerFlow_solution_rw_version2();
    }
}

bool PRONERDS::run_reconfiguration_method1_baran( bool simplified_distflow,
    POWER_FLOW_METHOD power_flow_solution) // default is BARAN_DIST_FLOW;
{
    // use Baran's method 1;
    _nerds_graph->run_reconfiguration_method1_baran(
        simplified_distflow, power_flow_solution);
}

bool PRONERDS::run_reconfiguration_method_new( bool simplified_distflow,
    POWER_FLOW_METHOD power_flow_solution) // default is BARAN_DIST_FLOW;
{
    // use new method;
    _nerds_graph->run_reconfiguration_method_new(
        simplified_distflow, power_flow_solution);
}

////////////////////////////////////////////////////////////////////////////////
//
// NERDS_GRAPH
//
////////////////////////////////////////////////////////////////////////////////

bool NERDS_GRAPH::run_DistFlow_solution_baran( bool simplified_distflow,
    bool use_super_feeder_node, bool used_during_partial_rw)
{
    // use_super_feeder_node is false by default, but is true then we 
    // create a super feeder node which will be the parent of all
    // normal feeder nodes; then every DistFlow solution will be
    // performed only with one call for one feeder tree, whose root
    // is the super-feeder node;
    // this is the so known DistFlow method for the power flow solution;
    // if simplified_distflow == true, then use the so called Simplified DistFlow
    // from Baran's paper; it's the one used actually in implementing
    // Method 1;
    // run the power flow calculations proposed in the 1989 paper
    // of Baran et al; this is an aproximation method; iterative;
    // forward and backward a couple of times;
    // for all feeder nodes call the backward and forward computations
    // iteratively until voltages converge: error is less than _delta_voltage;
    // Note: backward has to be the first iter all the time;
    //
    // if used_during_partial_rw is true (default is false), then, we
    // we do initial prob_intervals calculations as well;

    if ( use_super_feeder_node) {
        // super feeder node and arcs must have been added by this time...;
        run_super_feeder_DistFlow_solution_baran( simplified_distflow);
        return true;
    }
    bool debug_runtime = true; // runtime profiling;
    timeval time_val_1, time_val_2; // sec and microsec;
    if ( debug_runtime) { gettimeofday( &time_val_1, 0); }

    printf("\nRunning DistFlow power flow solution...");


    // (0) clean-up;
    clear_nodes();
    

    // (1) one time initial dry calculation of P,Q's from leaves towards root;
    for ( long f_i = 0; f_i < _feeders_count; f_i++) {
        NERDS_NODE *feeder_node = &_nodes[ f_i];
        initial_backward_P_Q_aggregation_feeder_tree( feeder_node);
    }
    if ( used_during_partial_rw) {
        // we do initial prob_intervals calculations as well, because
        // following updates will be done using partial random walks
        // performed only for some nodes; the others have to have already 
        // their prob_interval's calculated here;
        for ( long f_i = 0; f_i < _feeders_count; f_i++) {
            NERDS_NODE *feeder_node = &_nodes[ f_i];
            compute_probability_intervals( feeder_node); // mem alloc issue?
            feeder_node->set_lodging_type( HOME_NODE);
            feeder_node->set_lodging_cost( 1.0); // 1 V p.u.;
        }
        // debug;
        //if ( debug_runtime) {
        //  gettimeofday( &time_val_2, 0);
        //  printf ("\nRuntime compute_probability_intervals: %2.8f sec",
        //      time_val_2.tv_sec - time_val_1.tv_sec +
        //      double(time_val_2.tv_usec - time_val_1.tv_usec) / 1000000.0);
        //}
    }


    // (2) feeders are the first nodes all the time;
    long this_feeder_n_count = 0;
    for ( long f_i = 0; f_i < _feeders_count; f_i++) {
        NERDS_NODE *feeder_node = &_nodes[ f_i];
        //printf("\nFEEDER {%d/%d}", f_i, _feeders_count-1);

        // reset stuff to be used for a whole new "while" looping for
        // a new feeder id f_i;
        double sum_of_voltages_k_1 = 100; // prev iteration k-1;
        double sum_of_voltages_k = 0; // current iteration k;
        double delta_error = 100.;
        long iter_counter = 0;
        while ( delta_error > POWER_FLOW_SOLUTION_ERROR &&
                iter_counter < POWER_FLOW_ITERATION_BOUND) {

            //printf("\n iter %d B", iter_counter);
            backward_sweep( iter_counter, feeder_node, simplified_distflow);
            //backward_mixed_with_forward_sweep( iter_counter, feeder_node, simplified_distflow);
            //print_voltages();

            //printf("\n iter %d F", iter_counter);
            forward_sweep( iter_counter, feeder_node, simplified_distflow);
            //print_voltages();     

            sum_of_voltages_k = compute_sum_nodes_voltages( feeder_node,
                this_feeder_n_count); // curr iter;
            // use average Vi of all Vi's of the nodes of this feeder tree;
            sum_of_voltages_k = sum_of_voltages_k / double( max(1,this_feeder_n_count));
            //printf("\n---> vk = %.12f vk-1 = %.12f", sum_of_voltages_k,sum_of_voltages_k_1);
            // here one way is to simply take the difference between average Vi of last
            // two iterations and compare it to a hard coded number say 1E-2; but this
            // will depend on the number of nodes this feeder has; second way is to
            // see what fraction or percentage new delta is of the new average;
            delta_error = 100 * (sum_of_voltages_k - sum_of_voltages_k_1) /
                sum_of_voltages_k_1; // 100*(curr - prev)/prev;
            delta_error = fabs(delta_error); // work around oscillations too;
            sum_of_voltages_k_1 = sum_of_voltages_k; // curr will be come prev for next iter;
    
            //printf("\nIter %d: error=%.5f  losses=%.8f",iter_counter,delta_error,compute_total_losses());
            iter_counter++;
        }
    }
    // debug;
    if ( debug_runtime) {
        gettimeofday( &time_val_2, 0);
        _pronerds->add_to_runtime_power_flow(
            time_val_2.tv_sec - time_val_1.tv_sec + 
            double(time_val_2.tv_usec - time_val_1.tv_usec) / 1000000.0);
        //printf ("\nRuntime walltime of DistFlow power flow: %2.8f sec",
        //  time_val_2.tv_sec - time_val_1.tv_sec + 
        //  double(time_val_2.tv_usec - time_val_1.tv_usec) / 1000000.0);
    }
    //print_voltages();
    //print_losses();
}

bool NERDS_GRAPH::run_super_feeder_DistFlow_solution_baran( bool simplified_distflow)
{
    //
    // super feeder implementation
    //
    // create a super feeder node which will be the parent of all
    // normal feeder nodes; then every DistFlow solution will be
    // performed only with one call for one feeder tree, whose root
    // is the super-feeder node;
    // derived from run_DistFlow_solution_baran() above;

    bool debug_runtime = false; // runtime profiling;
    timeval time_val_1, time_val_2; // sec and microsec;
    if ( debug_runtime) { gettimeofday( &time_val_1, 0); }

    double sum_of_voltages_k_1 = 0; // prev iteration k-1;
    double sum_of_voltages_k = 0; // current iteration k;
    double delta_error = 1.;

    printf("\nRunning super feeder DistFlow power flow solution...");
    // (0) clean-up;
    clear_nodes();
    
    // (1) one time initial dry calculation of P,Q's from leaves towards root;
    NERDS_NODE *super_feeder_node = &_nodes[ _super_feeder_id];
    initial_backward_P_Q_aggregation_feeder_tree( super_feeder_node);

    // (2)
    long this_feeder_n_count = 0;
    long iter_counter = 0;
    while ( delta_error > POWER_FLOW_SOLUTION_ERROR &&
            iter_counter < POWER_FLOW_ITERATION_BOUND) {

        backward_sweep( iter_counter, super_feeder_node, simplified_distflow);
        //print_voltages();

        forward_sweep( iter_counter, super_feeder_node, simplified_distflow);
        //print_voltages();

        sum_of_voltages_k = compute_sum_nodes_voltages(
            super_feeder_node, this_feeder_n_count); // curr iter;
        sum_of_voltages_k = sum_of_voltages_k / double( max(1,this_feeder_n_count));
        delta_error = 100 * (sum_of_voltages_k - sum_of_voltages_k_1) /
            sum_of_voltages_k_1; // 100*(curr - prev)/prev;
        delta_error = fabs(delta_error); // work around oscillations too;
        sum_of_voltages_k_1 = sum_of_voltages_k; // curr will be come prev for next iter;
    
        //printf("\nIter %d: error=%.5f  losses=%.8f",iter_counter,delta_error,compute_total_losses());
        iter_counter++;
    }
    // debug;
    if ( debug_runtime) {
        gettimeofday( &time_val_2, 0);
        printf ("\nRuntime walltime of super feeder DistFlow power flow: %2.8f sec",
            time_val_2.tv_sec - time_val_1.tv_sec + 
            double(time_val_2.tv_usec - time_val_1.tv_usec) / 1000000.0);
    }
    //print_voltages();
    //print_losses();
}

void NERDS_GRAPH::backward_mixed_with_forward_sweep( long iter, 
    NERDS_NODE *feeder_node, bool simplified_distflow)
{
    //
    // currently not used: the whole thing becomes too much only forward updates;
    //
    // method to perform calculations from leaves towards the root 
    // of this feeder_tree; the feeder node is passed as argument, then
    // walk first towards leaves to find first leaf; start walking backward
    // till the first fork is met: here, launch a forward sweeps for all 
    // the other children of the fork's parent, and, then continue with the
    // backward walk; will this make the convergence depend on the leaf node we start?
    NERDS_NODE *upstream_node = feeder_node;
    NERDS_NODE *downstream_node = 0;
    NERDS_NODE *forward_node = 0;
    long upstream_id = upstream_node->id();
    long downstream_id = 0;
    long forward_id = 0; // id of nodes towards which we launch local forward sweeps;
    deque<long> f_queue; // forward queue;
    deque<long> b_queue; // backward queue;
    bool is_leaf_node = false;
    bool found_a_leaf_already = false;

    //printf("\nBackward sweep:");
    // (1) first, walk feeder_tree from root towards leaves to find first leaf;
    f_queue.push_back( upstream_id);
    while ( !f_queue.empty()) {

        // get node from queue and process;
        upstream_id = f_queue.back(); // get the last one: lifo;
        f_queue.pop_back(); // remove it;
        upstream_node = &_nodes[ upstream_id];
        // clean this node up and thus prepare it for the coming upstream sweep;
        upstream_node->reset_processed_children();
        // copy P Q values into k_1 (ie prev iter) storage;
        upstream_node->copy_P_Q_to_Pk_1_Qk_1();

        is_leaf_node = true;
        // look at all neighbors (adjanct nodes) and walk them thru
        // except the parent from where we came to this one;
        long neighbors_count = upstream_node->fans().size();
        long parent_id = upstream_node->parent_id();
        for ( long k = 0; k < neighbors_count; k++) {
            NERDS_ARC *arc = &_arcs[ upstream_node->fans( k)];
            if ( arc->type() != TIE) {
                downstream_id = ( arc->des_id() != upstream_id) ? 
                    arc->des_id() : arc->src_id();
                if ( downstream_id != parent_id) { // skip parent;
                    f_queue.push_back( downstream_id);
                    is_leaf_node = false; // has at least a child;
                }
            }
        }
        // found leaf?
        if ( is_leaf_node) {
            if ( found_a_leaf_already == false) {
                b_queue.push_back( upstream_id); // record only one leaf: starting point;
                found_a_leaf_already = true;
            }
            // if this is iteration 0 then Vk (current iter voltage) is
            // assigned value one; else Vk will be used as computed during
            // prev iter sweep;
            if ( iter == 0) {
                upstream_node->set_v_k( 1.0);
                upstream_node->set_P( 0.0);
                upstream_node->set_Q( 0.0);
            }
        }
    }

    // (2) second, by this time b_queue has just one leaf node;
    // now start backward/upstream mixed with local forward calculations;
    while ( !b_queue.empty()) {

        // (a) get node from queue and process;
        downstream_id = b_queue.back();
        b_queue.pop_back(); // remove it;
        downstream_node = &_nodes[ downstream_id];
        upstream_node = &_nodes[ downstream_node->parent_id()];
        upstream_id = upstream_node->id();
        //printf(" {<-%d}", downstream_id);

        // (b) calculate Pi_1 Qi_1 Vi_1 using eq. 2 from Baran paper;
        // this also records that one child of upstream_node (ie the 
        // downstream_node) has been processed as part of the upstream 
        // sweep (done by calling upstream_node->incr_processed_children());
        if ( downstream_node->type() != SUPER_FEEDER) {
            // ()
            backward_branch_equations( upstream_node, downstream_node,
                simplified_distflow, true); // is_for_lateral == true;

            // () if upstream_node is a fork node, then launch local forward
            // sweeps for the rest of children;
            bool is_fork_node = ( upstream_node->children_count() > 1);
            if ( is_fork_node) {
                long neighbors_count = upstream_node->fans().size();
                long parent_id = upstream_node->parent_id();
                for ( long k = 0; k < neighbors_count; k++) {
                    NERDS_ARC *arc = &_arcs[ upstream_node->fans( k)];
                    if ( arc->type() != TIE) {
                        forward_id = ( arc->des_id() != upstream_id) ? 
                            arc->des_id() : arc->src_id();
                        if ( forward_id != parent_id && // skip parent;
                            forward_id != downstream_id) { // skip node from where we came;
                            // launch local forward sweep:
                            // -- compute branch equations between fork node and child;
                            // -- call forward sweep for child;
                            forward_node = &_nodes[ forward_id];
                            forward_branch_equations( upstream_node, forward_node,
                                true, simplified_distflow, // is_fork_node = true;
                                true); // is_for_lateral = true;
                            forward_sweep( iter, forward_node,
                                simplified_distflow);

                            // also, accumulate this child's P,Q, albeit of prev
                            // iteration; but needed so that P,Q reflect all P,Q's
                            // during the normal forward sweep;
                            upstream_node->accumulate_P( arc->P());
                            upstream_node->accumulate_Q( arc->Q());
                        }
                    }
                }
            }

            // () do not neeed to record that one more child has been processed;
            // add upstream node to the b_queue directly to step closer towards
            // root of feeder_tree in this hybrid backward sweep;
            b_queue.push_back( upstream_node->id());    
        }
    }
}

void NERDS_GRAPH::backward_sweep( long iter, NERDS_NODE *feeder_node, 
    bool simplified_distflow)
{
    // method to perform calculations from leaves towards the root 
    // of this feeder_tree; the feeder node is passed as argument, then
    // walk first towards leaves to find them
    // once a leaf node is found, backwards (ie upstream) calculations
    // are started and propagated until lastly we get back to the 
    // feeder node; this whole routine traverses a feeder_tree twice;
    // how to improve?
    NERDS_NODE *upstream_node = feeder_node;
    NERDS_NODE *downstream_node = 0;
    long upstream_id = upstream_node->id();
    long downstream_id = 0;
    deque<long> f_queue; // forward queue;
    deque<long> b_queue; // backward queue;
    bool is_leaf_node = false;

    //printf("\nBackward sweep:");
    // (1) first, walk feeder_tree from root towards leaves to find the leaves;
    f_queue.push_back( upstream_id);
    while ( !f_queue.empty()) {

        // get node from queue and process;
        upstream_id = f_queue.back(); // get the last one: lifo;
        f_queue.pop_back(); // remove it;
        upstream_node = &_nodes[ upstream_id];
        // clean this node up and thus prepare it for the coming upstream sweep;
        upstream_node->reset_processed_children();
        // copy P Q values into k_1 (ie prev iter) storage;
        upstream_node->copy_P_Q_to_Pk_1_Qk_1();
        upstream_node->set_P_Q( 0.0, 0.0); // reset, to prepare for backwards accumulation;

        is_leaf_node = true;
        // look at all neighbors (adjanct nodes) and walk them thru
        // except the parent from where we came to this one;
        long neighbors_count = upstream_node->fans().size();
        long parent_id = upstream_node->parent_id();
        for ( long k = 0; k < neighbors_count; k++) {
            NERDS_ARC *arc = &_arcs[ upstream_node->fans( k)];
            if ( arc->type() != TIE) {
                downstream_id = ( arc->des_id() != upstream_id) ? 
                    arc->des_id() : arc->src_id();
                if ( downstream_id != parent_id) { // skip parent;
                    f_queue.push_back( downstream_id);
                    is_leaf_node = false; // has at least a child;
                }
            }
        }
        // found leaf?
        if ( is_leaf_node) {
            b_queue.push_back( upstream_id);
            // this should be commented out normally;
            if ( upstream_node->v_k() < 0.3) {
                printf("\nWarning: Power Flow numerical solution has small voltages.");
                printf("\n         This leads to instability and is due to small V_base.\n");
                exit(1);
            }
            //printf("(%d)", upstream_id);
            // if this is iteration 0 then Vk (current iter voltage) is
            // assigned value one; else Vk will be used as computed during
            // prev iter sweep;
            if ( iter == 0) {
                //upstream_node->set_v_k( 1.0);
                upstream_node->set_P( 0.0);
                upstream_node->set_Q( 0.0);
            }
        }
    }

    // (2) second, by this time b_queue has all leaves nodes of the feed-node;
    // now start backward/upstream calculations;
    while ( !b_queue.empty()) {

        // (a) get node from queue and process;
        downstream_id = b_queue.back();
        b_queue.pop_back(); // remove it;
        downstream_node = &_nodes[ downstream_id];
        upstream_node = &_nodes[ downstream_node->parent_id()];
        //printf(" {%d}", downstream_id);

        // (b) calculate Pi_1 Qi_1 Vi_1 using eq. 2 from Baran paper;
        // this also records that one child of upstream_node (ie the 
        // downstream_node) has been processed as part of the upstream 
        // sweep (done by calling upstream_node->incr_processed_children());
        if ( downstream_node->type() != SUPER_FEEDER) {
            backward_branch_equations( upstream_node, downstream_node, 
                simplified_distflow);

            // record that one more child has been processed;
            upstream_node->incr_processed_children();
            // check if all children of upstream_node were processed and
            // if so then add it to the b_queue to step closer towards
            // root of feeder_tree;
            if ( upstream_node->processed_children_all()) {
                b_queue.push_back( upstream_node->id());    
            }
        }
    }
}

void NERDS_GRAPH::forward_sweep( long iter, NERDS_NODE *feeder_node,
    bool simplified_distflow)
{
    // simplified_distflow == true, then we use simplified distflow of Baran;
    // method to perform calculations from root node of this feeder_tree; 
    // the feeder node is passed as argument, then juts walk first towards 
    // leaves;
    NERDS_NODE *upstream_node = feeder_node;
    NERDS_NODE *downstream_node = 0;
    long upstream_id = upstream_node->id();
    long downstream_id = 0;
    deque<long> f_queue; // forward queue;

    //printf("\nForward sweep:");
    // (1) simple forward walk of the tree;
    f_queue.push_back( upstream_id);
    while ( !f_queue.empty()) {

        // (a) get node from queue and process;
        upstream_id = f_queue.back(); // get the last one: lifo;
        f_queue.pop_back(); // remove it;
        upstream_node = &_nodes[ upstream_id];
        //printf(" {%d}", upstream_id);

        // (b) take every child, process, and then add it to queue;
        long neighbors_count = upstream_node->fans().size();
        bool is_fork_node = ( neighbors_count > 2);
        long parent_id = upstream_node->parent_id();
        for ( long k = 0; k < neighbors_count; k++) {
            NERDS_ARC *arc = &_arcs[ upstream_node->fans( k)];
            if ( arc->type() != TIE) { // good branch;
                downstream_id = ( arc->des_id() != upstream_id) ? 
                    arc->des_id() : arc->src_id();
                if ( downstream_id != parent_id) { // skip parent;
                    downstream_node = &_nodes[ downstream_id];
                    // calculate Pi Qi Vi using eq. 1 from Baran paper;
                    forward_branch_equations( upstream_node, downstream_node,
                        is_fork_node, simplified_distflow);
                    f_queue.push_back( downstream_id);
                }
            }
        }
    }
}

void NERDS_GRAPH::backward_branch_equations( NERDS_NODE *u_node,
    NERDS_NODE *d_node, bool simplified_distflow, bool is_for_lateral)
{
    // (1) simply use eq 2 from Baran paper;
    assert(d_node->parent_arc_id() >= 0 && d_node->parent_arc_id() < _arcs_count);
    NERDS_ARC *arc_i_1 = &_arcs[ d_node->parent_arc_id()]; // arc i-1;
    double r_i_1 = arc_i_1->r(); // ri-1;
    double x_i_1 = arc_i_1->x(); // xi-1;
    // Pi is what d_node sees downstream; for leaf node it should be zero;
    // Note: ctor of nodes sets it to one;
    double Pi = d_node->P();
    double Qi = d_node->Q();
    double Pi_prime = Pi + d_node->pl_kw();
    double Qi_prime = Qi + d_node->ql_kvar();
    double Vi_sq = ( d_node->v_k() * d_node->v_k());
    double I_sq = 0.0;
    double Pi_1 = 0.0;
    double Qi_1 = 0.0;
    double Vi_1 = 0.0;
    if ( simplified_distflow == false) { // normal DistFlow power flow;
        I_sq = ( Pi_prime*Pi_prime + Qi_prime*Qi_prime) / ( Vi_sq);
        Pi_1 = Pi + d_node->pl_kw() + r_i_1 * I_sq; // Pi-1;
        Qi_1 = Qi + d_node->ql_kvar() + x_i_1 * I_sq; // Qi-1;
        Vi_1 = sqrt( Vi_sq + 
                     2*(r_i_1 * Pi_prime + x_i_1 * Qi_prime) -
                     (r_i_1 * r_i_1 + x_i_1 * x_i_1)*I_sq); // Vi-1;
        // debug;
        //assert( Vi_sq + 2*(r_i_1 * Pi_prime + x_i_1 * Qi_prime) >= 
        //  (r_i_1 * r_i_1 + x_i_1 * x_i_1)*I_sq);
    } else {
        Pi_1 = Pi + d_node->pl_kw(); // Pi-1;
        Qi_1 = Qi + d_node->ql_kvar(); // Qi-1;
        Vi_1 = sqrt( Vi_sq + 2*(r_i_1 * Pi_prime + x_i_1 * Qi_prime)); // Vi-1;
    }
    // record Pi_1, Qi_1 on the arc between upstream and downstream node;
    // will be useful for "percentage" splitting between kids during forward
    // sweep;
    arc_i_1->set_P_Q( Pi_1, Qi_1);
    if ( is_for_lateral) {
        u_node->set_P( Pi_1);
        u_node->set_Q( Qi_1);
    } else { // used during normal/original backward sweep;
        u_node->accumulate_P( Pi_1);
        u_node->accumulate_Q( Qi_1);
    }
    //printf(" %d:Vi:%.4f Pi:%.4f Qi:%.4f",
    //  u_node->id(), u_node->v_k(), Vi_1, u_node->v_k());

    // record voltage of this iteration inside upstream node;
    // Note: here I may arrive from two diff leaves to the same node
    // and have two different voltages; which one to consider for the node?
    // compare to how we do it during forward sweep;
    u_node->set_v_k( min(Vi_1,u_node->v_k())); // EITHER set as smallest among kids...
    //u_node->set_v_k( Vi_1);    // OR set as only "most current"...
    //u_node->update_v_k( Vi_1); // OR store max between this and "most current";
}

void NERDS_GRAPH::forward_branch_equations( NERDS_NODE *u_node,
    NERDS_NODE *d_node, bool is_fork_node, bool simplified_distflow,
    bool is_for_lateral)
{
    // (1) simply use eq 1 from Baran paper;
    NERDS_ARC *arc_i = &_arcs[ d_node->parent_arc_id()]; // arc i;
    double r_i = arc_i->r(); // ri;
    double x_i = arc_i->x(); // xi;
    // Note: tricky part: at fork nodes we take a fraction of P and Q
    // calculated until now during this forward sweep; the fraction
    // is computed using prev-iter P,Q stored on the arc of each child;
    // Pk_1 is what u_node saw cumulatively downstream during prev iter;
    double Pi = 0.0;
    double Qi = 0.0;
    if ( is_for_lateral) {
        Pi = arc_i->P();
        Qi = arc_i->Q();
    } else {
        Pi = ( is_fork_node) ?
            ( u_node->P() * ( arc_i->P() / u_node->Pk_1())) : ( u_node->P());
        Qi = ( is_fork_node) ?
            ( u_node->Q() * ( arc_i->Q() / u_node->Qk_1())) : ( u_node->Q());

        //printf(" %s,%.12f,%.12f ",is_fork_node ? "y":"n",
        //     u_node->Pk_1(), u_node->P());
    }
    double Vi_sq = ( u_node->v_k() * u_node->v_k());
    double I_sq = 0.0;
    double Pi_1 = 0.0;
    double Qi_1 = 0.0;
    double Vi_1 = 0.0;
    if ( simplified_distflow == false) { // normal DistFlow power flow;
        I_sq = ( Pi*Pi + Qi*Qi) / ( Vi_sq);
        Pi_1 = Pi - d_node->pl_kw() - r_i * I_sq; // Pi+1;
        Qi_1 = Qi - d_node->ql_kvar() - x_i * I_sq; // Qi+1;
        Vi_1 = sqrt( Vi_sq -
                     2*(r_i * Pi + x_i * Qi) +
                     (r_i * r_i + x_i * x_i)*I_sq); // Vi+1;
    } else {
        Pi_1 = Pi - d_node->pl_kw(); // Pi+1;
        Qi_1 = Qi - d_node->ql_kvar(); // Qi+1;
        Vi_1 = sqrt( Vi_sq - 2*(r_i * Pi + x_i * Qi)); // Vi+1;
    }
    // before saving/overwriting, first copy old values into k_1 storage;
    // Pi_1,Qi_1 are the P,Q "flowing" downstream from this d_node; at
    // fork nodes this flowing is split between children; Note: in this case
    // (forward sweep) we do not store P,Q on arc because P,Q has meaning of
    // accumulated P,Q supplying all children of a node;
    d_node->copy_P_Q_to_Pk_1_Qk_1();
    d_node->set_P( Pi_1);
    d_node->set_Q( Qi_1);
    // record voltage of this iteration inside downstream node;
    d_node->set_v_k( Vi_1);

    //printf(" {%d-%d}:Vk:%.4f:%.4f:%.4f:%.4f ",
    //  u_node->id(),d_node->id(),u_node->v_k(),u_node->P(), Pi, Pi_1);
}

////////////////////////////////////////////////////////////////////////////////
//
// Baran Method 1
//
////////////////////////////////////////////////////////////////////////////////

bool NERDS_GRAPH::run_reconfiguration_method1_baran( bool simplified_distflow,
    POWER_FLOW_METHOD power_flow_solution) // default is BARAN_DIST_FLOW;
{
    // the outer main loop of the Baran Method 1;
    char msg[BUFFER_SIZE];
    GUI_GRAPHICS *gui = pronerds()->gui();
    double initial_losses = 1.0;

    // next arrays should be used only during collecting runtime data;
    // after that it should be commented out;
    //vector<double> losses_every_iteration;

    printf("\nRunning Baran method 1..."); printf("\n=========================");
    bool use_super_feeder_node = false;
    vector<int> nodes_for_rw_update;


    // (1) first time power flow solution;
    if ( power_flow_solution == BARAN_DIST_FLOW) {
        run_DistFlow_solution_baran( simplified_distflow, use_super_feeder_node);
    } else if ( power_flow_solution == RANDOM_WALK_POWER_FLOW_FULL) {
        run_PowerFlow_solution_rw_version1();
    } else if ( power_flow_solution == RANDOM_WALK_POWER_FLOW_PARTIAL) {
        // initial full/complete power flow solution is DistFlow because it's
        // faster than PowerFlow RW Full;
        run_DistFlow_solution_baran( simplified_distflow,use_super_feeder_node,true);
    }
    
    // (2) main algo;
    long iter_count = 0;
    long cumulated_be_count = 0;
    double this_loss_reduction = 1.0;
    double cumulated_loss_reduction=0.0, prev_cumulated_loss_reduction=0.0;
    int tie_i_k_2 = -1, sec_i_k_2 = -1; // two iteration ago; used to detect oscillation;
    int tie_i_k_1 = -1, sec_i_k_1 = -1; // one iteration ago;
    bool oscillation = false;
    while ( oscillation == false &&
            this_loss_reduction > RECONFIG_LOSS_REDUCTION_ERROR && // 1E-16
            iter_count < RECONFIG_ITERATION_BOUND) { // 100
        printf("\n================= ITER: %d", iter_count + 1);


        // () clear up storage; related to when we use rw based PowerFlow only;
        if ( power_flow_solution == RANDOM_WALK_POWER_FLOW_PARTIAL) {
            nodes_for_rw_update.clear();
            _sketch_affected_feeders.clear();
        }
        

        // () entertain user
        if ( _pronerds->use_gui()) {
            sprintf( msg, "Configuration before ITER %d", iter_count);
            gui->update_screen( PRIORITY_MAJOR, msg, NODES);
        }


        // () record some info;
        //print_graph(); // debug;
        if ( iter_count == 0) {
            initial_losses = compute_total_losses(); // record it;
        }
        bool debug_mode = false;
        if ( debug_mode) {
            double final_losses = compute_total_losses();
            printf ("\nInitial total losses:   %.16f p.u.", initial_losses);
            printf ("\nTotal losses at iter %d: %.16f p.u.", iter_count, final_losses);
            printf ("\nTotal losses red. after %d iters using Baran method 1: %.2f \%\n", 
                    iter_count, 100*(final_losses - initial_losses)/initial_losses);
        }


        // () search and record losses for all 1st order branch-exchanges;
        explore_branch_exchanges_baran1();


        // () sort exchanges;
        sort( _sketch_exchange_losses.begin(), _sketch_exchange_losses.end());
        //print_sketch_exchange_losses(); // debug;


        // () realize/implement the best exchange; needs graph restructuring;
        this_loss_reduction = 0.0;
        TUPLE_THREE &best_ex = _sketch_exchange_losses[ _sketch_exchange_losses.size() - 1];
        int tie_i = best_ex.tie_i();
        int sec_i = best_ex.sec_i();
        bool performed_at_least_one_be = false;
        if ( best_ex.loss() > 0) {
            if ( iter_count > 1 && tie_i_k_2 == tie_i && sec_i_k_2 == sec_i) {
                oscillation = true;
                break;
            }
            _arcs[ tie_i].set_type( SECTIONALIZING);
            _arcs[ sec_i].set_type( TIE);
            // record the node ids of the affected arcs; they will be used for
            // possible localized updates when running the random-walks partial 
            // technique;
            // Note: currently we may have duplicates recorded here!
            if ( power_flow_solution == RANDOM_WALK_POWER_FLOW_PARTIAL) {
                int s1 = _arcs[tie_i].src_id();
                int d1 = _arcs[tie_i].des_id();
                int s2 = _arcs[sec_i].src_id();
                int d2 = _arcs[sec_i].des_id();
                // also add feeder id's of these node to the affected list;
                add_feeder_id_of_this_node_to_affected_feeders( _nodes[ s1 ].feeder_id());
                nodes_for_rw_update.push_back( s1);
                add_feeder_id_of_this_node_to_affected_feeders( _nodes[ d1 ].feeder_id());
                nodes_for_rw_update.push_back( d1);
                if ( s2 != s1 && s2 != d1) {
                    add_feeder_id_of_this_node_to_affected_feeders( _nodes[ s2 ].feeder_id());
                    nodes_for_rw_update.push_back( s2);
                }
                if ( d2 != s1 && d2 != d1) {                
                    add_feeder_id_of_this_node_to_affected_feeders( _nodes[ d2 ].feeder_id());
                    nodes_for_rw_update.push_back( d2);
                }
                //printf ("\n----> Number of affected feeders: %d", _sketch_affected_feeders.size());
            }
            this_loss_reduction += best_ex.loss();
            // locally update the children counts of affected nodes;
            // () update the separate list of tie switches as well;
            for ( long i = 0; i < _ties_count; i ++) {
                if ( _ties[ i] == tie_i) { 
                    _ties[ i] = sec_i; 
                }
            }
            printf("\nClose:{%d}  Open:{%d}  Loss reduction: %.16f", tie_i, sec_i, best_ex.loss());
            performed_at_least_one_be = true;
        } else {
            printf("\n---> No branch exchange can lead to loss reduction!");
        }
        
        // () reflect the switch actions in the network graph;
        compute_children_counts();
        update_parents_of_all_nodes();
        //print_graph(); // debug;


        // () updates and cleanup to prepare for next iter;
        this_loss_reduction = best_ex.loss();
        _sketch_exchange_losses.clear();


        tie_i_k_2 = tie_i_k_1; sec_i_k_2 = sec_i_k_1;
        tie_i_k_1 = tie_i; sec_i_k_1 = sec_i;
        iter_count ++;
        cumulated_be_count ++;

        
        // () recompute the power flow solution after current iter reconfiguration:
        // "closing of tie switch" and "opening of sectionalizing switch";
        // basically prepare it for next iter;
        if ( performed_at_least_one_be == true) {
            if ( power_flow_solution == BARAN_DIST_FLOW) {
                run_DistFlow_solution_baran( simplified_distflow, use_super_feeder_node);
            } else if ( power_flow_solution == RANDOM_WALK_POWER_FLOW_FULL) {
                run_PowerFlow_solution_rw_version1();
            } else if ( power_flow_solution == RANDOM_WALK_POWER_FLOW_PARTIAL) {
                // Note: we count how many cumulated branch exchanges we
                // have (in Baran's method we do only one be in each iter)
                // and after every other 100 be's we do a full DistFlow for
                // better accuracy;
                if ( cumulated_be_count % 50 != 0) {
                    // PowerFlow updates based on partial RW's;
                    run_PowerFlow_solution_rw_update_subset_of_nodes( nodes_for_rw_update);
                } else { // do a full DistFlow from time to time;
                    run_DistFlow_solution_baran( simplified_distflow,use_super_feeder_node,true);   
                }
            }
        }
        

        // () another stopping criterion;
        // Note: for Baran's reconfiguration I do not stop if loss reduction 
        // between two consecutive iterations is less than 0.1%; that is because
        // the total cumulated loss reduction would be too small; if I did however,
        // runtime would be shorter;
        //prev_cumulated_loss_reduction = cumulated_loss_reduction;
        //cumulated_loss_reduction += this_loss_reduction;
        //if ( (cumulated_loss_reduction - prev_cumulated_loss_reduction)/
        //   prev_cumulated_loss_reduction < 0.001) {
        //  break;
        //}


        // debug and collecting data for every iteration; normaly this should
        // be commented out;
        //double current_losses = compute_total_losses();
        //losses_every_iteration.push_back(
        //  100*(current_losses - initial_losses)/initial_losses);
    }

    // (2) we do need to recompute the power flow solution after the last 
    // iter of the loop above;
    // Note: RW based PFS is slightly different than DistFlow, hence
    // estimation errors accumulate and the actual loss computation may not
    // be the actual one; so, here we do need one last DistFlow call;
    run_DistFlow_solution_baran( simplified_distflow, use_super_feeder_node);


    double final_losses = compute_total_losses();
    printf ("\n\nInitial total losses: %.16f p.u.", initial_losses);
    printf ("\nFinal total losses:   %.16f p.u.", final_losses);
    printf ("\nTotal losses reduction using Baran Method 1: %.2f \%\n", 
        100*(final_losses - initial_losses)/initial_losses);

    // debug and data colection; normaly this should be commented out;
    //for ( long j = 0; j < losses_every_iteration.size(); j++) {
    //  if ( j == 0) {
    //      printf("\n%d %.2f",j+1, -losses_every_iteration[j]);
    //  } else {
    //      printf("\n%d %.2f",j+1, -(losses_every_iteration[j] - losses_every_iteration[j-1]));
    //  }
    //}
}

bool NERDS_GRAPH::explore_branch_exchanges_baran1()
{
    // inner loop of  Baran Method 1;
    // take every tie switch/arc and explore branch exchange;
    // 1. if the two nodes of the tie are in the same feeder_tree then
    // the loop has to be found by walking backward the tree till
    // common ancestor node is found;
    // 2. if the two nodes of the tie are in different feeder_trees
    // then the super-node-zero is their ancestor;

    // () find all candidate branch exchanges and their loss reduction
    // or increase;
    // Note: we should record only those exchanges that lead to positive 
    // loss reduction;
    int tie_i = 0;
    for ( long i = 0; i < _ties_count; i ++) {
        tie_i = _ties[ i];
        NERDS_ARC *tie_arc = &_arcs[ tie_i];

        // () find the loop as two sets of arcs from nodes of tie switch
        // backwards to common ancestor;
        find_loop_as_left_right_arc_sets( tie_i);
        //print_left_right_loop( tie_i); // debug;

        // () estimate loss reduction due to branch exchanges in both
        // directions/sides of the tie switch; L->R R->L;
        // next call basically evaluates branch exchanges and records
        // them as 3-tuples in _sketch_exchange_losses;
        // here we can use either only exploration of 1st order
        // neighboring branches or all branches in the loop;
        estimate_loss_reduction_dueto_branch_exchanges_baran1_first_order( tie_i);
        //estimate_loss_reduction_dueto_branch_exchanges_baran1_entire_loop( tie_i);
    }
}

void NERDS_GRAPH::find_loop_as_left_right_arc_sets(int tie_i)
{
    // find the loop as two sets of arcs from nodes of tie switch
    // backwards to common ancestor;
    // update the info about this new loop inside _sketch_left_arc_ids and
    // _sketch_right_arc_ids; first elements in these sketch arrays
    // are the 1st order neighboring arcs; last ones in these arrays
    // are the arcs children of the common ancestor node;

    // () cleanup;
    _sketch_left_arc_ids.clear();
    _sketch_right_arc_ids.clear();

    NERDS_ARC *tie_arc = &_arcs[ tie_i];
    deque<long> b_queue; // backward queue;

    // () left array = left part of the loop;
    NERDS_NODE *downstream_node = 0;
    long downstream_id = tie_arc->src_id();
    b_queue.push_back( downstream_id);
    while ( !b_queue.empty()) {

        downstream_id = b_queue.back(); // get node from queue;
        b_queue.pop_back(); // remove it;
        downstream_node = &_nodes[ downstream_id];
        // add parent to queue;
        if ( downstream_node->type() != FEEDER) {
            // record this arc going backward towards root of feeder_tree;
            _sketch_left_arc_ids.push_back( downstream_node->parent_arc_id());
            b_queue.push_back( downstream_node->parent_id());
        }
    }

    // () right array = right part of the loop;
    downstream_id = tie_arc->des_id(); // the other side of the tie switch;
    b_queue.push_back( downstream_id);
    while ( !b_queue.empty()) {

        downstream_id = b_queue.back();
        b_queue.pop_back();
        downstream_node = &_nodes[ downstream_id];
        // add parent to queue;
        if ( downstream_node->type() != FEEDER) {
            // record this arc going backward towards root of feeder_tree;
            _sketch_right_arc_ids.push_back( downstream_node->parent_arc_id());
            b_queue.push_back( downstream_node->parent_id());
        }
    }

    // () the left right arrays may have common arcs in case their common ancestor
    // node was lower than the root node; remove the commmong arcs;
    long left_i = _sketch_left_arc_ids.size() - 1;
    long right_i = _sketch_right_arc_ids.size() - 1;
    while ( left_i >= 0 && right_i >= 0 &&
        _sketch_left_arc_ids[ left_i] == _sketch_right_arc_ids[ right_i]) {
        _sketch_left_arc_ids.pop_back();
        _sketch_right_arc_ids.pop_back();
        left_i --;
        right_i --;
    }
}

void NERDS_GRAPH::estimate_loss_reduction_dueto_branch_exchanges_baran1_first_order(int tie_i)
{
    // implement estimation eq 9 in Baran paper, part of his Method 1;
    // "left" are ids of arcs between nodes {o,... k-1,k} 
    // "right" are ids of arcs between nodes {o,... n-1,n,k}
    double delta_pl_bk = 0.; // loss reduction;
    double tr = 0.0; // sum of all r's in the loop;
    double drp = 0., drq = 0.;
    double sum_rP_left = 0.;
    double sum_rQ_left = 0.;
    double sum_rP_right = 0.;
    double sum_rQ_right = 0.;
    double r = 0.0;

    int left_count = _sketch_left_arc_ids.size();
    for ( long i = 0; i < left_count; i ++) {
        NERDS_ARC *arc = &_arcs[ _sketch_left_arc_ids[ i]];
        r = arc->r();
        sum_rP_left += r * arc->P();
        sum_rQ_left += r * arc->Q();        
        tr += r;
    }
    int right_count = _sketch_right_arc_ids.size();
    for ( long i = 0; i < right_count; i ++) {
        NERDS_ARC *arc = &_arcs[ _sketch_right_arc_ids[ i]];
        r = arc->r();
        sum_rP_right += r * arc->P();
        sum_rQ_right += r * arc->Q();       
        tr += r;
    }
    drp = sum_rP_left - sum_rP_right;
    drq = sum_rQ_left - sum_rQ_right;
    // Note: eq 9 has be computed carefully with respect to what is left and
    // what is right; when we investigate opening switches from the right "arm"
    // of the loop, then drp and drq have to be calculated considering that
    // "left is right and right is left";
    double drp_lrrl = -drp; // i.e., sum_rP_right - sum_rP_left;
    double drq_lrrl = -drq; // i.e., sum_rQ_right - sum_rQ_left;

    // (1) load transfer from Left TO Right: is realized by closing tie switch
    // and opening the first neighbor branch of it from the left array;
    int m_i = 0;
    NERDS_ARC *sectionalizing_arc = 0;
    double Pm = 0.0, Qm = 0.0;
    if ( _sketch_left_arc_ids.size() > 0) {
        m_i = _sketch_left_arc_ids[ 0];
        sectionalizing_arc = &_arcs[ m_i];
        Pm = sectionalizing_arc->P();
        Qm = sectionalizing_arc->Q();
        delta_pl_bk = 2*drp*Pm + 2*drq*Qm - tr*(Pm*Pm + Qm*Qm);
        // record the pair of switches(tie,sectionalizing) and its loss if realized;
        _sketch_exchange_losses.push_back( TUPLE_THREE( tie_i, m_i, delta_pl_bk));
    }
    
    // (2) load transfer from Right TO Left: is realized by closing tie switch
    // and opening the first neighbor branch of it from the right_a array;
    if ( _sketch_right_arc_ids.size() > 0) {
        m_i = _sketch_right_arc_ids[ 0];
        sectionalizing_arc = &_arcs[ m_i];
        Pm = sectionalizing_arc->P();
        Qm = sectionalizing_arc->Q();
        // here we have to use "left is right and right is left" values;
        delta_pl_bk = 2*drp_lrrl*Pm + 2*drq_lrrl*Qm - tr*(Pm*Pm + Qm*Qm);
        // record the pair of switches(tie,sectionalizing) and its loss if realized;
        _sketch_exchange_losses.push_back( TUPLE_THREE( tie_i, m_i, delta_pl_bk));
    }
}

void NERDS_GRAPH::estimate_loss_reduction_dueto_branch_exchanges_baran1_entire_loop(int tie_i)
{
    // implement estimation eq 9 in Baran paper, part of his Method 1;
    // "left" are ids of arcs between nodes {o,... k-1,k} 
    // "right" are ids of arcs between nodes {o,... n-1,n,k}
    double delta_pl_bk = 0.; // loss reduction;
    double tr = 0.; // sum of all r's in the loop;
    double drp = 0., drq = 0.;
    double sum_rP_left = 0.;
    double sum_rQ_left = 0.;
    double sum_rP_right = 0.;
    double sum_rQ_right = 0.;
    double r = 0.0;
    vector<TUPLE_THREE> losses_array;

    int left_count = _sketch_left_arc_ids.size();
    for ( long i = 0; i < left_count; i ++) {
        NERDS_ARC *arc = &_arcs[ _sketch_left_arc_ids[ i]];
        r = arc->r();
        sum_rP_left += r * arc->P();
        sum_rQ_left += r * arc->Q();        
        tr += r;
    }
    int right_count = _sketch_right_arc_ids.size();
    for ( long i = 0; i < right_count; i ++) {
        NERDS_ARC *arc = &_arcs[ _sketch_right_arc_ids[ i]];
        r = arc->r();
        sum_rP_right += r * arc->P();
        sum_rQ_right += r * arc->Q();       
        tr += r;
    }
    drp = sum_rP_left - sum_rP_right;
    drq = sum_rQ_left - sum_rQ_right;
    // Note: eq 9 has be computed carefully with respect to what is left and
    // what is right; when we investigate opening switches from the right "arm"
    // of the loop, then drp and drq have to be calculated considering that
    // "left is right and right is left";
    double drp_lrrl = -drp; // i.e., sum_rP_right - sum_rP_left;
    double drq_lrrl = -drq; // i.e., sum_rQ_right - sum_rQ_left;
    
    // (1) load transfer from Left TO Right: is realized by closing tie switch
    // and opening a branch from the left array;
    int m_i = 0;
    NERDS_ARC *sectionalizing_arc = 0;
    double Pm = 0.0, Qm = 0.0;
    for ( long i = 0; i < left_count; i ++) {
        m_i = _sketch_left_arc_ids[ i];
        sectionalizing_arc = &_arcs[ m_i];
        Pm = sectionalizing_arc->P();
        Qm = sectionalizing_arc->Q();
        delta_pl_bk = 2*drp*Pm + 2*drq*Qm - tr*(Pm*Pm + Qm*Qm);
        // record the pair of switches(tie,sectionalizing) and its loss if 
        // realized; of this loop;
        losses_array.push_back( TUPLE_THREE( tie_i, m_i, delta_pl_bk));
    }
    // (2) load transfer from Right TO Left: is realized by closing tie switch
    // and opening a branch from the right_a array;
    for ( long i = 0; i < right_count; i ++) {
        m_i = _sketch_right_arc_ids[ i];
        sectionalizing_arc = &_arcs[ m_i];
        Pm = sectionalizing_arc->P();
        Qm = sectionalizing_arc->Q();
        // here we have to use "left is right and right is left" values;
        delta_pl_bk = 2*drp_lrrl*Pm + 2*drq_lrrl*Qm - tr*(Pm*Pm + Qm*Qm);
        // record the pair of switches(tie,sectionalizing) and its loss if 
        // realized; of this loop;
        losses_array.push_back( TUPLE_THREE( tie_i, m_i, delta_pl_bk));
    }

    // () sort 'em all;
    sort( losses_array.begin(), losses_array.end());

    // debug;
    //printf("\n\nLoop exchange losses");
    //int ex_count = losses_array.size();
    //for ( long i = 0; i < ex_count; i++) {
    //  printf("\n %d -> %d  %.16f", losses_array[ i].tie_i(), losses_array[ i].sec_i(), losses_array[ i].loss());
    //}

    // () record in the sketch array of the big loop iteration only the largest
    // losss reduction involving branches of this loop;
    if ( losses_array.size() > 0) {
        _sketch_exchange_losses.push_back( losses_array[ losses_array.size() - 1]);
    }
}

////////////////////////////////////////////////////////////////////////////////
//
// New reconfiguration method: idea is based on min cost max flow formulation;
//
////////////////////////////////////////////////////////////////////////////////

bool NERDS_GRAPH::run_reconfiguration_method_new( bool simplified_distflow,
    POWER_FLOW_METHOD power_flow_solution) // default is BARAN_DIST_FLOW;
{
    // the outer main loop of the new min cost max flow based method;
    char msg[BUFFER_SIZE];
    GUI_GRAPHICS *gui = pronerds()->gui();
    double initial_losses = 1.0;

    // next arrays should be used only during collecting runtime data;
    // after that it should be commented out;
    //vector<double> losses_every_iteration;

    printf("\nRunning MCMF-based Method..."); printf("\n=========================");
    bool use_super_feeder_node = false;
    vector<int> nodes_for_rw_update;


    // (0) one time preparation of the generic indexing of the generic flow graph
    // and of the interface mcmf;
    // Note: _interface_mcmf will be cleared and recycled during every iteration;
    int max_outside_node_id = 
        initialize_indexing_of_generic_flow_graph();
    // calculate bound for max number of nodes in the flow graph; memory
    // will be allocated only as needed; this is 2F+S+4;
    int bound_nn = 2 * _feeders_count + _ties_count + 2;
    //printf("\n bound_nn: %d\n", bound_nn);
    // engine_typ can be: MCMF_ALGO_EDMONDS or MCMF_ALGO_CS2;
    MCMF_ALGORITHM engine_type = pronerds()->mcmf_algo();

    initialize_interface_mcmf( engine_type, max_outside_node_id, bound_nn);
    

    // (1) first time power flow solution;
    if ( power_flow_solution == BARAN_DIST_FLOW) {
        run_DistFlow_solution_baran( simplified_distflow, use_super_feeder_node);
    } else if ( power_flow_solution == RANDOM_WALK_POWER_FLOW_FULL) {
        run_PowerFlow_solution_rw_version1();
    } else if ( power_flow_solution == RANDOM_WALK_POWER_FLOW_PARTIAL) {
        // initial full/complete power flow solution is DistFlow because it's
        // faster than PowerFlow RW Full;
        run_DistFlow_solution_baran( simplified_distflow, use_super_feeder_node, true);
    }


    // (2) main algo;
    // sketch_exchange_losses_ids used to retreive mcmf solution; will store
    // indices of branch exchanges to be performed from the sketch array
    // _sketch_exchange_losses, taht stores all possible exchanges in current
    // iter flow graph;
    vector<int> sketch_exchange_losses_ids; 
    
    long iter_count = 0;
    long cumulated_be_count = 0;
    double this_loss_reduction = 1.0;
    double cumulated_loss_reduction=0.0, prev_cumulated_loss_reduction=0.0;
    int tie_i_k_2 = -1, sec_i_k_2 = -1; // two iteration ago; used to detect oscillation;
    int tie_i_k_1 = -1, sec_i_k_1 = -1; // one iteration ago;
    bool oscillation = false;
    while ( oscillation == false &&
            this_loss_reduction > RECONFIG_LOSS_REDUCTION_ERROR && // 1E-16
            iter_count < RECONFIG_ITERATION_BOUND) { // 120
        printf("\n================= ITER: %d", iter_count + 1);


        // () clear up storage; related to when we use rw based PowerFlow only;
        if ( power_flow_solution == RANDOM_WALK_POWER_FLOW_PARTIAL) {
            nodes_for_rw_update.clear();
        }


        // () entertain user
        if ( _pronerds->use_gui()) {
            sprintf( msg, "Configuration before ITER %d", iter_count);
            gui->update_screen( PRIORITY_MAJOR, msg, NODES);
        }


        // () record some info;
        //print_graph(); // debug;
        if ( iter_count == 0) {
            initial_losses = compute_total_losses(); // record it;
        }
        bool debug_mode = false;
        if ( debug_mode) {
            double final_losses = compute_total_losses();
            printf ("\nInitial total losses:   %.16f p.u.", initial_losses);
            printf ("\nTotal losses at iter %d: %.16f p.u.", iter_count, final_losses);
            printf ("\nTotal losses red. after %d iters using MCMF-based method: %.2f \%\n", 
                    iter_count, 100*(final_losses - initial_losses)/initial_losses);
        }


        // () search and record losses for all 1st order branch-exchanges;
        explore_branch_exchanges_mcmf(); // using new min cost max flow;
        // by this time _sketch_exchange_losses should contain branch-exchanges
        // that lead to positive loss reduction involving any tie switch; 
        // they also mean load-transfer between feeders that we would like to make;
        //print_sketch_exchange_losses(); // debug;


        // () build the mcmf flow graph, and solve the mcmf problem;
        populate_interface_mcmf_with_flow_graph_info();
        

        //_interface_mcmf.print_interface(); // debug;      
        //printf("\nBuilding and solving MCMF problem...");
        _interface_mcmf.run_mcmf( engine_type);
        //_interface_mcmf.print_flow( engine_type); // debug;


        // Note: the next call will first clear sketch_exchange_losses_ids,
        // before retreiving the solution of the mcmf problem;
        _interface_mcmf.get_branch_exhanges_to_be_realized(
            sketch_exchange_losses_ids, engine_type);

        
        // () realize/implement the branch-exchanges as dictated by the mcmf
        // solution; needs graph restructuring;
        this_loss_reduction = 0.0;
        int ex_count = sketch_exchange_losses_ids.size();
        // (a) monitor possible oscillations, when only one branch exchange 
        // is performed;
        if ( ex_count == 1) {
            TUPLE_THREE &best_ex = 
                _sketch_exchange_losses[ sketch_exchange_losses_ids[ 0] ];
            int tie_i = best_ex.tie_i();
            int sec_i = best_ex.sec_i();
            if ( iter_count > 1 && tie_i_k_2 == tie_i && sec_i_k_2 == sec_i) {
                oscillation = true;
                break;
            }
            tie_i_k_2 = tie_i_k_1; sec_i_k_2 = sec_i_k_1;
            tie_i_k_1 = tie_i; sec_i_k_1 = sec_i;
        }
        // (b) perform branch exchange;
        bool performed_at_least_one_be = false;
        for ( int be_i = 0; be_i < ex_count; be_i ++) {
            TUPLE_THREE &best_ex = 
                _sketch_exchange_losses[ sketch_exchange_losses_ids[be_i] ];
            int tie_i = best_ex.tie_i();
            int sec_i = best_ex.sec_i();
            _arcs[ tie_i].set_type( SECTIONALIZING);
            _arcs[ sec_i].set_type( TIE);
            // record the node ids of the affected arcs; they will be used for
            // possible localized updates when running the random-walks partial 
            // technique;
            // Note: currently we may have duplicates recorded here!
            if ( power_flow_solution == RANDOM_WALK_POWER_FLOW_PARTIAL) {
                int s1 = _arcs[tie_i].src_id();
                int d1 = _arcs[tie_i].des_id();
                int s2 = _arcs[sec_i].src_id();
                int d2 = _arcs[sec_i].des_id();
                // also add feeder id's of these node to the affected list;
                add_feeder_id_of_this_node_to_affected_feeders( _nodes[ s1 ].feeder_id());
                nodes_for_rw_update.push_back( s1);
                add_feeder_id_of_this_node_to_affected_feeders( _nodes[ d1 ].feeder_id());
                nodes_for_rw_update.push_back( d1);
                if ( s2 != s1 && s2 != d1) {
                    add_feeder_id_of_this_node_to_affected_feeders( _nodes[ s2 ].feeder_id());
                    nodes_for_rw_update.push_back( s2);
                }
                if ( d2 != s1 && d2 != d1) {                
                    add_feeder_id_of_this_node_to_affected_feeders( _nodes[ d2 ].feeder_id());
                    nodes_for_rw_update.push_back( d2);
                }
                //printf ("\n----> Number of affected feeders: %d", _sketch_affected_feeders.size());
            }
            this_loss_reduction += best_ex.loss();
            // update the separate list of tie switches as well;
            for ( long i = 0; i < _ties_count; i ++) {
                if ( _ties[ i] == tie_i) { 
                    _ties[ i] = sec_i; 
                }
            }
            printf("\nClose:{%d}  Open:{%d}  Loss reduction: %.16f", tie_i, sec_i, best_ex.loss());
            performed_at_least_one_be = true;
        }

        
        // () reflect the switch actions in the network graph;
        compute_children_counts();
        update_parents_of_all_nodes();


        // () updates and cleanup to prepare for next iter;
        _sketch_exchange_losses.clear();
        _sketch_feeder_pairs.clear(); // used by new method only;
        // deallocate memory arrays of MCMF_ALGO_CS2 engine if used;
        if ( engine_type == MCMF_ALGO_CS2 && !_interface_mcmf.is_empty()) {
            _interface_mcmf.deallocate_arrays_mcmf_cs2();
        }

        iter_count ++;
        cumulated_be_count ++;


        // () recompute the power flow solution after current iter reconfiguration:
        // "closing of tie switch" and "opening of sectionalizing switch";
        // basically prepare it for next iter;
        // do PFS updates only if we performed at least on branch exchange;
        if ( performed_at_least_one_be == true) {
            if ( power_flow_solution == BARAN_DIST_FLOW) {
                run_DistFlow_solution_baran( simplified_distflow, use_super_feeder_node);
            } else if ( power_flow_solution == RANDOM_WALK_POWER_FLOW_FULL) {
                run_PowerFlow_solution_rw_version1();
            } else if ( power_flow_solution == RANDOM_WALK_POWER_FLOW_PARTIAL) {
                // Note: we count how many cumulated branch exchanges we
                // have (in Baran's method we do only one be in each iter)
                // and after every other 100 be's we do a full DistFlow for
                // better accuracy;
                if ( cumulated_be_count % 50 != 0) {
                    // PowerFlow updates based on partial RW's;
                    run_PowerFlow_solution_rw_update_subset_of_nodes( nodes_for_rw_update);
                } else { // do a full DistFlow from time to time;
                    run_DistFlow_solution_baran( simplified_distflow,use_super_feeder_node,true);   
                }
            }
        }


        // () another stopping criterion;
        // if between two iterations there is less than 0.1% improvement in 
        // cumulated loss reductions, then stop;
        prev_cumulated_loss_reduction = cumulated_loss_reduction;
        cumulated_loss_reduction += this_loss_reduction;
        if ( (cumulated_loss_reduction - prev_cumulated_loss_reduction)/
             prev_cumulated_loss_reduction < 0.001) {
            //oscillation = true;
            break;
        }


        // debug and collecting data for every iteration; normaly this should
        // be commented out;
        //double current_losses = compute_total_losses();
        //losses_every_iteration.push_back(
        //  100*(current_losses - initial_losses)/initial_losses);
    }

    // (3) we do need to recompute the power flow solution after the last 
    // iter of the loop above;
    // Note: RW based PFS is slightly different than DistFlow, hence
    // estimation errors accumulate and the actual loss computation may not
    // be the actual one; so, here we do need one last DistFlow call;
    run_DistFlow_solution_baran( simplified_distflow, use_super_feeder_node);


    double final_losses = compute_total_losses();
    printf ("\n\nInitial total losses: %.16f p.u.", initial_losses);
    printf ("\nFinal total losses:   %.16f p.u.", final_losses);
    printf ("\nTotal losses reduction using new MCMF-based method: %.2f \%\n", 
        100*(final_losses - initial_losses)/initial_losses);

    // (4) deallocate memory of the Edmunds solver if was used;
    if ( engine_type == MCMF_ALGO_EDMONDS) {
        _interface_mcmf.free_arrays_mcmf_edmonds();
    }


    // debug and data colection; normaly this should be commented out;
    //for ( long j = 0; j < losses_every_iteration.size(); j++) {
    //  if ( j == 0) {
    //      printf("\n%d %.2f",j+1, -losses_every_iteration[j]);
    //  } else {
    //      printf("\n%d %.2f",j+1, -(losses_every_iteration[j] - losses_every_iteration[j-1]));
    //  }
    //}
}

int NERDS_GRAPH::initialize_indexing_of_generic_flow_graph()
{
    // creates one time only the generic indexing of the generic flow graph
    // whose subgraphs will constitute flow graphs at particular iterations;
    int counter = 2; // 0,1 are s and t;
    // add donor feeder nodes;
    for ( int i = 0; i < _feeders_count; i ++) {
        _sketch_donors_id.push_back( counter);
        counter ++;
    }
    for ( int i = 0; i < _feeders_count; i ++) {
        _sketch_acceptors_id.push_back( counter);
        counter ++;
    }
    for ( int i = 0; i < _ties_count; i ++) {
        _sketch_switches_id.push_back( counter);
        counter ++;
    }
    return counter;
}

void NERDS_GRAPH::initialize_interface_mcmf(
    MCMF_ALGORITHM engine_type,
    int max_outside_node_id, int bound_nn)
{
    // Note: this should be called only one time as it does 
    // memory allocation for the Edmonds engine;
    // this will also do some mem alloc and clearing;

    // (1) set bound_nn for efficient and memory allocation!
    if ( engine_type == MCMF_ALGO_EDMONDS) {
        _interface_mcmf.set_edmonds_bound_NN( bound_nn);
    }
    else if ( engine_type == MCMF_ALGO_CS2) {

    }

    // (2) others;
    _interface_mcmf.set_max_outside_node_id( max_outside_node_id);
}

void NERDS_GRAPH::populate_interface_mcmf_with_flow_graph_info()
{
    // put into interface itself (that is recycled) info about vertices, edges
    // capacity and cost for the flow graph;
    // Note: s and t will have node id as the two largest (additionally added)
    // node-of-power-system-graph-id +1 and +2;


    // (1) use this_flow_graph_counter to continue counting/indexing "switch
    // type" nodes in the generic flow graph; this local counter will
    // basically count all switch vertices of the actual mcmf problem;
    int this_flow_graph_counter = 2 * _sketch_donors_id.size() + 2;
    

    // (2) clear first;
    _interface_mcmf.clear_all();

    // get the max among all loss reductions; will be used for cost
    // computations for the arcs in the flow graph;
    int ex_count = _sketch_exchange_losses.size();
    double max_loss = 0.0;
    for ( int k = 0; k < ex_count; k++) {
        double loss = _sketch_exchange_losses[ k].loss();
        if ( loss > max_loss) max_loss = loss;
    }

    
    // (3) for every branch exchange stored in _sketch_exchange_losses
    // we have to create a path of four edges between s and t:
    // 1: s -> donor feeder
    // 2: donor feeder -> switch node
    // 3: switch node -> acceptor feeder
    // 4: acceptor feeder -> t
    // first we create edges 2,3 and then edges 1,4 because edges 1,4 are
    // common to many paths;

    _interface_mcmf.set_is_empty( true);
    int tie_i, sec_i;
    int donor_node_id, switch_node_id, acceptor_node_id;
    // (3.a)
    for ( int k = 0; k < ex_count; k++) {
        double loss = _sketch_exchange_losses[ k].loss();
        int cost = max(1, int(100*(max_loss - loss)/max_loss));
        tie_i = _sketch_exchange_losses[ k].tie_i();
        sec_i = _sketch_exchange_losses[ k].sec_i();
        donor_node_id = _sketch_donors_id[ _sketch_feeder_pairs[ k].from_i()];
        acceptor_node_id = _sketch_acceptors_id[ _sketch_feeder_pairs[ k].to_i()];
        switch_node_id = this_flow_graph_counter;
        this_flow_graph_counter ++; // prepare for the potentially next one;

        _interface_mcmf.add_arc(
            donor_node_id, switch_node_id, cost, // cost inv prop with loss;
            DONOR_FEEDER_VERTEX, -1,
            SWITCH_VERTEX, k); // k is the index of br ex in _sketch_exchange_losses;
        //printf("\n d:%d -> s:%d cost: %d", donor_node_id,switch_node_id,cost);

        _interface_mcmf.add_arc(
            switch_node_id, acceptor_node_id, 0, // cost 0;
            SWITCH_VERTEX, -1, 
            ACCEPTOR_FEEDER_VERTEX, -1);
        //printf("\n s:%d -> a:%d cost: %d", switch_node_id,acceptor_node_id,cost);
        _interface_mcmf.set_is_empty( false);
    }
    if ( _interface_mcmf.is_empty()) {
        //printf("\n---> Interface MCMF is empty! No MCMF problem!");
        return;
    }
    // (3.b)
    for ( int i = 0; i < _feeders_count; i ++) {
        donor_node_id = _sketch_donors_id[ i];
        if ( _interface_mcmf.is_node_added( donor_node_id)) {
            _interface_mcmf.inc_num_vertices_donors(); // record how many donor nodes we have;
            _interface_mcmf.add_arc(
                0, donor_node_id, // source S id in outside world is 0;
                0, // cost 0;
                SOURCE_S_VERTEX, -1, DONOR_FEEDER_VERTEX, -1);
            //printf("\n S:%d -> d:%d", 0, donor_node_id);
        }
        acceptor_node_id = _sketch_acceptors_id[ i];
        if ( _interface_mcmf.is_node_added( acceptor_node_id)) {
            _interface_mcmf.inc_num_vertices_acceptors(); // record how many acceptor nodes we have;
            _interface_mcmf.add_arc(
                acceptor_node_id, 1, // sink T id in outside world is 0;
                0, // cost 0;
                ACCEPTOR_FEEDER_VERTEX, -1, SINK_T_VERTEX, -1);
            //printf("\n a:%d -> T:%d", acceptor_node_id, 1);
        }
    }
}

bool NERDS_GRAPH::explore_branch_exchanges_mcmf()
{
    // inner loop of new method; uses min cost max flow idea;

    // () find all candidate branch exchanges and their loss reduction
    // or increase; we record only those exchanges that lead to positive 
    // loss reduction; for each such branch exchange we create node and
    // edges in the flow graph;
    int tie_i = 0;
    for ( long i = 0; i < _ties_count; i ++) {
        tie_i = _ties[ i];
        NERDS_ARC *tie_arc = &_arcs[ tie_i];

        // () find the loop as two sets of arcs from nodes of tie switch
        // backwards to common ancestor;
        find_loop_as_left_right_arc_sets( tie_i);

        // () estimate loss reduction due to branch exchanges in both
        // directions/sides of the tie switch; L->R R->L; for the positive
        // loss reduction we will augment the flow graph accordingly;
        estimate_loss_reduction_dueto_branch_exchanges_new_first_order( tie_i);
        //estimate_loss_reduction_dueto_branch_exchanges_new_entire_loop( tie_i);
    }
}

void NERDS_GRAPH::estimate_loss_reduction_dueto_branch_exchanges_new_first_order(int tie_i)
{
    // "new" meaning for the new/proposed method;
    // this is for tie switch id "tie_i";
    // closing tie_i switch and opening sec switch with id from the loop will
    // possibly lead to loss reduction by transfering load from one feeder to another;
    // when that is the case we create four edges between s and t nodes in the flow
    // graph;
    // loss reduction is computed using the estimation eq 9 in Baran paper;
    // "left" are ids of arcs between nodes {o,... k-1,k} 
    // "right" are ids of arcs between nodes {o,... n-1,n,k}
    double delta_pl_bk = 0.; // loss reduction;
    double tr = 0.; // sum of all r's in the loop;
    double drp = 0., drq = 0.;
    double sum_rP_left = 0.;
    double sum_rQ_left = 0.;
    double sum_rP_right = 0.;
    double sum_rQ_right = 0.;
    double r = 0.0;

    int left_count = _sketch_left_arc_ids.size();
    for ( long i = 0; i < left_count; i ++) {
        NERDS_ARC *arc = &_arcs[ _sketch_left_arc_ids[ i]];
        r = arc->r();
        sum_rP_left += r * arc->P();
        sum_rQ_left += r * arc->Q();        
        tr += r;
    }
    int right_count = _sketch_right_arc_ids.size();
    for ( long i = 0; i < right_count; i ++) {
        NERDS_ARC *arc = &_arcs[ _sketch_right_arc_ids[ i]];
        r = arc->r();
        sum_rP_right += r * arc->P();
        sum_rQ_right += r * arc->Q();       
        tr += r;
    }
    drp = sum_rP_left - sum_rP_right;
    drq = sum_rQ_left - sum_rQ_right;
    double drp_lrrl = -drp; // i.e., sum_rP_right - sum_rP_left;
    double drq_lrrl = -drq; // i.e., sum_rQ_right - sum_rQ_left;
    
    // (1) load transfer from Left TO Right: is realized by closing tie switch
    // and opening the first neighbor branch of it from the left array;
    double max_delta_pl_bk = -INT_MAX;
    int max_m_i = 0;
    int m_i = 0;
    NERDS_ARC *sectionalizing_arc = 0;
    double Pm = 0.0, Qm = 0.0;
    if ( _sketch_left_arc_ids.size() > 0) {
        m_i = _sketch_left_arc_ids[ 0];
        sectionalizing_arc = &_arcs[ m_i];
        Pm = sectionalizing_arc->P();
        Qm = sectionalizing_arc->Q();
        delta_pl_bk = 2*drp*Pm + 2*drq*Qm - tr*(Pm*Pm + Qm*Qm);

        //printf("\n %d -> %d  %.16f", tie_i, m_i, delta_pl_bk);
        if ( delta_pl_bk > max_delta_pl_bk) {
            max_delta_pl_bk = delta_pl_bk;
            max_m_i = m_i;
        }
    }
    
    // (2) load transfer from Right TO Left: is realized by closing tie switch
    // and opening the first neighbor branch of it from the right_a array;
    if ( _sketch_right_arc_ids.size() > 0) {
        m_i = _sketch_right_arc_ids[ 0];
        sectionalizing_arc = &_arcs[ m_i];
        Pm = sectionalizing_arc->P();
        Qm = sectionalizing_arc->Q();
        // here we have to use "left is right and right is left" values;
        delta_pl_bk = 2*drp_lrrl*Pm + 2*drq_lrrl*Qm - tr*(Pm*Pm + Qm*Qm);

        //printf("\n %d -> %d  %.16f", tie_i, m_i, delta_pl_bk);
        if ( delta_pl_bk > max_delta_pl_bk) {
            max_delta_pl_bk = delta_pl_bk;
            max_m_i = m_i;
        }
    }

    // (3) record the pair of switches(tie,sectionalizing) with max and positive 
    // loss reduction - if any and if realized;
    // record also the "from" and "to" feeder ids; will be used during 
    // flow graph construction using interface_mcmf;
    // Note: the sec switch that will be opened is the switch owned by
    // the feeder "from", which will have load transfered to the "to"
    // feeder, whose id can be found using the tie switch;
    int from_feeder_i = _nodes[ _arcs[max_m_i].src_id() ].feeder_id();
    int temp_i = _nodes[ _arcs[tie_i].src_id() ].feeder_id();
    int to_feeder_i = (temp_i != from_feeder_i) ? temp_i : _nodes[ _arcs[tie_i].des_id() ].feeder_id();
    // assert( from_feeder_i != to_feeder_i);
    if ( max_delta_pl_bk > 0) { // there is loss reduction;
        _sketch_exchange_losses.push_back( TUPLE_THREE( tie_i, max_m_i, max_delta_pl_bk));
        _sketch_feeder_pairs.push_back( PAIR_TWO( from_feeder_i, to_feeder_i));
    }
}

void NERDS_GRAPH::estimate_loss_reduction_dueto_branch_exchanges_new_entire_loop(int tie_i)
{
    double delta_pl_bk = 0.; // loss reduction;
    double tr = 0.; // sum of all r's in the loop;
    double drp = 0., drq = 0.;
    double sum_rP_left = 0.;
    double sum_rQ_left = 0.;
    double sum_rP_right = 0.;
    double sum_rQ_right = 0.;
    double r = 0.0;
    vector<TUPLE_THREE> losses_array;

    int left_count = _sketch_left_arc_ids.size();
    for ( long i = 0; i < left_count; i ++) {
        NERDS_ARC *arc = &_arcs[ _sketch_left_arc_ids[ i]];
        r = arc->r();
        sum_rP_left += r * arc->P();
        sum_rQ_left += r * arc->Q();        
        tr += r;
    }
    int right_count = _sketch_right_arc_ids.size();
    for ( long i = 0; i < right_count; i ++) {
        NERDS_ARC *arc = &_arcs[ _sketch_right_arc_ids[ i]];
        r = arc->r();
        sum_rP_right += r * arc->P();
        sum_rQ_right += r * arc->Q();       
        tr += r;
    }
    drp = sum_rP_left - sum_rP_right;
    drq = sum_rQ_left - sum_rQ_right;
    double drp_lrrl = -drp; // i.e., sum_rP_right - sum_rP_left;
    double drq_lrrl = -drq; // i.e., sum_rQ_right - sum_rQ_left;
    
    // (1) load transfer from Left TO Right: is realized by closing tie switch
    // and opening a branch from the left array;
    double max_delta_pl_bk = -INT_MAX;
    int max_m_i = 0;
    int m_i = 0;
    NERDS_ARC *sectionalizing_arc = 0;
    double Pm = 0.0, Qm = 0.0;
    for ( long i = 0; i < left_count; i ++) {
        m_i = _sketch_left_arc_ids[ i];
        sectionalizing_arc = &_arcs[ m_i];
        Pm = sectionalizing_arc->P();
        Qm = sectionalizing_arc->Q();
        delta_pl_bk = 2*drp*Pm + 2*drq*Qm - tr*(Pm*Pm + Qm*Qm);
        // record the pair of switches(tie,sectionalizing) and its loss if 
        // realized; of this loop;
        losses_array.push_back( TUPLE_THREE( tie_i, m_i, delta_pl_bk));
    }
    // (2) load transfer from Right TO Left: is realized by closing tie switch
    // and opening a branch from the right_a array;
    for ( long i = 0; i < right_count; i ++) {
        m_i = _sketch_right_arc_ids[ i];
        sectionalizing_arc = &_arcs[ m_i];
        Pm = sectionalizing_arc->P();
        Qm = sectionalizing_arc->Q();
        // here we have to use "left is right and right is left" values;
        delta_pl_bk = 2*drp_lrrl*Pm + 2*drq_lrrl*Qm - tr*(Pm*Pm + Qm*Qm);
        // record the pair of switches(tie,sectionalizing) and its loss if 
        // realized; of this loop;
        losses_array.push_back( TUPLE_THREE( tie_i, m_i, delta_pl_bk));
    }

    // (3) sort 'em all;
    sort( losses_array.begin(), losses_array.end());
    if ( losses_array.size() > 0) {
        max_delta_pl_bk = losses_array[ losses_array.size()-1].loss();
        max_m_i = losses_array[ losses_array.size()-1].sec_i();
    }

    // (4) record the pair of switches(tie,sectionalizing) with max and positive 
    // loss reduction - if any and if realized;
    // record also the "from" and "to" feeder ids; will be used during 
    // flow graph construction using interface_mcmf;
    // Note: the sec switch that will be opened is the switch owned by
    // the feeder "from", which will have load transfered to the "to"
    // feeder, whose id can be found using the tie switch;


    int from_feeder_i = _nodes[ _arcs[max_m_i].src_id() ].feeder_id();
    int temp_i = _nodes[ _arcs[tie_i].src_id() ].feeder_id();
    int to_feeder_i = (temp_i != from_feeder_i) ? temp_i : _nodes[ _arcs[tie_i].des_id() ].feeder_id();
    // assert( from_feeder_i != to_feeder_i);
    if ( max_delta_pl_bk > 0) { // there is loss reduction;
        _sketch_exchange_losses.push_back( TUPLE_THREE( tie_i, max_m_i, max_delta_pl_bk));
        _sketch_feeder_pairs.push_back( PAIR_TWO( from_feeder_i, to_feeder_i));
    }
}

////////////////////////////////////////////////////////////////////////////////
//
// more utilities;
//
////////////////////////////////////////////////////////////////////////////////

double NERDS_GRAPH::compute_total_losses_of_feeder_tree( NERDS_NODE *feeder_node)
{
    // compute sum of all losses in all arcs of this feeder_node tree;
    // done using a forward walk thru the tree;
    double result = 0.;
    NERDS_NODE *upstream_node = feeder_node;
    long upstream_id = upstream_node->id();
    long downstream_id = 0;
    double Pi = 0., Qi = 0., Vi = 1., ri = 0., loss_i = 0.;
    deque<long> f_queue; // forward queue;

    // (1) simple forward walk of the tree;
    f_queue.push_back( upstream_id);
    while ( !f_queue.empty()) {

        // (a) get node from queue and process;
        upstream_id = f_queue.back(); // get the last one: lifo;
        f_queue.pop_back(); // remove it;
        upstream_node = &_nodes[ upstream_id];

        // (b) take every child, process, and then add it to queue;
        long neighbors_count = upstream_node->fans().size();
        long parent_id = upstream_node->parent_id();
        Vi = upstream_node->v_k(); // voltage of parent node;
        for ( long k = 0; k < neighbors_count; k++) {
            NERDS_ARC *arc = &_arcs[ upstream_node->fans( k)];
            if ( arc->type() != TIE) { // good branch;
                downstream_id = ( arc->des_id() != upstream_id) ? 
                    arc->des_id() : arc->src_id();
                if ( downstream_id != parent_id) { // skip parent;
                    // calculate loss of this branch;
                    arc->get_P_Q( Pi, Qi); // P, Q thru this branch;
                    ri = arc->r();
                    loss_i = ri * ( (Pi*Pi + Qi*Qi)/(Vi*Vi));
                    result += loss_i; // accumulate;
                    // add the downstream node to queue;
                    f_queue.push_back( downstream_id);
                }
            }
        }
    }
    return result;
}

double NERDS_GRAPH::compute_total_losses()
{
    // done similarly to how print_losses is done: visiting each arc;
    // could be done also by calling compute_total_losses_of_feeder_tree()
    // for each of the feeder-nodes;
    double result = 0.;
    double Pi = 0., Qi = 0., Vi = 1., ri = 0., loss_i = 0.;
    int src_id = 0, des_id = 0;
    for ( long i = 0; i < _arcs_count; i ++) {
        if ( _arcs[i].type() != TIE) {
            _arcs[i].get_P_Q( Pi, Qi); // P, Q thru this branch;
            ri = _arcs[i].r();
            src_id = _arcs[ i].src_id();
            des_id = _arcs[ i].des_id();            
            if ( _nodes[ des_id].parent_id() == src_id) {
                Vi = _nodes[ src_id].v_k();
            } else {
                Vi = _nodes[ des_id].v_k();
            }
            assert( Vi > 0);
            loss_i = ri * ( (Pi*Pi + Qi*Qi)/(Vi*Vi));
            result += loss_i; // accumulate;
        }
    }
    return result;
}

void NERDS_GRAPH::update_parents_of_all_nodes()
{
    // walk through every feeder_tree and update parent node id of
    // each node; first reset the "processed" flag that is used for
    // this purpose;
    for ( long i = 0; i < _arcs_count; i++) {
        _arcs[ i].set_processed( false);
    }
    // "take a walk gentlemen" like pirates;
    // when we did not have super feeder node we used to do it
    // for every regular feeder node individually:
    // for ( long f_i = 0; f_i < _feeders_count; f_i++) {...}
    NERDS_NODE *super_feeder_node = &_nodes[ _super_feeder_id];
    super_feeder_node->set_parent_arc_id( -1); // super feeder node does not have parent arcs;
    super_feeder_node->set_parent_id( _super_feeder_id); // feeder node is its own parent;
    super_feeder_node->set_feeder_id( _super_feeder_id);
    update_parents_of_all_nodes_in_feeder_tree( super_feeder_node);
}
void NERDS_GRAPH::update_parents_of_all_nodes_in_feeder_tree( NERDS_NODE *feeder_node)
{
    int feeder_id = feeder_node->id();
    NERDS_NODE *from_node = feeder_node;
    long from_id = from_node->id();
    long to_id = 0;
    deque<long> f_queue; 
    f_queue.push_back( from_id);

    while ( !f_queue.empty()) {
        from_id = f_queue.back(); // get the last one: lifo;
        f_queue.pop_back(); // remove it;
        from_node = &_nodes[ from_id];
        long neighbors_count = from_node->fans().size();
        for ( long k = 0; k < neighbors_count; k++) {
            NERDS_ARC *arc = &_arcs[ from_node->fans( k)];
            if ( arc->type() != TIE) {
                to_id = ( arc->des_id() != from_id) ? arc->des_id() : arc->src_id();
                // this is the only place where the "processed" flag is used, 
                // because here is where we do need to update the parents, after
                // some structural changes took place, like closing and opening switches;
                if ( arc->processed() == false) {
                    f_queue.push_back( to_id);
                    arc->set_processed( true);
                    _nodes[ to_id].set_parent_arc_id( arc->id());
                    _nodes[ to_id].set_parent_id( from_id);
                    feeder_id = (_nodes[ to_id].type() == FEEDER) ? to_id : from_node->feeder_id();
                    _nodes[ to_id].set_feeder_id( feeder_id);
                    _nodes[ to_id].set_depth( from_node->depth() + 1);                  
                } else { // reset it to default: _processed = false;
                    arc->set_processed( false);
                }
            }
        }
    }
}

double NERDS_GRAPH::compute_sum_nodes_voltages( NERDS_NODE *feeder_node,
    long &tree_nodes_count)
{
    // simply add all voltages of all nodes for the tree rooted by
    // the given feeder node; use only for feeder nodes!
    // sum up voltages of iteration k+1 as being the latest ones;
    NERDS_NODE *upstream_node = feeder_node;
    double result = upstream_node->v_k();
    long upstream_id = upstream_node->id();
    long downstream_id = 0;
    deque<long> f_queue; 
    f_queue.push_back( upstream_id);
    tree_nodes_count = 0; // reset;
    
    while ( !f_queue.empty()) {
        upstream_id = f_queue.back(); // get the last one: lifo;
        f_queue.pop_back(); // remove it;
        upstream_node = &_nodes[ upstream_id];
        result += upstream_node->v_k(); // accumulate the voltage of this node;
        tree_nodes_count ++;
        
        long neighbors_count = upstream_node->fans().size();
        long parent_id = upstream_node->parent_id();
        for ( long k = 0; k < neighbors_count; k++) {
            NERDS_ARC *arc = &_arcs[ upstream_node->fans( k)];
            if ( arc->type() != TIE) {
                downstream_id = ( arc->des_id() != upstream_id) ? arc->des_id() : arc->src_id();
                if ( downstream_id != parent_id) {
                    f_queue.push_back( downstream_id);
                }
            }
        }
    }
    return result;
}

void NERDS_GRAPH::initial_backward_P_Q_aggregation_feeder_tree( NERDS_NODE *feeder_node)
{
    // estimate initially P,Q of all nodes by aggregating all P,Q
    // of downstream nodes; called only once at the beginning;
    NERDS_NODE *upstream_node = feeder_node;
    NERDS_NODE *downstream_node = 0;
    long upstream_id = upstream_node->id();
    long downstream_id = 0;
    deque<long> f_queue; // forward queue;
    deque<long> b_queue; // backward queue;
    bool is_leaf_node = false;

    // (1) first, walk feeder_tree from root towards leaves to find the leaves;
    f_queue.push_back( upstream_id);
    while ( !f_queue.empty()) {

        // get node from queue and process;
        upstream_id = f_queue.back(); // get the last one: lifo;
        f_queue.pop_back(); // remove it;
        upstream_node = &_nodes[ upstream_id];
        // clean this node up and thus prepare it for following sweeps;
        upstream_node->reset_processed_children();
        // copy P Q values into k_1 (ie prev iter) storage;
        upstream_node->copy_P_Q_to_Pk_1_Qk_1();
        upstream_node->set_P_Q( 0.0, 0.0); // reset, to prepare for backwards accumulation;

        is_leaf_node = true;
        // look at all neighbors (adjanct nodes) and walk them thru
        // except the parent from where we came to this one;
        long neighbors_count = upstream_node->fans().size();
        long parent_id = upstream_node->parent_id();
        for ( long k = 0; k < neighbors_count; k++) {
            NERDS_ARC *arc = &_arcs[ upstream_node->fans( k)];
            if ( arc->type() != TIE) {
                downstream_id = ( arc->des_id() != upstream_id) ? 
                    arc->des_id() : arc->src_id();
                if ( downstream_id != parent_id) { // skip parent;
                    f_queue.push_back( downstream_id);
                    is_leaf_node = false; // has at least a child;
                }
            }
        }
        // found leaf?
        if ( is_leaf_node) {
            b_queue.push_back( upstream_id);
            upstream_node->set_P( 0.0);
            upstream_node->set_Q( 0.0);
        }
    }

    // (2) second, by this time b_queue has all leaves nodes of the feeder-node;
    // now start backward/upstream aggregation;
    while ( !b_queue.empty()) {

        // (a) get node from queue and process;
        downstream_id = b_queue.back();
        b_queue.pop_back(); // remove it;
        downstream_node = &_nodes[ downstream_id];
        upstream_node = &_nodes[ downstream_node->parent_id()];

        // (b) accumulate P,Q;
        if ( downstream_node->type() != SUPER_FEEDER) {
            double temp_P = downstream_node->pl_kw() + downstream_node->P();
            double temp_Q = downstream_node->ql_kvar() + downstream_node->Q();
            NERDS_ARC *arc = &_arcs[ downstream_node->parent_arc_id()];
            arc->set_P_Q( temp_P, temp_Q);
            upstream_node->accumulate_P( temp_P);
            upstream_node->accumulate_Q( temp_Q);
            // record that one more child has been processed;
            upstream_node->incr_processed_children();
            // check if all children of upstream_node were processed and
            // if so then add it to the b_queue to step closer towards
            // root of feeder_tree;
            if ( upstream_node->processed_children_all()) {
                b_queue.push_back( upstream_node->id());    
            }
        }
    }
}

void NERDS_GRAPH::compute_feeders_depth()
{
    // visit all leaves of any feeder and record the maximum depth;
    // forward sweep to get all leaves; Note: depends on the fact that 
    // parent_id's are already updated;
    deque<long> f_queue; // forward queue;

    for ( long f_i = 0; f_i < _feeders_count; f_i++) {
        NERDS_NODE *feeder_node = &_nodes[ f_i];
        NERDS_NODE *upstream_node = feeder_node;
        NERDS_NODE *downstream_node = 0;
        long upstream_id = upstream_node->id();
        long downstream_id = 0;
        int this_feeder_max_depth = -1;
        bool is_leaf_node = false;

        // walk feeder_tree from root towards leaves to find the leaves;
        f_queue.push_back( upstream_id);
        while ( !f_queue.empty()) {
            // get node from queue and process;
            upstream_id = f_queue.back(); // get the last one: lifo;
            f_queue.pop_back(); // remove it;
            upstream_node = &_nodes[ upstream_id];
            is_leaf_node = true;
            // look at all neighbors (adjanct nodes) and walk them thru
            // except the parent from where we came to this one;
            long neighbors_count = upstream_node->fans().size();
            long parent_id = upstream_node->parent_id();
            for ( long k = 0; k < neighbors_count; k++) {
                NERDS_ARC *arc = &_arcs[ upstream_node->fans( k)];
                if ( arc->type() != TIE) {
                    downstream_id = ( arc->des_id() != upstream_id) ? 
                        arc->des_id() : arc->src_id();
                    if ( downstream_id != parent_id) { // skip parent;
                        f_queue.push_back( downstream_id);
                        is_leaf_node = false; // has at least a child;
                    }
                }
            }
            // found leaf;
            if ( is_leaf_node) {
                if ( this_feeder_max_depth <= upstream_node->depth()) {
                    this_feeder_max_depth = upstream_node->depth();
                }
            }
        }
        // now that this feeder has been processed, record its max depth;
        _feeders_depth[ f_i] = this_feeder_max_depth;
    }
}

////////////////////////////////////////////////////////////////////////////////
//
// debugs and prints;
//
////////////////////////////////////////////////////////////////////////////////

void NERDS_GRAPH::testing_distflow( bool simplified_distflow)
{
    // this is using hard coded sets f switches to close and to open;
    // these are solutions proposed/found by various papers; purpose
    // is to see how much loss reduction is indeed achieved; can use to
    // calibrate my implementation of DistFlow...;
    // (1)
    //print_graph(); // debug;
    bool use_super_feeder_node = true;
    run_DistFlow_solution_baran( simplified_distflow, use_super_feeder_node);
    print_voltages();
    double initial_losses = compute_total_losses(); // record it;
    // (2)
    vector<int> switches_to_close;
    vector<int> switches_to_open;
    int testcase = 1;
    // (a) civanlar
    if ( testcase == 1) {
        switches_to_close.push_back( 13);
        switches_to_close.push_back( 14);
        switches_to_open.push_back( 7);
        switches_to_open.push_back( 6);
    }
    // (b) bus_33
    else if ( testcase == 2) {
        // random
        //switches_to_close.push_back( 34);
        //switches_to_open.push_back( 10);
        //switches_to_close.push_back( 35);
        //switches_to_open.push_back( 31);
        // malaysia solution
        switches_to_close.push_back( 34);
        switches_to_open.push_back( 8);
        switches_to_close.push_back( 33);
        switches_to_open.push_back( 13);
        switches_to_close.push_back( 32);
        switches_to_open.push_back( 5);
        switches_to_close.push_back( 35);
        switches_to_open.push_back( 31);    
    }
    // (a) su
    else if ( testcase == 3) {
        switches_to_close.push_back( 83);
        switches_to_close.push_back( 84);
        switches_to_close.push_back( 86);
        switches_to_close.push_back( 87);
        switches_to_close.push_back( 90);
        switches_to_close.push_back( 92);
        switches_to_close.push_back( 93);
        switches_to_close.push_back( 94);
        switches_to_close.push_back( 95);
        switches_to_open.push_back( 54);
        switches_to_open.push_back( 6);
        switches_to_open.push_back( 71);
        switches_to_open.push_back( 12);
        switches_to_open.push_back( 82);
        switches_to_open.push_back( 38);
        switches_to_open.push_back( 33);
        switches_to_open.push_back( 40);
        switches_to_open.push_back( 61);
    }
    // (3)
    int count = switches_to_close.size();
    for ( int s_i = 0; s_i < count; s_i++) {
        int tie_i = switches_to_close[ s_i];
        int sec_i = switches_to_open[ s_i];
        _arcs[ tie_i].set_type( SECTIONALIZING);
        _arcs[ sec_i].set_type( TIE);
        // update the separate list of tie switches as well;
        for ( long i = 0; i < _ties_count; i++) {
            if ( _ties[ i] == tie_i) { 
                _ties[ i] = sec_i; 
            }
        }
        printf("\nClose:{%d}  Open:{%d}", tie_i, sec_i);
    }
    // (4) reflect the switch actions in the network graph;
    compute_children_counts();
    update_parents_of_all_nodes();

    // (5)
    run_DistFlow_solution_baran( simplified_distflow, use_super_feeder_node);
    double final_losses = compute_total_losses();
    printf ("\n\nInitial total losses: %.16f p.u.", initial_losses);
    printf ("\nFinal total losses:   %.16f p.u.", final_losses);
    printf ("\nTotal losses reduction using hard-coded solution: %.2f \%\n", 
        100*(final_losses - initial_losses)/initial_losses);
}

void NERDS_GRAPH::testing_hitting_times( bool simplified_distflow)
{
    // here I collect agerage numbers for the number of stepts till a RW
    // terminates;
    // runtime variation with the number of nodes updated;


    // (1) first run the initial DistFlow to compute once all voltages;
    run_DistFlow_solution_baran(
        simplified_distflow, 
        false, // use_super_feeder_node
        true); // used_during_partial_rw
    //print_graph(); // debug;

    // (2) then, run the RW based PowerFlow partially on a certain number of
    // nodes;
    timeval time_val_1, time_val_2; // sec and microsec;
    vector<int> nodes_for_rw_update;
    _sketch_affected_feeders.clear();

    nodes_for_rw_update.push_back( 45);
    add_feeder_id_of_this_node_to_affected_feeders( _nodes[ 15 ].feeder_id());

    gettimeofday( &time_val_1, 0);
    run_PowerFlow_solution_rw_update_subset_of_nodes( nodes_for_rw_update);
    gettimeofday( &time_val_2, 0);
    printf ("\n>>> Runtime of partial RW based PowerFlow: %2.8f sec",
        time_val_2.tv_sec - time_val_1.tv_sec + 
        double(time_val_2.tv_usec - time_val_1.tv_usec) / 1000000.0);

    //print_voltages();
}

void NERDS_GRAPH::print_tie_switches()
{
    // for user's entertainment;
    printf("Final tie switches: ");
    for ( long j = 0; j < _arcs_count; j ++) {
        if ( _arcs[j].type() == TIE)
            printf(" %d", j);
    }
    printf("\n");
}

void NERDS_GRAPH::print_graph()
{
    // debug;
    printf("\nNetwork graph");
    printf("\nNodes: %d", _nodes.size());
    for ( long i = 0; i < _nodes_count; i ++) {
        printf("\n%d %s (%d,%d) fid:%d par:%d par_arc:%d kids:%d PL:%.12f QL:%.12f Vk-1:%.2f Vk:%.2f",
               i, (_nodes[i].lodging_type() == HOME_NODE ? "H" : "M"),
               _nodes[i].x_coord(), _nodes[i].y_coord(),
               _nodes[i].feeder_id(),
               _nodes[i].parent_id(), _nodes[i].parent_arc_id(), 
               _nodes[i].children_count(),
               _nodes[i].pl_kw(), _nodes[i].ql_kvar(),
               _nodes[i].v_k_1(), _nodes[i].v_k());
        printf("  arcs(%d): ",_nodes[i].fans().size());
        for ( long k = 0; k < _nodes[i].fans().size(); k ++) {
            printf(" %d:%d->%d", _nodes[i].fans()[k],
            _arcs[ _nodes[i].fans()[k]].src_id(),
            _arcs[ _nodes[i].fans()[k]].des_id());
        }
        printf("  prob_intervals(%d): {",_nodes[i].prob_intervals().size());
        for ( long k = 0; k < _nodes[i].prob_intervals().size(); k ++) {
            printf(" %d", _nodes[i].prob_intervals()[k].to_i());
        } printf(" }");
    }
    printf("\nArcs: %d", _arcs.size());
    for ( long j = 0; j < _arcs_count; j ++) {
        printf("\n%d %d->%d type:%s", j,
            _arcs[j].src_id(), _arcs[j].des_id(),
            (_arcs[j].type() == TIE) ? "TIE" : "S_NS");
        printf("  r:%.5f x:%.5f g:%.5f ", _arcs[j].r(), _arcs[j].x(), _arcs[j].g());
    }
    printf("\n");
}

void NERDS_GRAPH::print_voltages()
{
    //printf("\nNode voltages: \n");
    //for ( long i = 0; i < _nodes_count; i ++) { printf("%d ", i); }
    //printf("\nPk:\n");
    //for ( long i = 0; i < _nodes_count; i ++) {
    //  printf("%d:%.4f:%.4f   ", i, _nodes[i].P(), _nodes[i].Q()); 
    //}
    printf("\nVk:\n");
    for ( long i = 0; i < _nodes_count; i ++) {
        if ( _nodes[i].type() != SUPER_FEEDER)
            //if ( _nodes[i].feeder_id() == 5)
            printf("%d: %.8f \n", i, _nodes[i].v_k());
            //printf("%.8f \n", _nodes[i].v_k());   
    }
    //printf("\nVk**2:\n");
    //for ( long i = 0; i < _nodes_count; i ++) {
    //  printf("%d:%.8f ", i, _nodes[i].v_k()*_nodes[i].v_k()); 
    //}
}

void NERDS_GRAPH::print_feeders_depth()
{
    // for user's entertainment;
    printf("Feeder depths: ");
    for ( long j = 0; j < _feeders_count; j ++) {
        printf(" %d", _feeders_depth[ j]);
    }
    printf("\n");
}

void NERDS_GRAPH::print_losses()
{
    // in each arc; to be called after running power flow estimation;
    printf("\nLosses in arcs:\n");
    double Pi = 0., Qi = 0., Vi = 1., ri = 0., loss_i = 0.;
    int src_id = 0, des_id = 0;
    for ( long i = 0; i < _arcs_count; i ++) {
        if ( _arcs[i].type() != TIE) {
            _arcs[i].get_P_Q( Pi, Qi); // P, Q thru this branch;
            ri = _arcs[i].r();
            src_id = _arcs[ i].src_id();
            des_id = _arcs[ i].des_id();            
            if ( _nodes[ des_id].parent_id() == src_id) {
                Vi = _nodes[ src_id].v_k();
            } else {
                Vi = _nodes[ des_id].v_k();
            }
            loss_i = ri * ( (Pi*Pi + Qi*Qi)/(Vi*Vi));
            printf("%d:%.4f ", i, loss_i);
        }
    }
}

void NERDS_GRAPH::print_left_right_loop(int tie_i)
{
    printf("\n\nLoop for tie switch %d:", tie_i);
    printf("\nLeft arcs: ");
    int left_count = _sketch_left_arc_ids.size();
    for ( long i = 0; i < left_count; i ++) {
        printf(" %d", _sketch_left_arc_ids[ i]);
    }
    printf("\nRight arcs: ");
    int right_count = _sketch_right_arc_ids.size();
    for ( long i = 0; i < right_count; i ++) {
        printf(" %d", _sketch_right_arc_ids[ i]);
    }
}

void NERDS_GRAPH::print_sketch_exchange_losses()
{
    int ex_count = _sketch_exchange_losses.size();
    if ( ex_count > 0) {
        printf("\nExchange losses:");
        for ( long i = 0; i < ex_count; i ++) {
            printf("\n %d -> %d  %.16f", _sketch_exchange_losses[ i].tie_i(),
                   _sketch_exchange_losses[ i].sec_i(), _sketch_exchange_losses[ i].loss());
        }
    } else {
        printf("\n---> Exchange losses array is empty!");
    }
}


////////////////////////////////////////////////////////////////////////////////
//
// NERDS_GRAPH: New power-flow solution: idea is based on random-walks;
//
////////////////////////////////////////////////////////////////////////////////

// control what is the minimum prob interval used to artificialy correct
// nodes with adjacent arcs with r,x normal and also with r,x very high that
// would translate into zero prob interval, whcih would bloc the random
// walk to ascend toward root;
#define MINIMUM_PROB_INTERVAL 5
// controls the accuracy runtime tradeoff;
#define RANDOM_WALKS_COUNT_LIMIT 1 // default 1;

bool NERDS_GRAPH::run_PowerFlow_solution_rw_version1()
{
    // this implements PowerFlow solution based on random-walks;
    // Version 1:
    // -- select first all leaf nodes and compute their voltages
    // -- convert them to HOMEs
    // -- process the remaining nodes;
    // issue: leaf nodes are far away from root -> longer walks ->runtime
    // it's intended to replace DistFlow method of Baran and to be faster;
    // for details see Qian and Sapatnekar paper;

    bool debug_runtime = true; // runtime profiling;
    timeval time_val_1, time_val_2; // sec and microsec;
    if ( debug_runtime) { gettimeofday( &time_val_1, 0); }
    printf("\nRunning PowerFlow random walks solution...");


    // (0) clean-up;
    clear_nodes();


    // (1)
    // (a) one time initial dry calculation of P,Q's from leaves towards root;
    for ( long f_i = 0; f_i < _feeders_count; f_i++) {
        NERDS_NODE *feeder_node = &_nodes[ f_i];
        initial_backward_P_Q_aggregation_feeder_tree( feeder_node);
    }
    // (b) one time initial calculation of probability intervals that
    // model transition probabilities; also set a first set of HOME's;
    // Note: this does mem allocations and deallocation for vector; should
    // address this;
    for ( long f_i = 0; f_i < _feeders_count; f_i++) {
        NERDS_NODE *feeder_node = &_nodes[ f_i];
        compute_probability_intervals( feeder_node); // mem alloc issue?
        feeder_node->set_lodging_type( HOME_NODE);
        feeder_node->set_lodging_cost( 1.0);
        feeder_node->set_v_k( 1.0); // 1 V p.u.;
    }
    //print_graph(); // debug;
    

    // (2) for every feeder tree select first all leaf nodes and compute their 
    // voltages by the random-walks method; for each of the nodes perform
    // a bounded number of random walks;
    vector<int> nodes_for_rw;
    for ( long f_i = 0; f_i < _feeders_count; f_i++) {

        // (a) get the leaf nodes;
        nodes_for_rw.clear(); // mem alloc issue?
        NERDS_NODE *feeder_node = &_nodes[ f_i];
        get_leaf_nodes( feeder_node, nodes_for_rw);

        // (b) now we should have all leaf nodes in nodes_for_rw; process them
        // and transform them into HOME's to speed up the precedure;
        // estimate voltage of every leaf node - average of a num of random walks;
        long count = nodes_for_rw.size();
        for ( long j = 0; j < count; j ++) {
            long this_node_id = nodes_for_rw[j];
            NERDS_NODE *this_node = &_nodes[ this_node_id];
            estimate_voltage_by_random_walks( this_node_id); // RANDOM WALKS;
            // printf(" estimated voltage: %f", this_node->v_k());
            // once a leaf node is processed convert it to a HOME;
            this_node->set_lodging_type( HOME_NODE);
        }
    }


    // (3) now take the ramaining nodes of the feeder tree - possibly in
    // random order - and process them;
    // Note: this should be done only for remaining nodes of selected feeders
    // in a localized manner throughout the power system; should be included 
    // in the for loop above;
    nodes_for_rw.clear();
    for ( long i = _feeders_count; i < _nodes_count; i ++) {
        NERDS_NODE *this_node = &_nodes[ i];
        if ( this_node->type() != SUPER_FEEDER && this_node->lodging_type() != HOME_NODE) {
            estimate_voltage_by_random_walks( this_node->id()); // RANDOM WALKS;
            this_node->set_lodging_type( HOME_NODE);
        }
    }


    // (4) compute Pi, Qi for all feeder trees; at this time all voltages
    // are computed using random walks; this step is just on simple forward sweep;
    for ( long f_i = 0; f_i < _feeders_count; f_i++) {
        NERDS_NODE *feeder_node = &_nodes[ f_i];
        forward_sweep_rw( feeder_node); // compute Pi, Qi;
    }


    // debug;
    if ( debug_runtime) {
        gettimeofday( &time_val_2, 0);
        _pronerds->add_to_runtime_power_flow( 
            time_val_2.tv_sec - time_val_1.tv_sec + 
            double(time_val_2.tv_usec - time_val_1.tv_usec) / 1000000.0);
        //printf ("\nRuntime walltime of PowerFlow random walks: %2.8f sec",
        //  time_val_2.tv_sec - time_val_1.tv_sec + 
        //  double(time_val_2.tv_usec - time_val_1.tv_usec) / 1000000.0);
    }
    //print_voltages();
    return true;
}

bool NERDS_GRAPH::run_PowerFlow_solution_rw_version2()
{
    // Version 2:
    // -- forward sweep and process all nodes
    // -- convert them to HOMEs as we go
    // this should make on average walks be shorter?

    bool debug_runtime = true; // runtime profiling;
    timeval time_val_1, time_val_2; // sec and microsec;
    if ( debug_runtime) { gettimeofday( &time_val_1, 0); }
    printf("\nRunning PowerFlow random walks solution...");

    // (0) clean-up;
    clear_nodes();

    // (1)
    // (a) one time initial dry calculation of P,Q's from leaves towards root;
    for ( long f_i = 0; f_i < _feeders_count; f_i++) {
        NERDS_NODE *feeder_node = &_nodes[ f_i];
        initial_backward_P_Q_aggregation_feeder_tree( feeder_node);
    }
    // (b) one time initial calculation of probability intervals that
    // model transition probabilities; also set a first set of HOME's;
    // Note: this does mem allocations and deallocation for vector; should
    // address this;
    for ( long f_i = 0; f_i < _feeders_count; f_i++) {
        NERDS_NODE *feeder_node = &_nodes[ f_i];
        compute_probability_intervals( feeder_node); // mem alloc issue?
        feeder_node->set_lodging_type( HOME_NODE);
        feeder_node->set_lodging_cost( 1.0);
        feeder_node->set_v_k( 1.0); // 1 V p.u.;
    }
    //print_graph(); // debug;

    // (2) just walk thru nodes and process; by construction of testcases nodes
    // are sort of consecutive anyway, so we'll advance from feeder nodes toward
    // leaves in this way; approximately;
    for ( long i = _feeders_count; i < _nodes_count; i ++) {
        NERDS_NODE *this_node = &_nodes[ i];
        if ( this_node->type() != SUPER_FEEDER && this_node->lodging_type() != HOME_NODE) {
            estimate_voltage_by_random_walks( this_node->id()); // RANDOM WALKS;
            this_node->set_lodging_type( HOME_NODE);
        }
    }   

    // (3) compute Pi, Qi for all feeder trees; at this time all voltages
    // are computed using random walks; this step is just on simple forward sweep;
    for ( long f_i = 0; f_i < _feeders_count; f_i++) {
        NERDS_NODE *feeder_node = &_nodes[ f_i];
        forward_sweep_rw( feeder_node); // compute Pi, Qi;
    }

    // debug;
    if ( debug_runtime) {
        gettimeofday( &time_val_2, 0);
        _pronerds->add_to_runtime_power_flow( 
            time_val_2.tv_sec - time_val_1.tv_sec + 
            double(time_val_2.tv_usec - time_val_1.tv_usec) / 1000000.0);
        //printf ("\nRuntime walltime of PowerFlow random walks: %2.8f sec",
        //  time_val_2.tv_sec - time_val_1.tv_sec + 
        //  double(time_val_2.tv_usec - time_val_1.tv_usec) / 1000000.0);
    }
    //print_voltages();
    return true;
}

bool NERDS_GRAPH::run_PowerFlow_solution_rw_update_subset_of_nodes(
    vector<int> &nodes_for_rw_update)
{
    // developed starting from the above run_PowerFlow_solution_rw_version1();
    // this implements *partial* PowerFlow solution based on random-walks;
    // it takes only a list of nodes and does the voltage calculations only
    // for these, keeping the other nodes at their voltages from the previous
    // run of the PowerFlow solutions; here we also set as HOME only the root;
    // we should also experiment with keeping some of the other nodes (for 
    // which we do not recalculate voltages) also as HOME?
    // simplified steps:
    // -- process all nodes from the provided list
    // -- convert them to HOMEs as they are processed

    bool debug_runtime = true; // runtime profiling;
    timeval time_val_1, time_val_2; // sec and microsec;
    if ( debug_runtime) { gettimeofday( &time_val_1, 0); }
    printf("\nRunning *update* PowerFlow random walks solution...");


    // (0) clean-up; this time we only make all nodes as motels and only the 
    // nodes from subset have their prob intervals cleared; everything else
    // Pk,Qk,Vk for all nodes and P,Q for all arcs are kept;
    clear_nodes_partially( nodes_for_rw_update);
    
    
    // (1)
    // (a) one time initial dry calculation of P,Q's from leaves towards root;
    // only for affected feeders;
    long affected_feeders_count = _sketch_affected_feeders.size();
    for ( long f_i = 0; f_i < affected_feeders_count; f_i++) {
        NERDS_NODE *feeder_node = &_nodes[ _sketch_affected_feeders[f_i] ];
        initial_backward_P_Q_aggregation_feeder_tree( feeder_node);
    }
    // (b) one time initial calculation of probability intervals that
    // model transition probabilities; also set a first set of HOME's;
    // Note: this does mem allocations and deallocation for vector; should
    // address this; also every node processed here should have had its
    // prob intervals cleared; here they are allocated and computed again;
    long count = nodes_for_rw_update.size();
    for ( long j = 0; j < count; j ++) {
        NERDS_NODE *this_node = &_nodes[ nodes_for_rw_update[j] ];
        compute_probability_intervals_one_node( this_node);
    }
    for ( long f_i = 0; f_i < _feeders_count; f_i++) {
        NERDS_NODE *feeder_node = &_nodes[ f_i];
        feeder_node->set_lodging_type( HOME_NODE);
        feeder_node->set_lodging_cost( 1.0);
        feeder_node->set_v_k( 1.0); // 1 V p.u.;
    }
    // (c) consider as home nodes also all nodes at 1/2 max depth in the
    // affected feeder-trees; this way we shothen the hitting times of
    // random walk;
    //convert_to_home_nodes_at_half_depth(); // if use this uncomment stuff inside clear_nodes_partially();
    convert_to_home_all_nodes_but_2nd_order_neighbours( nodes_for_rw_update);
    //print_graph(); // debug;
    //exit(1);
    

    // (2) get each node from the passed list of nodes to be updated, and
    // process to compute their voltages by the random-walks method; for each 
    // of the nodes perform a bounded number of random walks;
    for ( long j = 0; j < count; j ++) {
        long this_node_id = nodes_for_rw_update[j];
        NERDS_NODE *this_node = &_nodes[ this_node_id];
        estimate_voltage_by_random_walks( this_node_id); // RANDOM WALKS;
        // once a node is processed convert it to a HOME;
        this_node->set_lodging_type( HOME_NODE);
    }


    // (3) compute Pi, Qi for all "affected" feeder trees; only feeder trees
    // which had nodes to be updated are "affected"; at this time all voltages 
    // of all nodes are known: some were from previous runs of the PowerFlow 
    // solutions, and some were just updated; this step is just on simple 
    // forward sweep;
    for ( long f_i = 0; f_i < affected_feeders_count; f_i++) {
        NERDS_NODE *feeder_node = &_nodes[ _sketch_affected_feeders[f_i] ];
        forward_sweep_rw( feeder_node); // compute Pi, Qi;
    }
    

    // debug;
    if ( debug_runtime) {
        gettimeofday( &time_val_2, 0);
        _pronerds->add_to_runtime_power_flow( 
            time_val_2.tv_sec - time_val_1.tv_sec + 
            double(time_val_2.tv_usec - time_val_1.tv_usec) / 1000000.0);
        //printf ("\nRuntime walltime of *update* PowerFlow random walks: %2.8f sec",
        //  time_val_2.tv_sec - time_val_1.tv_sec + 
        //  double(time_val_2.tv_usec - time_val_1.tv_usec) / 1000000.0);
    }
    //print_voltages();
    return true;
}

void NERDS_GRAPH::convert_to_home_nodes_at_half_depth()
{
    // this is to speed up the partial RW based PowerFlow solution;
    deque<long> f_queue; // forward queue;

    long affected_feeders_count = _sketch_affected_feeders.size();
    for ( long f_i = 0; f_i < affected_feeders_count; f_i++) {
        NERDS_NODE *feeder_node = &_nodes[ _sketch_affected_feeders[f_i] ];
        long feeder_half_depth = _feeders_depth[ _sketch_affected_feeders[f_i] ];
        feeder_half_depth = ( feeder_half_depth >> 1);

        NERDS_NODE *upstream_node = feeder_node;
        NERDS_NODE *downstream_node = 0;
        long upstream_id = upstream_node->id();
        long downstream_id = 0;

        // walk feeder_tree from root towards leaves to find nodes at 1/2 depth;
        f_queue.push_back( upstream_id);
        while ( !f_queue.empty()) {
            // get node from queue and process;
            upstream_id = f_queue.back(); // get the last one: lifo;
            f_queue.pop_back(); // remove it;
            upstream_node = &_nodes[ upstream_id];
            // look at all neighbors (adjanct nodes) and walk them thru
            // except the parent from where we came to this one;
            long neighbors_count = upstream_node->fans().size();
            long parent_id = upstream_node->parent_id();
            for ( long k = 0; k < neighbors_count; k++) {
                NERDS_ARC *arc = &_arcs[ upstream_node->fans( k)];
                if ( arc->type() != TIE) {
                    downstream_id = ( arc->des_id() != upstream_id) ? 
                        arc->des_id() : arc->src_id();
                    if ( downstream_id != parent_id) { // skip parent;
                        if ( upstream_node->depth() == feeder_half_depth) {
                            upstream_node->set_lodging_type( HOME_NODE);
                            //printf(" ---> converted to home %d", upstream_node->id());
                        } else {
                            f_queue.push_back( downstream_id);
                        }
                    }
                }
            }
        }
    }
}

void NERDS_GRAPH::convert_to_home_all_nodes_but_2nd_order_neighbours(
    vector<int> &nodes_for_rw_update)
{
    // this is to speed up the partial RW based PowerFlow solution;
    // we actually do not convert to homes, but we go and convert to models
    // nodes to be updated and their neighbors; this is because I commented
    // out stuff inside clear_nodes_partially(); do not forget that!
    deque<long> b_queue; // forward queue;

    for ( long i = 0; i < _nodes_count; i ++) {
        NERDS_NODE *this_node = &_nodes[ i];
        this_node->set_lodging_type( HOME_NODE);
    }

    long count = nodes_for_rw_update.size();
    for ( long j = 0; j < count; j ++) {
        long this_node_id = nodes_for_rw_update[j];
        b_queue.push_back( this_node_id);
        long this_depth = _nodes[ this_node_id].depth();

        NERDS_NODE *node = 0;
        NERDS_NODE *n_node = 0; // neighbour;
        long id = 0;
        // for now, for simplicity I look only at neighbors going up (backward);
        // however we could have neighbors the other children of this node as well
        // as possible children of this node (this node is the initial node from the
        // list of nodes to be updated);
        while ( !b_queue.empty()) {
            // get node;
            id = b_queue.back();
            b_queue.pop_back(); // remove;
            node = &_nodes[ id];
            // convert to motel;
            node->set_lodging_type( MOTEL_NODE);
            // get its parent, add to queue if within neighborhood;
            n_node = &_nodes[ node->parent_id()];
            if ( n_node->type() != FEEDER) {
                if ( abs(n_node->depth() - this_depth) <= 2)
                    b_queue.push_back( n_node->id());
            }
        }
    }
}

void NERDS_GRAPH::get_leaf_nodes( NERDS_NODE *feeder_node, vector<int> &nodes_for_rw)
{
    // forward sweep to get all leaves;
    NERDS_NODE *upstream_node = feeder_node;
    NERDS_NODE *downstream_node = 0;
    long upstream_id = upstream_node->id();
    long downstream_id = 0;
    deque<long> f_queue; // forward queue;
    bool is_leaf_node = false;

    // walk feeder_tree from root towards leaves to find the leaves;
    f_queue.push_back( upstream_id);
    while ( !f_queue.empty()) {

        // get node from queue and process;
        upstream_id = f_queue.back(); // get the last one: lifo;
        f_queue.pop_back(); // remove it;
        upstream_node = &_nodes[ upstream_id];

        is_leaf_node = true;
        // look at all neighbors (adjanct nodes) and walk them thru
        // except the parent from where we came to this one;
        long neighbors_count = upstream_node->fans().size();
        long parent_id = upstream_node->parent_id();
        for ( long k = 0; k < neighbors_count; k++) {
            NERDS_ARC *arc = &_arcs[ upstream_node->fans( k)];
            if ( arc->type() != TIE) {
                downstream_id = ( arc->des_id() != upstream_id) ? 
                    arc->des_id() : arc->src_id();
                if ( downstream_id != parent_id) { // skip parent;
                    f_queue.push_back( downstream_id);
                    is_leaf_node = false; // has at least a child;
                }
            }
        }
        // found leaf?
        if ( is_leaf_node) {
            nodes_for_rw.push_back( upstream_id);
        }
    }
}

double NERDS_GRAPH::estimate_voltage_by_random_walks( long node_id)
{
    double total_v = 0.0;
    long rw_counter = 0;
    long steps_count = 0;
    long total_steps_count = 0;
    
    NERDS_NODE *node = &_nodes[ node_id];
    if ( node->lodging_type() == HOME_NODE)
        { return node->v_k(); } // do nothing, it's already a HOME;

    double voltage = 0.0;
    // (1) do a bounded number of random walks starting from this node;
    // if a walk terminates at a HOME aggregate the result; take the average 
    // at the end;
    for ( long rw_i = 0; rw_i < RANDOM_WALKS_COUNT_LIMIT; rw_i++) {

        voltage = 0.0; // reset;
        if ( perform_random_walk_from_node( node, voltage, steps_count)) {
            // accumulate it only if this random walk terminated at a HOME;
            total_v += voltage;
            total_steps_count += steps_count;
            rw_counter ++;
        }
    }

    // (2) average;
    total_v /= max(1, rw_counter);
    total_steps_count /= max(1, rw_counter);
    //printf("\n>>> Average number of steps per RW: %d  Number of RW's: %d  Termination steps-limit: %d",
    //     total_steps_count, rw_counter, 20 * ( _nodes_count / _feeders_count));
    //printf("\n>>> Voltage of node at the end of RW's: %.8f", total_v);
    // record it too here internally;
    node->set_v_k( total_v);

    return total_v;
}

bool NERDS_GRAPH::perform_random_walk_from_node( NERDS_NODE *node, 
    double &v, long &steps_count)
{
    // return the node to which it's decided to randomly walk, as a walk-step
    // from "node" to one of its fans as dicttaed by transistion probabilities;

    // if a walk takes more than this limit, it's terminated;
    const long STEPS_COUNT_LIMIT = 20 * ( _nodes_count / _feeders_count);

    bool terminated_at_home = true;
    double total_price = 0.0;
    int rand_int = 0;
    int next_arc_id = 0, next_node_id = 0;

    if ( node->lodging_type() == HOME_NODE)
        { return node->v_k(); } // do nothing, it's already a HOME;
    //printf("\nRandom walk: {%d} ", node->id());
    
    // (1) get the reward or penalty for lodging at this node;
    // done because if this is a HOME then we'll return without
    // performing any walk at all;
    total_price = node->lodging_cost();


    // (2) jump from one node to another till we hit a HOME;
    long counter = 0;
    NERDS_NODE *this_node = node;
    while ( this_node->lodging_type() != HOME_NODE && 
        counter < STEPS_COUNT_LIMIT) {
        //assert( this_node->lodging_cost() <= 0); // debug;
        // (a) get arc id to walk;
        rand_int = my_rand_int( 100);
        next_arc_id = arc_id_for_walk( this_node, rand_int);
        //assert( next_arc_id >= 0);
        // (b) get the node id to actual node we jump to;
        NERDS_ARC *arc = &_arcs[ next_arc_id];
        next_node_id = ( arc->des_id() != this_node->id()) ?
            arc->des_id() : arc->src_id();
        this_node = &_nodes[ next_node_id];
        // (c) accumulate penalty or final reward for lodging at this node;
        total_price += (this_node->lodging_type() != HOME_NODE) ? 
            this_node->lodging_cost() : this_node->v_k();
        //printf("^%d", next_node_id);
        counter ++;
    }

    if ( this_node->lodging_type() != HOME_NODE) {
        //printf(" {RW not finished!}");
        // if the random walk wandered around for too long an was terminated,
        // still "assume" the it reached a home and add the "reward" of 1;
        // we introduce errors this way but it's better than leaving the voltage 
        // to something negative very small;
        total_price += 1;
        terminated_at_home = true;
    }
    // upload the final result, which is valid or "artificially valid" depending
    // if the walk terminated at one of the HOMEs or not; see above comments;
    v = total_price;
    steps_count = counter;
    //printf(" {%d %f}", node->id(), v);

    return terminated_at_home;
}

inline int NERDS_GRAPH::arc_id_for_walk( NERDS_NODE *node, int rand_int) 
{
    // get the entry PAIR_TWO from "_prob_intervals" that is the
    // first with prob > rand_int;

    int arc_id = -1;
    long count = node->prob_intervals().size();
    for ( long k = 0; k < count; k ++) {
        if ( rand_int <= node->prob_intervals( k).to_i() ) {
            // node->prob_intervals( k).from_i() stores the index "k" 
            // of element in _fans that represents the arc to be walked;
            arc_id = node->fans( node->prob_intervals( k).from_i() );
            break;
        }
    }
    //assert( arc_id >= 0);
    return arc_id;
}

void NERDS_GRAPH::compute_probability_intervals( NERDS_NODE *feeder_node)
{
    // should be called only once per PowerFlow solution calculation!
    // this performs a forward walk on this feeder tree and computes and
    // popultaes the probability and probbability intervals for each node;
    // description of meaning of probability intervals; assume some node n4:
    //      |
    //      n3
    //      |
    // n2...n4--n7            n4...n2 is a tie switch
    //     / \
    //    n5  n6
    //  
    // that has "fan" list {n3 n2 n5 n6 n7} via arcs with conductances "g", as
    // stored inside the tree-graph already; node n4 will have "_prob_intervals" 
    // {(0,20) (2,40) (3,50) (4,100)} which means
    // that: prob(n4->n3)=1/5, prob(n4->n5)=1/5, prob(n4->n6)=1/10, and
    // prob(n4->n7)=1/2; when we'll "random-walk", we'll generate a rand number
    // between 0-100 and walk arc n4->n6 is that number is between 40-50;
    // note that "_prob_intervals" does not contain an entry for index "1"
    // of _fans because that is an arc for a tie switch; we use PAIR_TWO
    // elements (index in fans, prob) to directly retreive index for fans;
    // 
    // Note: there are some pathological examples (such as 202 bus testcase)
    // that has 4--->7--->8, with 7->8 having r=0.0000001 which would make
    // the probability range for going from 7 to 4 (hence toward root) zero;
    // to get arround these I artificially inflate it to MINIMUM_PROB_INTERVAL %;
    // this way a random walk has a better chance to go back to root and terminate;
    //
    // here we also compute current loads/sinks at each node; under some
    // simplifying assumptions;

    NERDS_NODE *upstream_node = feeder_node;
    NERDS_NODE *downstream_node = 0;
    long upstream_id = upstream_node->id();
    long downstream_id = 0;
    deque<long> f_queue; // forward queue;
    double g_total = 0.0;
    long prob_total = 0;
    long delta_prob = 0;
    double current_I = 0.0;
    double m_x = 0.0;
    //printf("\n m_x's:");

    // (1) simple forward walk of the tree;
    f_queue.push_back( upstream_id);
    while ( !f_queue.empty()) {

        // (a) get node from queue and process;
        upstream_id = f_queue.back(); // get the last one: lifo;
        f_queue.pop_back(); // remove it;
        upstream_node = &_nodes[ upstream_id];

        // (b) 
        // () compute the total resistance of all arcs of this node;
        long neighbors_count = upstream_node->fans().size();
        g_total = 0.0;
        for ( long k = 0; k < neighbors_count; k++) {
            NERDS_ARC *arc = &_arcs[ upstream_node->fans( k)];
            if ( arc->type() != TIE) { // good branch;
                g_total += arc->g(); // dc, ac case;
            }
        }
        // () compute the current I and then the penalty/reward for lodging;
        // we assume sink/load-current I = P + j Q / V, where V is the voltage
        // at this node, which ideally is 1; after simplification:
        // I = mag(I) = sqrt(P^2 + Q^2); neglect phase for now;
        current_I = sqrt( 
            upstream_node->pl_kw() * upstream_node->pl_kw() + 
            upstream_node->ql_kvar() * upstream_node->ql_kvar());
        m_x = current_I / g_total;
        upstream_node->set_lodging_cost( -m_x);
        //printf(" %.16f", -m_x);

        // (c) go one more thru all fan arcs of this node and compute
        // "transition" probabilities;
        prob_total = 0;
        long parent_id = upstream_node->parent_id();
        for ( long k = 0; k < neighbors_count; k++) {
            NERDS_ARC *arc = &_arcs[ upstream_node->fans( k)];
            if ( arc->type() != TIE) {
                delta_prob = long(100 * ( arc->g() / g_total)); // dc, ac case;
                // correct any possible degenerate situations that would
                // result into zero prob interval;
                delta_prob = max( delta_prob, MINIMUM_PROB_INTERVAL);
                prob_total += delta_prob;
                prob_total = min( prob_total, 100);
                upstream_node->add_prob_intervals( PAIR_TWO(k, prob_total));
                // keep track and add children to processing queue;
                downstream_id = ( arc->des_id() != upstream_id) ? 
                    arc->des_id() : arc->src_id();
                if ( downstream_id != parent_id) { // skip parent;
                    downstream_node = &_nodes[ downstream_id];
                    f_queue.push_back( downstream_id);
                }
            }
        }
        // correct the last value inside _prob_intervals to be 100, as
        // it my be less due to roundoff errors;
        upstream_node->set_prob_intervals(
            upstream_node->prob_intervals().size() - 1, long(100));
    }
}

void NERDS_GRAPH::compute_probability_intervals_one_node( NERDS_NODE *node)
{
    // self explanatory; used when doing only partial rw PowerFlow;
    double g_total = 0.0;
    long prob_total = 0;
    long delta_prob = 0;
    double current_I = 0.0;
    double m_x = 0.0;
    //printf("\n m_x's:");

    // () compute the total resistance of all arcs of this node;
    long neighbors_count = node->fans().size();
    g_total = 0.0;
    for ( long k = 0; k < neighbors_count; k++) {
        NERDS_ARC *arc = &_arcs[ node->fans( k)];
        if ( arc->type() != TIE) { // good branch;
            g_total += arc->g(); // dc, ac case;
        }
    }

    // () compute the current I and then the penalty/reward for lodging;
    // we assume sink/load-current I = P + j Q / V, where V is the voltage
    // at this node, which ideally is 1; after simplification:
    // I = mag(I) = sqrt(P^2 + Q^2); neglect phase for now;
    current_I = sqrt( 
                     node->pl_kw() * node->pl_kw() + 
                     node->ql_kvar() * node->ql_kvar());
    m_x = current_I / g_total;
    node->set_lodging_cost( -m_x);
    //printf(" %.16f", -m_x);

    // () go one more thru all fan arcs of this node and compute
    // "transition" probabilities;
    prob_total = 0;
    for ( long k = 0; k < neighbors_count; k++) {
        NERDS_ARC *arc = &_arcs[ node->fans( k)];
        if ( arc->type() != TIE) {
            delta_prob = long(100 * ( arc->g() / g_total)); // dc, ac case;
            delta_prob = max( delta_prob, MINIMUM_PROB_INTERVAL);
            prob_total += delta_prob;
            prob_total = min( prob_total, 100);
            node->add_prob_intervals( PAIR_TWO(k,prob_total));
        }
    }
    // correct the last value inside _prob_intervals to be 100, as
    // it my be less due to roundoff errors;
    node->set_prob_intervals( node->prob_intervals().size() - 1, long(100));
}

void NERDS_GRAPH::forward_sweep_rw( NERDS_NODE *feeder_node)
{
    // this is needed after the voltages are computed using random walks;
    // this is just one forward sweep to compute Pi's and Qi's;
    NERDS_NODE *upstream_node = feeder_node;
    NERDS_NODE *downstream_node = 0;
    long upstream_id = upstream_node->id();
    long downstream_id = 0;
    deque<long> f_queue; // forward queue;

    // (1) simple forward walk of the tree;
    f_queue.push_back( upstream_id);
    while ( !f_queue.empty()) {

        // (a) get node from queue and process;
        upstream_id = f_queue.back(); // get the last one: lifo;
        f_queue.pop_back(); // remove it;
        upstream_node = &_nodes[ upstream_id];

        // (b) take every child, process, and then add it to queue;
        long neighbors_count = upstream_node->fans().size();
        bool is_fork_node = ( neighbors_count > 2);
        long parent_id = upstream_node->parent_id();
        for ( long k = 0; k < neighbors_count; k++) {
            NERDS_ARC *arc = &_arcs[ upstream_node->fans( k)];
            if ( arc->type() != TIE) { // good branch;
                downstream_id = ( arc->des_id() != upstream_id) ? 
                    arc->des_id() : arc->src_id();
                if ( downstream_id != parent_id) { // skip parent;
                    downstream_node = &_nodes[ downstream_id];
                    // calculate Pi Qi using eq. 1 from Baran paper; however,
                    // note that here Vi's are computed already by random walks;
                    forward_branch_equations_rw( upstream_node, downstream_node,
                        is_fork_node);
                    f_queue.push_back( downstream_id);
                }
            }
        }
    }
}

void NERDS_GRAPH::forward_branch_equations_rw( NERDS_NODE *u_node,
    NERDS_NODE *d_node, bool is_fork_node)
{
    // used during the just one forward sweep to compute Pi's and Qi's,
    // after Vi's were computed using random walks;

    // (1) simply use eq 1 from Baran paper;
    NERDS_ARC *arc_i = &_arcs[ d_node->parent_arc_id()]; // arc i;
    double r_i = arc_i->r(); // ri;
    double x_i = arc_i->x(); // xi;
    // Note: tricky part: at fork nodes we take a fraction of P and Q
    // calculated until now during this forward sweep; the fraction
    // is computed using prev-iter P,Q stored on the arc of each child;
    // Pk_1 is what u_node saw cumulatively downstream during prev iter;
    double Pi = 0.0;
    double Qi = 0.0;
    Pi = ( is_fork_node) ?
        ( u_node->P() * ( arc_i->P() / u_node->Pk_1())) : ( u_node->P());
    Qi = ( is_fork_node) ?
        ( u_node->Q() * ( arc_i->Q() / u_node->Qk_1())) : ( u_node->Q());
    
    double Vi_sq = ( u_node->v_k() * u_node->v_k()); // computed using random walks;
    double I_sq = 0.0;
    double Pi_1 = 0.0;
    double Qi_1 = 0.0;
    I_sq = ( Pi*Pi + Qi*Qi) / ( Vi_sq);
    Pi_1 = Pi - d_node->pl_kw() - r_i * I_sq; // Pi+1;
    Qi_1 = Qi - d_node->ql_kvar() - x_i * I_sq; // Qi+1;
    // before saving/overwriting, first copy old values into k_1 storage;
    // Pi_1,Qi_1 are the P,Q "flowing" downstream from this d_node; at
    // fork nodes this flowing is split between children; Note: in this case
    // (forward sweep) we do not store P,Q on arc because P,Q has meaning of
    // accumulated P,Q supplying all children of a node;
    d_node->copy_P_Q_to_Pk_1_Qk_1();
    d_node->set_P( Pi_1);
    d_node->set_Q( Qi_1);
}
