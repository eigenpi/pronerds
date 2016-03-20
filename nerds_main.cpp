#include "nerds_gui.h"
#include "nerds.h"
#include "nerds_mcmf.h"
#include <assert.h>
#include <time.h>

using namespace std;

////////////////////////////////////////////////////////////////////////////////
//
// launching point;
//
////////////////////////////////////////////////////////////////////////////////

// next two are for debug only;
int test_mcmf_edmonds();
int test_mcmf_cs2();

int main( int argc, char *argv[])
{
    // welcome;
    char welcome[] =
        "\n--------------------------------------------------\n"
        "PRO-NERDS\n"
        "cristinel.ababei@ndsu.edu, Compiled "__DATE__" \n"
        "--------------------------------------------------\n";
    printf("%s", welcome);

    timeval time_val_1, time_val_2; // sec and microsec;
    clock_t start_clock, end_clock;
    clock_t diff_clock;
    char msg[BUFFER_SIZE];
 

    // (1) create the network;
    PRONERDS pronerds;

    
    // (2) create the network graph from the network_file.pos;
    pronerds.parse_command_arguments( argc, argv);
    pronerds.create_network_graph( argc, argv);
    pronerds.print_initial_stats( argc, argv); // entertain user;


    // (3) create empty gui object and populate if required by user;
    GUI_GRAPHICS gui( &pronerds);
    pronerds.set_gui( &gui); // innitially gui is empty;
    if ( pronerds.use_gui()) { // if user asked to use the gui;
        // mark flag that we are gonna use the gui; set is_gui_usable
        // and wait_for_user_input_automode; then build;
        gui.set_graphics_state( true, 1);
        gui.build();
        gui.init_draw_coords( 32.0, 168.0);
    } else { // gui is not usable;
        gui.set_graphics_state( false, 1);
    }


    // (4) GUI stuff: entertain user before;
    if ( pronerds.use_gui()) {
        sprintf( msg, "Initial configuration %d", 0);
        pronerds.gui()->update_screen( PRIORITY_MAJOR, msg, NODES);
    }


    // (6) work;
    gettimeofday( &time_val_1, 0);
    start_clock = clock();
    assert(start_clock != (clock_t)(-1));



    // run the thing using the requested PFS method;
    // pfs_algo can be:
    //     BARAN_DIST_FLOW (default, Baran's)
    //     RANDOM_WALK_POWER_FLOW_PARTIAL
    //     RANDOM_WALK_POWER_FLOW_FULL
    // reconfig_algo can be:
    //     RECONFIG_BARAN (default, Baran's)
    //     RECONFIG_MCMF  (MCMF based, proposed)
    // mcmf_algo can be:
    //     MCMF_ALGO_EDMONDS (default, Edmonds-Karp)
    //     MCMF_ALGO_CS2     (Goldberg's)
    POWER_FLOW_METHOD pfs_algo = pronerds.pfs_algo();
    
    // (a) baran1 method: reference;
    if ( pronerds.reconfig_algo() == RECONFIG_BARAN) {
        pronerds.run_reconfiguration_method1_baran(
            false, // simplified_distflow;
            pfs_algo); // power flow solution;
    }
    // (b) new/proposed method;
    else if ( pronerds.reconfig_algo() == RECONFIG_MCMF) { // default;
        pronerds.run_reconfiguration_method_new(
            false, // simplified_distflow;
            pfs_algo); // power flow solution;
    }



    // debug only;
    //test_mcmf_edmonds();
    //test_mcmf_cs2();
    //pronerds.nerds_graph()->testing_distflow( false); // testing and debug;
    //pronerds.nerds_graph()->testing_hitting_times( false); // testing and debug;
    //pronerds.run_DistFlow_solution_baran( false); // simplified_distflow = true default;
    //pronerds.run_PowerFlow_solution_rw( 1); // version 1 or 2;
    //pronerds.nerds_graph()->print_voltages();
    

    gettimeofday( &time_val_2, 0);
    end_clock = clock();
    assert(end_clock != (clock_t)(-1));
    diff_clock = end_clock - start_clock;


    // (7) GUI stuff: entertain user after;
    if ( pronerds.use_gui()) {
        sprintf( msg, "Final configuration %d", 1);
        pronerds.gui()->update_screen( PRIORITY_MAJOR, msg, NODES);
    }
    //pronerds.nerds_graph()->print_voltages(); // debug;

    printf ("\n");
    pronerds.nerds_graph()->print_feeders_depth();
    pronerds.nerds_graph()->print_tie_switches();
    //printf ("cputime : processor time used = %2.4f sec\n", (double)diff_clock/CLOCKS_PER_SEC);
    double diff_sec_usec = time_val_2.tv_sec - time_val_1.tv_sec + 
        double(time_val_2.tv_usec - time_val_1.tv_usec) / 1000000.0;
    pronerds.add_to_runtime_all( diff_sec_usec); // record it too;
    printf ("walltime : fraction spent on power flow = %2.8f sec ( %.2f\% )\n",
            pronerds.runtime_power_flow(),
            (100 * pronerds.runtime_power_flow() / diff_sec_usec));
    printf ("walltime : elapsed (wall clock) time = %2.8f sec\n", diff_sec_usec);
    printf ("\n");

        
    // (9) do some clean-up;
    if ( gui.is_gui_usable()) {
        gui.close_graphics(); // close down X Display;
    }
}


////////////////////////////////////////////////////////////////////////////////
//
// testbenches;
//
////////////////////////////////////////////////////////////////////////////////

int test_mcmf_edmonds()
{
    // use on a small cooked example; debugging purposes only;
    int num_vertices = 18;
    int num_edges = 27;
    MCMF_EDMONDS my_mcmf_problem( num_vertices);

    int bound_nn = 2*3 + 15 + 2 ; // 2 * _feeders_count + _ties_count + 2;
    my_mcmf_problem.set_NN_and_allocate_arrays( bound_nn);

    my_mcmf_problem.set_edge(0, 1, 1, 94);
    my_mcmf_problem.set_edge(1, 2, 1, 0);
    my_mcmf_problem.set_edge(0, 3, 1, 66);
    my_mcmf_problem.set_edge(3, 2, 1, 0);
    my_mcmf_problem.set_edge(0, 4, 1, 35);
    my_mcmf_problem.set_edge(4, 2, 1, 0);
    my_mcmf_problem.set_edge(0, 5, 1, 1);
    my_mcmf_problem.set_edge(5, 2, 1, 0);
    my_mcmf_problem.set_edge(0, 6, 1, 26);
    my_mcmf_problem.set_edge(6, 2, 1, 0);

    my_mcmf_problem.set_edge(7, 8, 1, 78);
    my_mcmf_problem.set_edge(8, 2, 1, 0);
    my_mcmf_problem.set_edge(7, 9, 1, 80);
    my_mcmf_problem.set_edge(9, 2, 1, 0);

    my_mcmf_problem.set_edge(10, 11, 1, 87);
    my_mcmf_problem.set_edge(11, 12, 1, 0);
    my_mcmf_problem.set_edge(10, 13, 1, 41);
    my_mcmf_problem.set_edge(13, 12, 1, 0);
    my_mcmf_problem.set_edge(10, 14, 1, 68);
    my_mcmf_problem.set_edge(14, 12, 1, 0);
    my_mcmf_problem.set_edge(10, 15, 1, 59);
    my_mcmf_problem.set_edge(15, 12, 1, 0);

    my_mcmf_problem.set_edge(16, 0, 1, 0);
    my_mcmf_problem.set_edge(16, 7, 1, 0);
    my_mcmf_problem.set_edge(16, 10, 1, 0);

    my_mcmf_problem.set_edge(2, 17, 1, 0);
    my_mcmf_problem.set_edge(12, 17, 1, 0);


    int flow_cost = 0;
    int max_flow = my_mcmf_problem.run_edmonds( num_vertices, 16, 17, flow_cost); // S,T;
    printf("\n flow: %d", max_flow);
    printf("\n cost: %d\n", flow_cost);
    for ( int i = 0; i < num_vertices; i++) {
        for ( int j = 0; j < num_vertices; j++) {
            int this_flow = my_mcmf_problem.get_fnet(i,j);          
            if ( this_flow)
                printf("\n %d -> %d: %d cost: %d", i, j, 
                       my_mcmf_problem.get_fnet(i,j), my_mcmf_problem.get_cost(i,j));
        }
    }

    return 0;
}


int test_mcmf_cs2()
{
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

    //my_mcmf_problem.set_arc( 1, 2, 0, 4, 1);
    //my_mcmf_problem.set_arc( 1, 3, 0, 8, 5);
    //my_mcmf_problem.set_arc( 2, 3, 0, 5, 0);
    //my_mcmf_problem.set_arc( 3, 5, 0, 10, 1);
    //my_mcmf_problem.set_arc( 5, 4, 0, 8, 0);
    //my_mcmf_problem.set_arc( 5, 6, 0, 8, 9);
    //my_mcmf_problem.set_arc( 4, 2, 0, 8, 1);
    //my_mcmf_problem.set_arc( 4, 6, 0, 8, 1);
    //my_mcmf_problem.set_supply_demand_of_node( 1, 10);
    //my_mcmf_problem.set_supply_demand_of_node( 6, -10);

    // pos_873_7 iter=1;
    int num_nodes = 37;
    int num_arcs = 58;
    MCMF_CS2 my_mcmf_problem( num_nodes, num_arcs);

    my_mcmf_problem.set_arc(0, 1, 0, 1, 99);
    my_mcmf_problem.set_arc(1, 2, 0, 1, 0);

    my_mcmf_problem.set_arc(3, 4, 0, 1, 89);
    my_mcmf_problem.set_arc(4, 5, 0, 1, 0);
    my_mcmf_problem.set_arc(3, 10, 0, 1, 38);
    my_mcmf_problem.set_arc(10, 11, 0, 1, 0);
    my_mcmf_problem.set_arc(3, 14, 0, 1, 98);
    my_mcmf_problem.set_arc(14, 15, 0, 1, 0);

    my_mcmf_problem.set_arc(6, 7, 0, 1, 76);
    my_mcmf_problem.set_arc(7, 5, 0, 1, 0);
    my_mcmf_problem.set_arc(6, 25, 0, 1, 87);
    my_mcmf_problem.set_arc(25, 11, 0, 1, 0);
    my_mcmf_problem.set_arc(6, 28, 0, 1, 88);
    my_mcmf_problem.set_arc(28, 11, 0, 1, 0);
    my_mcmf_problem.set_arc(6, 30, 0, 1, 92);
    my_mcmf_problem.set_arc(30, 21, 0, 1, 0);
    my_mcmf_problem.set_arc(6, 31, 0, 1, 72);
    my_mcmf_problem.set_arc(31, 21, 0, 1, 0);
    my_mcmf_problem.set_arc(6, 32, 0, 1, 79);
    my_mcmf_problem.set_arc(32, 15, 0, 1, 0);
    my_mcmf_problem.set_arc(6, 33, 0, 1, 57);
    my_mcmf_problem.set_arc(33, 2, 0, 1, 0);

    my_mcmf_problem.set_arc(8, 9, 0, 1, 1);
    my_mcmf_problem.set_arc(9, 2, 0, 1, 0);
    my_mcmf_problem.set_arc(8, 29, 0, 1, 98);
    my_mcmf_problem.set_arc(29, 11, 0, 1, 0);

    my_mcmf_problem.set_arc(12, 13, 0, 1, 85);
    my_mcmf_problem.set_arc(13, 2, 0, 1, 0);
    my_mcmf_problem.set_arc(12, 16, 0, 1, 96);
    my_mcmf_problem.set_arc(16, 2, 0, 1, 0);
    my_mcmf_problem.set_arc(12, 17, 0, 1, 87);
    my_mcmf_problem.set_arc(17, 18, 0, 1, 0);
    my_mcmf_problem.set_arc(12, 19, 0, 1, 97);
    my_mcmf_problem.set_arc(19, 18, 0, 1, 0);
    my_mcmf_problem.set_arc(12, 20, 0, 1, 96);
    my_mcmf_problem.set_arc(20, 21, 0, 1, 0);
    my_mcmf_problem.set_arc(12, 24, 0, 1, 94);
    my_mcmf_problem.set_arc(24, 11, 0, 1, 0);
    my_mcmf_problem.set_arc(12, 27, 0, 1, 93);
    my_mcmf_problem.set_arc(27, 11, 0, 1, 0);
    my_mcmf_problem.set_arc(12, 34, 0, 1, 65);
    my_mcmf_problem.set_arc(34, 21, 0, 1, 0);

    my_mcmf_problem.set_arc(22, 23, 0, 1, 99);
    my_mcmf_problem.set_arc(23, 11, 0, 1, 0);
    my_mcmf_problem.set_arc(22, 26, 0, 1, 94);
    my_mcmf_problem.set_arc(26, 11, 0, 1, 0);

    my_mcmf_problem.set_arc(35, 0, 0, 1, 0);
    my_mcmf_problem.set_arc(35, 3, 0, 1, 0);
    my_mcmf_problem.set_arc(35, 6, 0, 1, 0);
    my_mcmf_problem.set_arc(35, 8, 0, 1, 0);
    my_mcmf_problem.set_arc(35, 12, 0, 1, 0);
    my_mcmf_problem.set_arc(35, 22, 0, 1, 0);

    my_mcmf_problem.set_arc(2, 36, 0, 1, 0);
    my_mcmf_problem.set_arc(5, 36, 0, 1, 0);
    my_mcmf_problem.set_arc(15, 36, 0, 1, 0);
    my_mcmf_problem.set_arc(11, 36, 0, 1, 0);
    my_mcmf_problem.set_arc(18, 36, 0, 1, 0);
    my_mcmf_problem.set_arc(21, 36, 0, 1, 0);

    my_mcmf_problem.set_supply_demand_of_node( 35, 5);
    my_mcmf_problem.set_supply_demand_of_node( 36, -5);

    my_mcmf_problem.run_cs2();

    return 0;
}
