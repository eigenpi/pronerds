#ifndef _NERDS_UTILS_
#define _NERDS_UTILS_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

using namespace std;

#define BUFFER_SIZE 256 // used during parsing files or gui stuff;
#define TOKENS " \t\n" // used during input file parsing;
const double PI = 3.141592658979323846;

////////////////////////////////////////////////////////////////////////////////
//
// generic functions used for dirty work;
//
////////////////////////////////////////////////////////////////////////////////

FILE *my_fopen( const char *fname, char *mode, int prompt);
char *my_fgets( int &line_number, bool &continue_on_next_line,
                char *buff, int max_size, FILE *fp);
char *my_strtok ( int &line_number, bool &continue_on_next_line,
                  char *ptr, char *tokens, FILE *fp, char *buff);
int my_atoi (const char *str);
int my_get_int( int &line_number, bool &continue_on_next_line,
                char *ptr, FILE *fp, char *buff, int min_val);
double my_get_double( int &line_number, bool &continue_on_next_line,
                      char *ptr, FILE *fp, char *buff,
                      double low_lim, double upp_lim);
int my_read_int_option(int argc, char *argv[], int i);

void *my_malloc(size_t size);
void *my_realloc(void *ptr, size_t size);

////////////////////////////////////////////////////////////////////////////////
//
// RANDOM_NUMBER_GENERATOR
//
////////////////////////////////////////////////////////////////////////////////

class RANDOM_NUMBER_GENERATOR
{
    private:
        long _seed;

    private:
        double sflat01();
        double gauss01();
    public:
        RANDOM_NUMBER_GENERATOR(long seed);
        RANDOM_NUMBER_GENERATOR();
        ~RANDOM_NUMBER_GENERATOR() {}

        double flat_d( double low, double high);
        long flat_l( long low, long high);
        unsigned long flat_ul( unsigned long low, unsigned long high);
        unsigned long long flat_ull( unsigned long long low, unsigned long long high);
        double gauss_mean_d( double mean, double variance);
        long gauss_mean_l( long mean, double variance);
        unsigned long gauss_mean_ul( unsigned long mean, double variance);
        unsigned long long gauss_mean_ull( unsigned long long mean, double variance);
        void set_seed( long seed);
};

#endif
