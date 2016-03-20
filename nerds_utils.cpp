#include "nerds_utils.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

////////////////////////////////////////////////////////////////////////////////
//
// Utility functions used during input files processing;
//
////////////////////////////////////////////////////////////////////////////////

FILE *my_fopen( const char *fname, char *mode, int prompt) 
{
    FILE *fp;
    while (1) {
        if ( prompt) 
            scanf( "%s", fname);
        if ( (fp = fopen(fname, mode)) != NULL)
            break;
        printf("Error opening file %s for %s access.\n", fname, mode);
        if ( !prompt) exit(1);
        printf("Please enter another filename.\n");
    }
    return fp;
}

char *my_fgets( int &line_number, bool &continue_on_next_line,
                char *buff, int max_size, FILE *fp) 
{
    // get an entire line from the file, update the line number and prune
    // any comment; a \ at the end of a line with no comment part (#) means
    // to continue;

    char *line_val;
    int i = 0;
 
    continue_on_next_line = false;
    line_number ++;
    
    line_val = fgets(buff, max_size, fp);

    if (line_val == NULL) return (line_val);
    
    // check the length of this line;
    for ( i = 0; i < max_size; i++) {
        if ( buff[i] == '\n') 
            break;
        if ( buff[i] == '\0') {
            printf("Error. Line number %d is too long while reading.\n", line_number);
            printf("All lines must be at most %d characters long.\n", BUFFER_SIZE - 2);
            exit (1);
        }
    }
    
    
    for ( i = 0; i < max_size && buff[i] != '\0'; i++) {
        if ( buff[i] == '#') {
             buff[i] = '\0';
             break;
        }
    }
    
    if ( i < 2 ) return (line_val);
    if ( buff[i-1] == '\n' && buff[i-2] == '\\') { 
        continue_on_next_line = true;
        buff[i-2] = '\n';
        buff[i-1] = '\0';
    }
    return (line_val);
}

char *my_strtok( int &line_number, bool &continue_on_next_line,
                  char *ptr, char *tokens, FILE *fp, char *buff) 
{
    char *line_val;
    line_val = strtok( ptr, tokens);
    
    while (1) {
        if ( line_val != NULL || continue_on_next_line == false)
            return (line_val);
        if ( my_fgets( line_number, continue_on_next_line,
                       buff, BUFFER_SIZE, fp) == NULL) 
            return (NULL);
        line_val = strtok( buff, tokens);
    }
}

int my_atoi( const char *str) 
{
    if ( str[0] < '0' || str[0] > '9') 
        return (-1);
    return atoi(str);
}

int my_get_int( int &line_number, bool &continue_on_next_line,
                char *ptr, FILE *fp, char *buff, int min_val) 
{
    int val = 0;

    ptr = my_strtok( line_number, continue_on_next_line, NULL, TOKENS, fp, buff);
    if ( ptr == NULL) {
        printf("Error: missing value on line %d.\n", line_number);
        exit(1);
    }
    val = my_atoi(ptr);
    if ( val < min_val) {
        printf("Error: Incorrect value on line %d. Should be at least %d\n",
               line_number, min_val);
        exit(1);
    }   
    ptr = my_strtok( line_number, continue_on_next_line, NULL, TOKENS, fp, buff);
    if ( ptr != NULL) {
        printf("Error: Extra characters at end of line %d.\n", line_number);
        exit (1);
    }
    
    return val;
}

double my_get_double( int &line_number, bool &continue_on_next_line,
                      char *ptr, FILE *fp, char *buff,
                      double low_lim, double upp_lim) 
{
    double val = 0.0;

    ptr = my_strtok( line_number, continue_on_next_line, NULL, TOKENS, fp, buff);
    if ( ptr == NULL) {
        printf("Error: missing value on line %d.\n", line_number);
        exit(1);
    }
    val = atof(ptr);
    if ( val <= low_lim || val > upp_lim) {
        printf("Error: Incorrect value on line %d. Should be within (%f, %f]\n",
               line_number, low_lim, upp_lim);
        exit(1);
    }
    
    ptr = my_strtok( line_number, continue_on_next_line, NULL, TOKENS, fp, buff);
    if ( ptr != NULL) {
        printf("Error: Extra characters at end of line %d.\n", line_number);
        exit (1);
    }
    
    return val;
}

int my_read_int_option( int argc, char *argv[], int i) 
{
    int value = 0, num_read = 0;
    if (argc > i + 1) {
        num_read = sscanf(argv[i + 1], "%d", &value);
    }
    if (num_read != 1) {       
        printf("Error:  %s option requires an integer value.\n\n", argv[i]);
        exit(1);
    }
    return (value);
}

void *my_malloc(size_t size) 
{
    void *ret;
    if ((ret = malloc (size)) == NULL) {
        fprintf(stderr,"Error:  Unable to malloc memory.  Aborting.\n");
        abort ();
        exit (1);
    }
    return (ret);
}

void *my_realloc(void *ptr, size_t size) 
{
    void *ret;
    if ((ret = realloc (ptr,size)) == NULL) {
        fprintf(stderr,"Error:  Unable to realloc memory.  Aborting.\n");
        exit (1);
    }
    return (ret);
}

////////////////////////////////////////////////////////////////////////////////
//
// RANDOM
//
////////////////////////////////////////////////////////////////////////////////

RANDOM_NUMBER_GENERATOR::RANDOM_NUMBER_GENERATOR(long seed) :_seed(seed)
{
    srandom( _seed);
}

RANDOM_NUMBER_GENERATOR::RANDOM_NUMBER_GENERATOR() : _seed(1)
{
    srandom( _seed);
}

double RANDOM_NUMBER_GENERATOR::sflat01()
{
    double val = random() * 1.0 / RAND_MAX;
    return val;
}

void RANDOM_NUMBER_GENERATOR::set_seed(long seed) 
{
    _seed = seed;
    srandom( seed);
}

double RANDOM_NUMBER_GENERATOR::flat_d(double low, double high) 
{
    assert( low < high);
    return (high - low) * sflat01() + low;
}

long RANDOM_NUMBER_GENERATOR::flat_l(long low, long high)
{
    assert( low < high);
    long val = long((high - low) * sflat01() + low);
    assert( val >= low && val < high);
    return val;
}

unsigned long RANDOM_NUMBER_GENERATOR::flat_ul(unsigned long low, 
    unsigned long high) 
{
    assert( low < high);
    unsigned long val = (unsigned long)((high-low) * sflat01() + low);
    assert( val >= low && val < high);
    return val;
}

unsigned long long RANDOM_NUMBER_GENERATOR::flat_ull(unsigned long long low,
    unsigned long long high) 
{
    assert( low < high);
    unsigned long long val = (unsigned long long)((high-low) * sflat01() + low);
    assert( val >= low && val < high);
    return val;
}

double RANDOM_NUMBER_GENERATOR::gauss01()
{
    // mean = 0, variance = 1;
    double modifier;
    double compile_b;
    double in_a, in_b;
    double out_a;
    static double out_b;
    static int recalc = 1;

    if (recalc) {
        // Range from (0:1], not [0:1). Had to change this to prevent log(0).
        in_a = 1.0 - sflat01();
        in_b = sflat01();

        modifier = sqrt( -2.0 * log(in_a));
        compile_b = 2.0 * PI * in_b;

        out_a = modifier * cos(compile_b);
        out_b = modifier * sin(compile_b);

        recalc = 0;
        return out_b;
    }

    recalc = 1;
    return out_b;
}

double RANDOM_NUMBER_GENERATOR::gauss_mean_d(double mean, double variance) 
{
    double ret = gauss01() * variance + mean;
    return ret;
}

long RANDOM_NUMBER_GENERATOR::gauss_mean_l(long mean, double variance) 
{
    return ((long)(gauss01() * variance + mean));
}

unsigned long RANDOM_NUMBER_GENERATOR::gauss_mean_ul(unsigned long mean, 
    double variance) 
{
    return((unsigned long)(gauss01() * variance + mean));
}

unsigned long long RANDOM_NUMBER_GENERATOR::gauss_mean_ull(unsigned long long mean, 
    double variance) 
{
    return ((unsigned long long)(gauss01() * variance + mean));
}
