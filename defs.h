/******************************************************************************\
*								 Definitions							 *
\******************************************************************************/
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <list> 
#define CHAR_LEN 1000
#define EPS1 1.0e-12

using namespace std; 

// Data structures
typedef int allele; 										// type allele
typedef struct {
			allele *chromosome;								
			double fitness;									
} individual;												// data structure individual
typedef struct {			
			individual *ind;
			double sum_fitness;
			double mean_fitness;
			double max_fitness;
			int best_individual;		
} population;												// data structure population


// Global variables and parameters
extern int chrom_size;										// size of the chromosome (dimension of the problem) 
extern int gen;												// generation
extern int local_search;									// with (1) or without (0) local search 		
extern int n_instances;										// number of instances (for NK landscape problem)
extern int n_runs_per_instance;								// number of runs for each instance (for NK landscape problem)
extern int popsize;											// size of the population 
extern int save_datagen_flag;								// flag for saving data for generation in the first run
extern int tau_reset;										// number of generations for reseting part of the population and applying local search
extern int tournament_size;									// size of the pool for tournament selection 
extern int type_crossover;									// crossover type: 1-2PointX; 2-UX; 3-PX, 4-simple ePX, 5-ePX, 6-Multiple 2PointX, 7-Multiple UX								
extern long int max_gen;									// maximum number of generations (when stop criterion is time, it is used for saving some statistics data)
extern long int max_rec_comp, tot_rec_comp;					// Statistics: maximum number and total number of recombination components (for PX and rPX)
extern long int sum_sucRate, sum_impRate, sum_cross;		// Statistics: sum of sucessful recombinations, improvements, crossovers
extern double p_cross;										// crossover rate	
extern double p_mut;										// mutation rate
extern double resetpop_rate;								// controls the percentage of new individuals reset
extern population popold , popnew;							// population (new and old)
// Vectors
extern int *vsort_aux;										// auxiliar sorted vector of integers (used in different methods)
extern int *file_gen;										// data to be stored: number of generations								
extern double *file_best_fitness, *time_run;				// data to be stored: best fitness, runtime
extern double *file_best_fitness_gen;						// data to be stored: best fitness over the generations for run 0
extern double *file_impRate, *file_sucRate;					// data to be stored: improvement and successful crossover rates 
// Matrices
extern int **File_best_ind;									// data to be stored: best individual
extern double **File_recComp;								// data to be stored: recombination partition information
	
// Function declaration
//  aux_functions.cpp
void desaloc_matrixd(double **Matrix , int lines);
void desaloc_matrixi(int **Matrix , int lines);
void rand_perm_size(int *inp, int *out, int size_inp, int size_out);
int isVectorEqual(int *v1, int *v2, int size_v);
int random_int(int L_range, int H_range);
int *aloc_vectori(int lines);
int **aloc_matrixi(int lines , int collums);
double random_dou(void);
double *aloc_vectord(int lines);
double **aloc_matrixd(int lines , int collums);
individual *aloc_vectorind(int lines);
// file_man.cpp
void file_output(int N, int K, int prob_type, int total_runs);
void read_problem(int N, int K, int instance, list<int> *PHI, double **Fl);
// statistics.cpp
void statistics(population *pop, int n_run);
// selection.cpp
int selection(population *pop);
// transformation.cpp
void mutation (allele *parent, allele *offspring);
void Point2X(allele *parent1, allele *parent2, allele *offspring1);
void UX(allele *parent1, allele *parent2, allele *offspring1  );

