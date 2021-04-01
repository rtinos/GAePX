#include "defs.h"

// Global variables and parameters
int chrom_size;										// size of the chromosome (dimension of the problem) 
int gen;											// generation
int local_search=1;									// with (1) or without (0) local search 
int n_instances = 10;								// number of instances (for NK landscape problem)
int n_runs_per_instance = 10;						// number of runs for each instance (for NK landscape problem)
int popsize=100;									// size of the population 
int save_datagen_flag=0;							// flag for saving data for generation in the first run
int tau_reset=50;									// number of generations for reseting part of the population and applying local search
int tournament_size=3;								// size of the pool for tournament selection 
int type_crossover;									// crossover type: 1-2PointX; 2-UX; 3-PX, 4-simple ePX, 5-ePX, 6-Multiple 2PointX, 7-Multiple UX		
long int max_gen;									// maximum number of generations (when stop criterion is time, it is used for saving some statistics data)
long int max_rec_comp, tot_rec_comp;				// Statistics: maximum number and total number of recombination components (for PX and rPX)
long int sum_sucRate, sum_impRate, sum_cross;		// Statistics: sum of sucessful recombinations, improvements, crossovers
double p_cross=0.6;									// crossover rate											
double p_mut;										// mutation rate
double resetpop_rate=0.99;							// controls the percentage of new individuals reset
population popold , popnew;							// population (new and old)
// Vectors
int *vsort_aux;										// auxiliar sorted vector of integers (used in different methods)
int *file_gen;										// data to be stored: number of generations								
double *file_best_fitness, *time_run;				// data to be stored: best fitness, runtime
double *file_best_fitness_gen;						// data to be stored: best fitness over the generations for run 0
double *file_impRate, *file_sucRate;				// data to be stored: improvement and successful crossover rates 
// Matrices
int **File_best_ind;								// data to be stored: best individual
double **File_recComp;								// data to be stored: recombination partition information
