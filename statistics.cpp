/******************************************************************************\
*								 Statistics						 *
\******************************************************************************/
#include "defs.h"


/******************************************************************************\
*				  	Population Statistics 									   *
\******************************************************************************/
void statistics(population *pop, int n_run)
{	

	pop->sum_fitness = pop->ind[0].fitness; 			// sum of the fitness in the population
	pop->max_fitness = pop->ind[0].fitness;   			// maximum fitness in the population
	pop->best_individual = 0;							// best individual in the population
	if (pop->ind[0].fitness > file_best_fitness[n_run]){
			file_best_fitness[n_run] = pop->ind[0].fitness;
			for(int i=0;i<chrom_size;i++) 
				File_best_ind[n_run][i]=pop->ind[0].chromosome[i];
	}		
	for(int j=1;j<popsize;j++) {
		pop->sum_fitness = pop->sum_fitness + pop->ind[j].fitness;
		if (pop->ind[j].fitness > pop->max_fitness )	{	
			pop->max_fitness = pop->ind[j].fitness; 
			pop->best_individual = j;			
		}
		if (pop->ind[j].fitness > file_best_fitness[n_run]){		
			file_best_fitness[n_run] = pop->ind[j].fitness;
			for(int i=0;i<chrom_size;i++)
				File_best_ind[n_run][i]=pop->ind[j].chromosome[i];
		}
	}

	pop->mean_fitness = pop->sum_fitness / popsize; 	// mean fitness in the population
	
	// Save data for fitness along the generations: only for the first run
	if (save_datagen_flag==1 && n_run==0){			
		if (gen<max_gen)			
			file_best_fitness_gen[gen]=pop->max_fitness;				
	}
	
	// Save data for crossover statistics (only for first max_gen generations)
	if (gen==max_gen-1){
		if (sum_cross>0){
			file_sucRate[n_run]=((double) sum_sucRate)/sum_cross;
			file_impRate[n_run]=((double) sum_impRate)/sum_cross;
		}
		else{
			file_sucRate[n_run]=0.0;
			file_impRate[n_run]=0.0;		
		}
		if (type_crossover==3 || type_crossover==4 || type_crossover==5){
			if (sum_cross>0)
				File_recComp[n_run][0]=((double) tot_rec_comp)/sum_cross;
			else
				File_recComp[n_run][0]=0.0;
			File_recComp[n_run][1]=(double) max_rec_comp;
		}
	}	
	
}

