/*******************************************************************************\
*  	Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 		  	   	*
*																				*
* 	Genetic Algorithm with epsilon-Partition Crossover (ePX)			*
* 	Test problem: NK Landscapes						 							*
*						 														*
* 	Copyright (C) 2021  Renato Tinos <rtinos@ffclrp.usp.br>						*
* 						 														*
* 	Reference: Tinos, R.; Whitley, D.; Chicano, F. & Ochoa, G. (2021), 		 	*                    
* 				"Partition Crossover for Continuous Optimization: ePX",		 	*
*				Proc. of GECCO'2021.											*
*																				*
* 	ga_epx_nkland is free software: you can redistribute it and/or modify it 	*
* 		under the terms of the GNU General Public License as published by the	*
* 		Free Software Foundation, either version 3 of the License, or			*
* 		(at your option) any later version.						 				*
* 						 														*
* 	ga_epx_nkland is distributed in the hope that it will be useful, but		*
* 		WITHOUT ANY WARRANTY; without even the implied warranty of				*
* 		MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.					*
* 		See the GNU General Public License for more details.					*
* 																				*
* 	You should have received a copy of the GNU General Public License along		*
* 		with this program.  If not, see <http://www.gnu.org/licenses/>.			*		
\*******************************************************************************/

#include <time.h>
#include "defs.h"
#include "Mk.h"


/******************************************************************************\
*				  	Print population					 			 			*
\******************************************************************************/
void print_data(population *pop, int n_run){

	cout <<"Generation:"<< gen << ", run: "<<n_run<<endl;
	cout <<"Best individual:"<< pop->best_individual << endl;
	cout <<"Fitness of the best individual:"<< pop->max_fitness << endl;
	cout <<"Mean fitness: "<< pop->mean_fitness << endl;	
	/*for (int i=0;i<popsize ;i++) {	
		cout <<"("<< pop->ind[i].fitness<<") " ;
		for (int gene=0;gene<chrom_size ;gene++) 
			cout << pop->ind[i].chromosome[gene]<<" ";
		cout << endl;
	}*/
	
}


/******************************************************************************\
*								Fitness (real) Computation    				   *
\******************************************************************************/
double compFitness( Mk *Mk_instance, allele *ind ){
		double Fitness;

		Fitness = Mk_instance->compFitness(ind);		// fitness function for Mk landscapes

		return Fitness;
}


/***********************************************************************************************************\
*								 Crossover	 													     		*
* Crossover type: 1-2PointX; 2-UX; 3-PX, 4-simple ePX, 5-ePX, 6-Multiple 2PointX, 7-Multiple UX			    *
\***********************************************************************************************************/
double crossover(Mk *Mk_instance, allele *parent1, allele *parent2, allele *offspring)
{
	int *offspring_aux, n_rep=10; // number of repetitions for types 6 and 7
	double fit_offspring, fit_offspring_aux, alpha=0.95;

	if (type_crossover==1){
		// 2-point Crossover
		Point2X(parent1,parent2,offspring);
		fit_offspring=compFitness(Mk_instance,offspring);	
	}
	else if (type_crossover==2){
		// Uniform Crossover
		UX(parent1,parent2,offspring);
		fit_offspring=compFitness(Mk_instance,offspring);	
	}	
	else if (type_crossover==3){
		// Partition Crossover
		fit_offspring=Mk_instance->PX(parent1,parent2,offspring);	
	}
	else if (type_crossover==4){
		// simple epsilon-Partition Crossover (i.e., without removing common vertices)
		fit_offspring=Mk_instance->sePX(parent1,parent2,offspring,alpha);	
	}
	else if (type_crossover==5){
		// epsilon-Partition Crossover
		fit_offspring=Mk_instance->ePX(parent1,parent2,offspring,alpha);
	}	
	else if (type_crossover==6){
		// Multiple 2-point Crossover
		offspring_aux=aloc_vectori(chrom_size);
		for (int i=0;i<n_rep;i++){
			if (i==0){
				Point2X(parent1,parent2,offspring);
				fit_offspring=compFitness(Mk_instance,offspring);
			}
			else{
				Point2X(parent1,parent2,offspring_aux);
				fit_offspring_aux=compFitness(Mk_instance,offspring_aux);
				if (i==0 || (fit_offspring_aux>fit_offspring)){
					fit_offspring=fit_offspring_aux;
					for(int gene=0;gene<chrom_size;gene++)
						offspring[gene]=offspring_aux[gene];
						
				}
			}
		}
		delete [] offspring_aux;	
	}
	else if (type_crossover==7){
		// Multiple Uniform Crossover
		offspring_aux=aloc_vectori(chrom_size);
		for (int i=0;i<n_rep;i++){
			if (i==0){
				UX(parent1,parent2,offspring);
				fit_offspring=compFitness(Mk_instance,offspring);
			}
			else{
				UX(parent1,parent2,offspring_aux);
				fit_offspring_aux=compFitness(Mk_instance,offspring_aux);
				if (i==0 || (fit_offspring_aux>fit_offspring) ){
					fit_offspring=fit_offspring_aux;
					for(int gene=0;gene<chrom_size;gene++)
						offspring[gene]=offspring_aux[gene];
						
				}
			}
		}
		delete [] offspring_aux;	
	}	
	
	if ( type_crossover==3 || type_crossover==4 || type_crossover==5 ){
		if (Mk_instance->n_rec_comp > max_rec_comp)
			max_rec_comp=Mk_instance->n_rec_comp;
		tot_rec_comp+=Mk_instance->n_rec_comp;
	}
	
	// Test	
	/*double aux=compFitness(Mk_instance,offspring);
		if ( fabs(aux -fit_offspring)>EPS1){
			cout<<"PROBLEM: FITNESS INCORRECT!"<<endl;	
			cout<<"Fitness: "<<fit_offspring<<", "<<aux<<endl;	
			exit(1);
		}	*/

	return fit_offspring;
	
}


/******************************************************************************\
*								 Generation of the GA						   *
\******************************************************************************/
void generation(Mk *Mk_instance, int n_run){
	int j=0 , parent1, parent2;
	
	do {			
		parent1=selection( &popold );			// Selection of the first parent 
		if ( random_dou () < p_cross ){			
			// Crossover
			parent2=selection( &popold );		// Selection of the second parent 
			if (isVectorEqual(popold.ind[parent1].chromosome,popold.ind[parent2].chromosome,chrom_size)==1){
				// parents are equal
				for (int gene=0;gene<chrom_size;gene++)
					popnew.ind[j].chromosome[gene]=popold.ind[parent1].chromosome[gene];
				popnew.ind[j].fitness=popold.ind[parent1].fitness;
			}
			else {
				// parents are different
				popnew.ind[j].fitness=crossover(Mk_instance, popold.ind[parent1].chromosome , popold.ind[parent2].chromosome,  popnew.ind[j].chromosome);										
				//cout<<"Fitness: i) parents ("<< popold.ind[parent1].fitness<<", "<< popold.ind[parent2].fitness<<"), ii) offspring ("<<popnew.ind[j].fitness<<")"<<endl;				
				// Statistics: crossover rates
				sum_cross++;	// total number of recombinations
				if ( (popnew.ind[j].fitness-popold.ind[parent1].fitness)>EPS1 && (popnew.ind[j].fitness-popold.ind[parent2].fitness)>EPS1 )
					sum_sucRate++;	// Record successful recombination rate
				if ( (popnew.ind[j].fitness-file_best_fitness[n_run])>EPS1 )	
					sum_impRate++; // Record improvement rate
			}				
		}
		else {			
			// Mutation				
			mutation(popold.ind[j].chromosome, popnew.ind[j].chromosome);	
			popnew.ind[j].fitness=compFitness(Mk_instance, popnew.ind[j].chromosome);			
		}
		
		// check if improved best individual
		if (popnew.ind[j].fitness > file_best_fitness[n_run]){		
			file_best_fitness[n_run] = popnew.ind[j].fitness;
			for(int i=0;i<chrom_size;i++)
				File_best_ind[n_run][i]=popnew.ind[j].chromosome[i];
		}
				
		j = j + 1;	
	} while ( j < popsize-1);

	// Elitism (for individual with index j=popsize-1)
	for (int gene=0;gene<chrom_size;gene++)
		popnew.ind[j].chromosome[gene]=popold.ind[popold.best_individual].chromosome[gene];
	popnew.ind[j].fitness=popold.ind[popold.best_individual].fitness;	
	
}	


/******************************************************************************\
*				  	Random individual					 				 	   *
\******************************************************************************/
void randomInd(Mk *Mk_instance, int num_ind){
	
	for (int gene=0;gene<chrom_size;gene++) 
     	popold.ind[num_ind].chromosome[gene] = random_int(0,1);

    popold.ind[num_ind].fitness = compFitness(Mk_instance, popold.ind[num_ind].chromosome);											
	
}


/******************************************************************************\
*				  	Initiate Population 					 				 *
\******************************************************************************/
void initiatePop(Mk *Mk_instance, int n_run){
				
	// Dynamic allocation: populations
	popold.ind = aloc_vectorind(popsize);
	popnew.ind = aloc_vectorind(popsize);

	for (int num_ind=0;num_ind<popsize;num_ind++){
		// Dynamic allocation: chromosomes	
		popold.ind[num_ind].chromosome = aloc_vectori(chrom_size);
		popnew.ind[num_ind].chromosome = aloc_vectori(chrom_size);

		// Random Initialization
		randomInd(Mk_instance,num_ind);	 	
      	
      	if (local_search==1){
      		// Applying Local Search 
			//cout<<"Initialization: Individual- "<<num_ind<<", fitness before local search:"<<popold.ind[num_ind].fitness<<endl;
			popold.ind[num_ind].fitness = Mk_instance->LsFi(popold.ind[num_ind].chromosome,popold.ind[num_ind].fitness); 		// local search: first improvement (for Mk landscapes)				
			//cout<<"...after local search: "<<popold.ind[num_ind].fitness<<endl; 
		}
	}
	file_best_fitness[n_run]=popold.ind[0].fitness;
	for(int i=0;i<chrom_size;i++)
     	File_best_ind[n_run][i]=popold.ind[0].chromosome[i];
	statistics(&popold, n_run);
	//print_data(&popold, n_run);
	
}


/******************************************************************************\
*				  	Copy Population							 			 	   *
\******************************************************************************/
void copy_pop( void ){
		
	for (int ind=0;ind<popsize;ind++) {	
		popold.ind[ind].fitness=popnew.ind[ind].fitness;
		for (int gene=0;gene<chrom_size;gene++) 
			popold.ind[ind].chromosome[gene]=popnew.ind[ind].chromosome[gene];
	}
	
}


/******************************************************************************\
*				  	Run of the GA 			 								   *
\******************************************************************************/
void ga(Mk *Mk_instance, int n_run ){
	
	int j, gen_init=0;
	double time_aux; 
	clock_t time_start;
	
	// Statistics of the number of partitions
	tot_rec_comp=0;
	max_rec_comp=0;
	// Statistics: crossover
	sum_sucRate=0;
	sum_impRate=0;
	sum_cross=0;
		
	// Initialization
	time_start=clock();	
	gen=0;
	initiatePop(Mk_instance, n_run);				// initiating population
	
	do {
		gen++; 										// generation index
		if ((gen-gen_init)>tau_reset){			
			//cout<<"gen: "<<gen<<endl;	
			//cout<<"Reseting part of the population and applying local search...."<<endl;
			gen_init=gen;
			// Elitism
			for (int gene=0;gene<chrom_size;gene++)
				popold.ind[0].chromosome[gene]=popold.ind[popold.best_individual].chromosome[gene];
			popold.ind[0].fitness=popold.ind[popold.best_individual].fitness;
			popold.best_individual=0; 
			// Reseting part of the population
			j=(int) (resetpop_rate*popsize) + 1;
			if (j>popsize)
				j=popsize;
			for (int num_ind=1;num_ind<j;num_ind++)
				randomInd(Mk_instance, num_ind);			
			if (local_search==1){	
				// Applying local search
				for (int num_ind=0;num_ind<popsize;num_ind++){
				    //cout<<"Initialization: Individual- "<<num_ind<<", fitness before local search:"<<popold.ind[num_ind].fitness<<endl;
					popold.ind[num_ind].fitness = Mk_instance->LsFi(popold.ind[num_ind].chromosome,popold.ind[num_ind].fitness); 		// local search: first improvement (for Mk landscapes)				
					//cout<<"...after local search: "<<popold.ind[num_ind].fitness<<endl; 
				}
			}
			statistics(&popold,n_run);
		}
		
		generation(Mk_instance,n_run);
			
		copy_pop();						// popold=popnew
		
		statistics(&popold,n_run);	
		//print_data(&popold,n_run);
					
		time_aux = ( (double) ( clock() - time_start ) ) / ( (double) CLOCKS_PER_SEC);

	} while ( time_aux < ((double) (Mk_instance->M*Mk_instance->k)/10.0) );
	//} while ( gen < max_gen );
	
	//time_aux = ( (double) ( clock() - time_start ) ) / ( (double) CLOCKS_PER_SEC);
	// Data to be saved
	time_run[n_run]=time_aux;
	file_gen[n_run]=gen;
	if (gen<max_gen-1){
		// Save data for crossover statistics (if not saved function statistics() yet )
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
	
	// Deleting population
	for (int num_ind=0;num_ind<popsize;num_ind++){
		delete [] popold.ind[num_ind].chromosome;
		delete [] popnew.ind[num_ind].chromosome;
	}
	delete [] popold.ind;
	delete [] popnew.ind;

}


/******************************************************************************\
*				  	NK model  												   *
\******************************************************************************/
void NK(int N, int K, int model_neig, int instance, list<int> *PHI, double **Fl) 
{
	int i, n_comb_fl, *v_aux, q_NK_1; 
	
	n_comb_fl=(int) pow(2.0,K+1);   					// Fl is a matrix N x n_comb_fl

	if (model_neig==1){
		// adjacent NK model
		read_problem(N, K, instance, PHI, Fl);			// read Fl and VIG for adjacent model file
	}
	else{
		if (model_neig==2 || model_neig==3){
			// random NKq model
			if (model_neig==2)
				q_NK_1=n_comb_fl/2;
			else
				q_NK_1=n_comb_fl;
			q_NK_1=q_NK_1-1;			
			// matrix Fl
			for (int l=0;l<N;l++)
				for (int j=0;j<n_comb_fl;j++)
					Fl[l][j]=(double) random_int(0, q_NK_1)/q_NK_1;	// random fl 
		}
		else{
			// random NK model
			// matrix Fl
			for (int l=0;l<N;l++)
				for (int j=0;j<n_comb_fl;j++)
					Fl[l][j]=random_dou();					// random fl 
		}
		
	    // vector of lists PHI (PHI[l] indicates the variables that influences subfunction fl)
		v_aux=aloc_vectori(K+1);
		for (int l=0;l<N;l++){
			PHI[l].push_back(l); 						// index l is always in PHI[l]
			rand_perm_size(vsort_aux, v_aux, N, K+1);	
			i=0;	
			for (int j=0;j<K;j++){
				if (v_aux[i]==l)
					i++;				
				PHI[l].push_back(v_aux[i]); 			// random index in PHI[l]
				i++;
			}
		}
		delete [] v_aux;
	}

}


/******************************************************************************\
*				  	Main													   *
\******************************************************************************/
int main(int argc , char *argv[])
{
	int total_runs, n_run=0;
	int N_NK, K_NK, model_NK;					// variables for the NK model
	double **Fl_NK;								// 	matrix with the values for the subfunctions fl in the NK model
	list<int> *PHI_NK;							// 	vector of lists PHI (PHI[l] indicates the variables that influences subfunction fl)
	
	// Arguments
	if( argc < 5) {
		cout<<"Insufficient number of arguments!"<<endl;
		cout<<"Call: ga_epx_nkland <N> <K> <model_NK> <crossover_type>"<<endl;
		exit(1);
	}
	else{
		N_NK=atoi(argv[1]);
		K_NK=atoi(argv[2]);	
		model_NK=atoi(argv[3]);
		type_crossover=atoi(argv[4]);
		if ( N_NK<1 || K_NK<0 || K_NK>N_NK-1 || model_NK<0 || model_NK>3 || type_crossover<1 || type_crossover>7 ) {
			cout<<"Incorrect arguments!"<<endl;
			cout<<"Call: ga_epx_nkland  <N> (N>0) <K> (0<=K<=N) <model_NK> (0-random NK; 1-regular NK; 2,3-random NKq) <type_crossover> (1<=type_crossover<=7)"<<endl;
			exit(1);
		}
	}	
	
	// Parameters
	chrom_size=N_NK;													// size of the chromosome
	p_mut=3.0/chrom_size;												// mutation rate
	total_runs=n_instances*n_runs_per_instance;							// number of runs
	max_gen=chrom_size*5;												// maximum number of generations (when stop criterion is time, it is used for saving some statistics data)	
			
	// Allocation of vectors and matrices
	Fl_NK = aloc_matrixd ( N_NK , (int) pow(2.0,K_NK+1) );
	file_best_fitness_gen=aloc_vectord(max_gen);
	file_best_fitness=aloc_vectord(total_runs);
	file_gen=aloc_vectori(total_runs);
	time_run=aloc_vectord(total_runs);
	file_sucRate=aloc_vectord(total_runs);
	file_impRate=aloc_vectord(total_runs);
	File_best_ind=aloc_matrixi(total_runs,chrom_size);
	File_recComp=aloc_matrixd(total_runs,2);
	vsort_aux=aloc_vectori(chrom_size);									// Auxiliar sorted vector of integers (used in different methods)
	for (int i=0;i<chrom_size;i++)
		vsort_aux[i]=i;
	
	cout << "\n ***** Hybrid Genetic Algorithm ****" << endl;
	cout << "N="<<N_NK <<", K="<<K_NK<< endl;
	cout << "type_crossover="<<type_crossover<< endl;

	for (int instance=0;instance<n_instances;instance++) {	
		cout <<"Instance: "<< instance << endl;
		srand(instance);												// random seed  (for instance)		
		// Creating the NK Model
		PHI_NK = new list<int>[N_NK]; 
		NK(N_NK, K_NK, model_NK, instance, PHI_NK, Fl_NK);				// create matrices VIG_NK and Fl_NK

		// Creating the Mk Model	
		Mk *Mk_instance = new Mk(N_NK, N_NK, K_NK+1, PHI_NK, Fl_NK);	// from class Mk (Mk.h)
		//Mk_instance->print();											// print graph and costs of the elements for the NK Landscape
		for (int i_run=0;i_run<n_runs_per_instance;i_run++){
			cout <<"Run:"<< n_run <<", random seed: " << n_instances+n_run << endl;
			srand( n_instances+n_run );									// random seed  (for run)
			ga(Mk_instance, n_run);					    				// run GA
			n_run++;
		}

		delete [] PHI_NK;
		delete Mk_instance;
	}
	
	file_output(N_NK,K_NK,model_NK,total_runs);					// save data

	// Desallocation of vectors and matrices	
	desaloc_matrixd (Fl_NK,N_NK);
	desaloc_matrixd (File_recComp,total_runs);
	desaloc_matrixi (File_best_ind,total_runs);
	delete [] time_run;
	delete [] file_gen;
	delete [] file_best_fitness;
	delete [] file_best_fitness_gen;
	delete [] file_sucRate;
	delete [] file_impRate;	
	delete [] vsort_aux;
		
	return 0;
}
