/******************************************************************************\
*				  				 Files Manipulation							 *
\******************************************************************************/

#include "defs.h"
#include<cstring>
#include<fstream>


/*******************************************************************\
*	 Read the problem instance  								    *
* Obs.: the order (big endian) for the variables in the 		    *
*		file is different from the order (little endian) used here  *
\*******************************************************************/
void read_problem(int N, int K, int instance, list<int> *PHI, double **Fl){
	  int k=0, n_sub, *v_aux;	  
	  char line[CHAR_LEN], *keywords, Delimiters[] = " :=\n\t\r\f\v", name[CHAR_LEN], *name_p;
		
	  v_aux=aloc_vectori(K+1);
	  n_sub=(int) pow(2.0,K+1);
	  name_p = name;
	  sprintf(name,"prob/N%d/K%d/solver_fitcont_%d.dat",N,K,instance);
	
	  ifstream fin(name_p);
	  while((fin.getline(line, CHAR_LEN-1))){
			if(!(keywords = strtok(line, Delimiters)))
	  			continue;
			if(!strcmp(keywords, "m")){
				for(int i=0; i<K+1; i++)
					v_aux[i] = atoi(strtok(NULL, Delimiters));
				for(int i=0; i<K+1; i++)
					PHI[k].push_back(v_aux[K-i]);				
	  			for(int i=0; i<n_sub; i++)
					fin>>Fl[k][i];
				k++;
			}
			
	  }
	  fin.close();
	  
	  delete [] v_aux;
}


/******************************************************************************\
* 					Save  : end of the simulation						   *
\******************************************************************************/
void file_output(int N, int K, int model_NK, int total_runs)
{
	char *name_p;
	char name[CHAR_LEN];
	FILE *Bestfit_file, *Bestind_file, *Time_file, *Gen_file, *sucRate_file, *impRate_file;

    name_p = name;
		
  	// Best fitness in each generation for run 0
  	if (save_datagen_flag==1){
  		FILE *Bfg_file;
		sprintf(name,"bfg_N%d_K%d_p%d_c%d.dat",N,K,model_NK,type_crossover);
		if ((Bfg_file = fopen(name_p,"w"))==NULL) {
			puts("The file bfg to be saved cannot be open \n");
			exit(1);
		}
		for (int i=0;i<max_gen;i++) {
			fprintf(Bfg_file,"%.3f ",file_best_fitness_gen[i]);
		}
		fclose(Bfg_file);
	}

    // Best fitness 
	sprintf(name,"bfi_N%d_K%d_p%d_c%d.dat",N,K,model_NK,type_crossover);
	if ((Bestfit_file = fopen(name_p,"w"))==NULL) {
		puts("The file bfi to be saved cannot be open \n");
		exit(1);
	}
	for (int i=0;i<total_runs;i++) {
		fprintf(Bestfit_file,"%.14f ",file_best_fitness[i]);
	}
	fclose(Bestfit_file);
		
	 // Best individuals
	sprintf(name,"bind_N%d_K%d_p%d_c%d.dat",N,K,model_NK,type_crossover);
	if ((Bestind_file = fopen(name_p,"w"))==NULL) {
		puts("The file bind to be saved cannot be open \n");
		exit(1);
	}
	for (int i=0;i<total_runs;i++) {
		for (int gene=0;gene<chrom_size;gene++)
			fprintf(Bestind_file,"%d ",File_best_ind[i][gene]);
		fprintf(Bestind_file,"\n");
	}
	fclose(Bestind_file);
	
  	// Time for each run
	sprintf(name,"time_N%d_K%d_p%d_c%d.dat",N,K,model_NK,type_crossover);
	if ((Time_file = fopen(name_p,"w"))==NULL) {
		puts("The file time to be saved cannot be open \n");
		exit(1);
	}
	for (int i=0;i<total_runs;i++) {
		fprintf(Time_file,"%.2f ",time_run[i]);
	}
	fclose(Time_file);
	
	// Number of generations for each run
	sprintf(name,"gen_N%d_K%d_p%d_c%d.dat",N,K,model_NK,type_crossover);
	if ((Gen_file = fopen(name_p,"w"))==NULL) {
		puts("The file gen to be saved cannot be open \n");
		exit(1);
	}
	for (int i=0;i<total_runs;i++) {
		fprintf(Gen_file,"%d ",file_gen[i]);
	}
	fclose(Gen_file);		
		
	// Crossover Statistics: successful recombination rate 
	sprintf(name,"sucRate_N%d_K%d_p%d_c%d.dat",N,K,model_NK,type_crossover);
	if ((sucRate_file = fopen(name_p,"w"))==NULL) {
		puts("The file sucRate to be saved cannot be open \n");
		exit(1);
	}
	for (int i=0;i<total_runs;i++)		
		fprintf(sucRate_file,"%.8f ",file_sucRate[i]);		
	fclose(sucRate_file);
	
	// Crossover Statistics: improvement rate 
	sprintf(name,"impRate_N%d_K%d_p%d_c%d.dat",N,K,model_NK,type_crossover);
	if ((impRate_file = fopen(name_p,"w"))==NULL) {
		puts("The file impRate to be saved cannot be open \n");
		exit(1);
	}
	for (int i=0;i<total_runs;i++)		
		fprintf(impRate_file,"%.8f ",file_impRate[i]);		
	fclose(impRate_file);
	
	if (type_crossover==3 || type_crossover==4 || type_crossover==5){
			FILE *recComp_file;
			sprintf(name,"recComp_N%d_K%d_p%d_c%d.dat",N,K,model_NK,type_crossover);
			if ((recComp_file = fopen(name_p,"w"))==NULL) {
				puts("The file recComp to be saved cannot be open \n");
				exit(1);
			}
			for (int i=0;i<total_runs;i++)	
				fprintf(recComp_file,"%.3f ",File_recComp[i][0]);
			fprintf(recComp_file,"\n");
			for (int i=0;i<total_runs;i++)	
				fprintf(recComp_file,"%.0f ",File_recComp[i][1]);				
			fclose(recComp_file);		
	}
	
}
