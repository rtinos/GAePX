#include "defs.h"


/******************************************************************************\
*								 2-point Crossover	 						 *
\******************************************************************************/
void Point2X(allele *parent1, allele *parent2, allele *offspring){
		int p1, p2, aux;
		double aux_d;

		// defining the crossover points
		p1 =random_int(0,chrom_size-1);		// point 1
		p2 =random_int(0,chrom_size-1);		// point 2
		while (p1 == p2)
			p2 =random_int(0,chrom_size-1);	// point 2
		if (p1>p2) {
			aux=p1;
			p1=p2;
			p2=aux;
		}											 
		// generating the offspring
		aux_d=random_dou();
		if (aux_d<0.5){
			for (int gene=0;gene<chrom_size;gene++) {
				if (gene<p1 || gene>=p2) 
					offspring[gene] = parent1[gene];	
				else
					offspring[gene] = parent2[gene];
			}
		}
		else{
			for (int gene=0;gene<chrom_size;gene++) {
				if (gene<p1 || gene>=p2) 
					offspring[gene] = parent2[gene];	
				else
					offspring[gene] = parent1[gene];
			}
		}
		
}


/******************************************************************************\
*								 Uniform Crossover	 						 *
\******************************************************************************/
void UX(allele *parent1, allele *parent2, allele *offspring1  ){
		int aux;

		for (int gene=0;gene<chrom_size;gene++){
			if (parent1[gene]==parent2[gene])
				aux=0;
			else
				aux=random_int (0,1);				// mask: define if the gene comes from parent 1 (0) or 2 (1)				
			if (aux==0)
				offspring1[gene] = parent1[gene];			
			else
				offspring1[gene] = parent2[gene];
		}	
				
}


/******************************************************************************\
*								 Mutation														   *
\******************************************************************************/
void mutation (allele *parent, allele *offspring){
	
	for (int gene=0;gene<chrom_size;gene++){
		if ( random_dou () < p_mut ){
			if (parent[gene]==0)
				offspring[gene]=1;
			else
				offspring[gene]=0;
		}
		else{
			offspring[gene]=parent[gene];
		}
	}	

}


