/******************************************************************************\
*								 Selection									 *
\******************************************************************************/
#include "defs.h"


/******************************************************************************\
*								Tournament Selection							 *
\*******************************************************************************/
int selection_tournament(population *pop)
{
	int individual_rand, individual_chosen;
	
	individual_chosen=random_int(0, (popsize-1) );
	for (int i=1;i<tournament_size;i++){
		individual_rand=random_int(0, (popsize-1) );
		if ( pop->ind[individual_rand].fitness > pop->ind[individual_chosen].fitness  )
			individual_chosen=individual_rand;
	}
		
	return individual_chosen;
}


/******************************************************************************\
*								Selection										 *
\*******************************************************************************/
int selection(population *pop) 
{
	int individual_chosen;
	
	individual_chosen = selection_tournament( pop );

 	return individual_chosen;
}



