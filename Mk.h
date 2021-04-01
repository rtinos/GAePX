/******************************************************************************\
*								 Class: Mk Landscapes					 	   *
* Problem: NK Landscapes													   *
\******************************************************************************/
#include <cmath>
#include <cstdlib>
#include "graph.h"																// Graph class
#define EPS2 1.0e-12

class Mk {
	private:	
		list<int> *PHI; 														// Vector of lists PHI (PHI[l]: indicates the variables that influences subfunction fl)	
		list<int> *PSI; 														// Vector of lists PSI (PSI[i]: indicates the indices of the subfunction fl influenced by the i-th variable)	
		list<int> *VIG; 														// Vector of lists VIG (VIG[i]: indicates the variables that interact with i-th variable in subfunctions fl)	
		double **Fl;															// Matrix with the contributions for each subfunction fl
		typedef struct{
			int size;	   															// number of vertices in the recombining component
			double cost_parent1;													// cost of the candidate partition (if it is a true recombining partition) - blue solution
			double cost_parent2;													// cost of the candidate partition (if it is a true recombining partition) - red solution
			list<int> var_id; 														// List that indicates variables in each component	
		} recomb_comp;															// recombining components for partition crossover
		void dec2binvec(int number, int *x, int size_x);						// Transform an integer into a binary vector with lenght size_x 	
		void dec2binvec_aux(int number, int *x, int index); 					// Auxiliar: dec2binvec		
		void randPerm(int *v, int size_v);										// Random permutation 	
		int binvec2dec(int *v, int size_v);										// Transform a binary vector into a decimal (big endian notation)
		int randomInt(int L_range, int H_range);			 					// Random integer between L_range and H_range
		double dfElement (int *x, int var_index);								// Fitness difference caused by changing bit
	public:
	    int N;																	// Number of decision variables
	    int M;																	// Number of subfunctions
		int k;																	// Epistasis degree
		int n_rec_comp;															// number of recombining components
		Mk(int N, int M, int k, list<int> *PHI, double **Fl);
	    ~Mk(void);
	    void print(void);														// Print Mk information
	    double compFitness(int *x);												// Compute Fitness
		double LsFi(int *x, double f);											// Local search: first improvement	
		double PX(int *parent1, int *parent2, int *offspring);					// Partition crossover (PX)
		double sePX(int *parent1, int *parent2, int *offspring, double alpha);	// Simple  epsilon - partition crossover (sePX)
		double ePX(int *parent1, int *parent2, int *offspring, double alpha);	// epsilon - partition crossover (ePX)
};


/******************************************************************************\
*								Constructor									   *
\******************************************************************************/
Mk::Mk(int N, int M, int k, list<int> *PHI, double **Fl){
	int i, i_next, aux_flag, var_index, var_index_next, *v_aux;
	list<int>::iterator ii;
	
	// Parameters for Mk landscape
	this->N = N;
	this->k = k;
	this->M = M;
				
	// Vector of lists PHI (PHI[l]: indicates the variables that influences subfunction fl)		
	this->PHI = PHI;

	// Vector of lists PSI (PSI[i]: indicates the index of the subfunction fl influenced by the i-th variable)
	PSI = new list<int>[N]; 
	for (int l=0;l<M;l++){
		for(list<int>::iterator j = PHI[l].begin(); j != PHI[l].end(); j++){
			var_index=*j;
			PSI[var_index].push_back(l);
		}
	}	
		
	// Vector of lists VIG (VIG[i]: indicates the variables that interact with i-th variable in subfunctions fl)
	VIG = new list<int>[N]; 
	v_aux=new int[k];
	for (int l=0;l<M;l++){
		i=0;
		for(list<int>::iterator j = PHI[l].begin(); j != PHI[l].end(); j++){
			v_aux[i]=*j;
			i++;
		}
		for (i=0;i<k;i++){
			var_index=v_aux[i];
			i_next=i+1;
			while (i_next<k){
				var_index_next=v_aux[i_next];
				aux_flag=0;
				if ( !VIG[var_index].empty() && !VIG[var_index_next].empty() ){
					ii=VIG[var_index].begin();
					while (ii != VIG[var_index].end() && *ii != var_index_next)					
						ii++;
					if (ii != VIG[var_index].end())	
						aux_flag=1;
				}
				if (aux_flag==0){
					VIG[var_index].push_back(var_index_next); 
					VIG[var_index_next].push_back(var_index); 	
				}
				i_next++;
			}
		}
	}
	delete [] v_aux;
				
	// Matrix with the contributions for each f_l
	this->Fl = Fl;
		
}


/******************************************************************************\
*								 Destructor									   *
\******************************************************************************/
Mk::~Mk(void){
	
	delete [] VIG; 
	delete [] PSI;  
	          
}  


/******************************************************************************\
* Transform a binary vector into a decimal 									   *
\******************************************************************************/
int Mk::binvec2dec(int *v, int size_v){
	int y=0, base=1;
	
	for (int i=0;i<size_v;i++){	
		y += v[i]*base;
		base=base*2;	
	}
	
	return( y );  
}


/******************************************************************************\
*							Compute Fitness									   *
\******************************************************************************/
double Mk::compFitness (int *x){
	int i, *v_aux;
	double f=0.0;

	v_aux = new int[k];
	for (int l=0;l<M;l++){
		i=0;
		for(list<int>::iterator j = PHI[l].begin(); j != PHI[l].end(); j++){	
			v_aux[i]=x[*j];
			i++;
		}
		f = f + Fl[l][binvec2dec(v_aux,k)];
	}
	f = f/N;

	delete [] v_aux;
	return f;
}


/******************************************************************************\
*			Auxiliar: dec2binvec				   							   *
\******************************************************************************/
void Mk::dec2binvec_aux(int number, int *x, int index)
{
	int remainder;

	if (number <= 1) 	{
		x[index]=number;
		return;
	}
	remainder = number%2;	
	x[index]=remainder;
	dec2binvec_aux(number >> 1, x, index+1);
	
}


/******************************************************************************\
*			  Transform an integer into a binary vector with lenght size_x 	   *
\******************************************************************************/
void Mk::dec2binvec(int number, int *x, int size_x)
{
	
	for (int i=0;i<size_x;i++)
		x[i]=0;

	dec2binvec_aux(number, x, 0);

}


/******************************************************************************\
*						 Random Integer between L_range and H_range			   *
\******************************************************************************/
int Mk::randomInt(int L_range, int H_range)
{
	return(  (int) ( (rand()/(RAND_MAX+1.0))*(H_range-L_range+1)+L_range ) );  	// random integer beteween [L_range and H_range]
}


/******************************************************************************\
* 			Random Permutation 				 								   *
\******************************************************************************/
void Mk::randPerm(int *v, int size_v)
{
	int j, aux, *v_aux;
	
	v_aux=new int[size_v];
	
	for(int i=0;i<size_v;i++)
		v_aux[i]=i;
	
	for(int i=0;i<size_v-1;i++) {
		j= randomInt(i,size_v-1);  		
		aux=v_aux[i];
		v_aux[i]=v_aux[j];
		v_aux[j]=aux;
		v[i]=v_aux[i];
	}
	v[size_v-1]=v_aux[size_v-1];
	
	delete [] v_aux;
}


/******************************************************************************\
*		Fitness difference caused by changing bit 							   *
\******************************************************************************/
double Mk::dfElement (int *x, int var_index){
	int i, l, i_sel, *v_aux;
	double df=0.0;

	if (!PSI[var_index].empty()){
		v_aux = new int[k];
		
		for(list<int>::iterator li = PSI[var_index].begin(); li != PSI[var_index].end(); li++){
			l=*li;
			i=0;
			for(list<int>::iterator j = PHI[l].begin(); j != PHI[l].end(); j++){	
				if (*j==var_index)
					i_sel=i;				// index of var_index in comp_v
				v_aux[i]=x[*j];
				i++;
			}
			df=df-Fl[l][binvec2dec(v_aux,k)];			// removing contribution
			// Flipping the bit
			if (v_aux[i_sel]==0)
				v_aux[i_sel]=1;
			else
				v_aux[i_sel]=0;
			df=df+Fl[l][binvec2dec(v_aux,k)];			// adding contribution	
		}
		
		delete [] v_aux;
	}

	return df;
}


/******************************************************************************\
*					Local Search: First Improvement							   *
\******************************************************************************/
double Mk::LsFi (int *x, double f){
	int i, j=0, *order_v, nr=0;
	double df;

	// Randomly defining the order of the search
	order_v=new int[N];	
	randPerm(order_v, N);
	// Local Search		
	while(j<N){
		i=order_v[j];						
		df=dfElement(x,i);		// Computing the difference in the fitness when flipping the i-th bit							
		if (df>EPS2){
			// Flipping the i-th bit
			if (x[i]==0)
				x[i]=1;
			else
				x[i]=0;
			f=f+df/N;				// new fitness	
			j=0;
			nr++;
			if (nr>N)
				j=N;				// a limit for reseting j is here defined 				
		}
		else
			j++;		 	
	}

	
	delete [] order_v;
	
	return f;

}


/******************************************************************************\
*								 Partition Crossover						   *
* Reference: Tinós R, Whitley D, Chicano F. (2015). "Partition crossover for   *
*    pseudo-boolean optimization",  In Proc. of the 2015 ACM Conference on     *
*    Foundations of Genetic Algorithms XIII, pp. 137-149.                      *
\******************************************************************************/
double Mk::PX(int *parent1, int *parent2, int *offspring){
	int i, j, c, n_comp, n_vert=0, var_index, var_index_next, *comp_id, *v_aux, *f_comp, *map_var2vert, *map_vert2var;
	double fit_offspring=0.0, cost_remaining=0.0;									
	list<int>::iterator ii;
	Graph *G_rec;	
	recomb_comp *comp;
	
	comp_id=new int[N];									// indicates the recombining component associated with each variable
	f_comp=new int[M];									// indicates the recombining component associated with each subfunction fl
	map_var2vert=new int[N];							// vector used to map variables to vertices
	map_vert2var=new int[N];							// vector used to map vertices to variables
	v_aux=new int[k];									// auxiliary vector for the values of the variables in a subfunction fl
	
	n_rec_comp=0;

	// Assigning common variables
	for (i=0;i<N;i++){
		if ( parent1[i]!=parent2[i] ){
			map_var2vert[i]=n_vert;
			map_vert2var[n_vert]=i;
			n_vert++;
		}		
		else{
			map_var2vert[i]=-1;
			offspring[i]=parent1[i];
		}			
	}
	// Initializing f_comp
	for (int l=0;l<M;l++)
		f_comp[l]=-1;									// indicates the component associated to fl	(when equal to -1, no component is associated to fl)
	
	// Step 1: Creating the recombination graph and assigning common variables
	G_rec=new Graph(n_vert);				// recombination graph 
	for (i=0;i<N; i++){
		for(list<int>::iterator jj = VIG[i].begin(); jj != VIG[i].end(); jj++){
			j=*jj;
			var_index=map_var2vert[i];
			var_index_next=map_var2vert[j];
			if ( parent1[i]!=parent2[i] && parent1[j]!=parent2[j] && G_rec->isThereEdge(var_index,var_index_next)==0 )		
					G_rec->addUndirectedEdge(var_index,var_index_next);			
		}
	}	
		
	// Step 2: Finding the connected components of G_rec and computing its properties
	//G_rec->printGraph();
	G_rec->connectedComponents(comp_id);
	n_comp=0;											// number of components
	for (i=0;i<n_vert;i++)
		if (comp_id[i]>n_comp)
			n_comp=comp_id[i];	
	n_comp++;	 	// remember that the first component has label 0 									
	comp=new recomb_comp[n_comp];						// recombining components		
	for (i=0;i<n_comp;i++){
		comp[i].size=0;
		comp[i].cost_parent1=0.0;
		comp[i].cost_parent2=0.0;
	}
	for (i=0;i<n_vert;i++)
		comp[ comp_id[i] ].size += 1;
	for (i=0;i<n_vert;i++){
		j=map_vert2var[i];
		c=comp_id[i];
		comp[c].var_id.push_back(j);
		ii=PSI[j].begin();
		while (ii != PSI[j].end()){		
			f_comp[*ii]=c;		// indicates that component c is associated to fl
			ii++;
		}
	}			
		
	// Step 3: Generating offspring			
	for (int l=0;l<M;l++){
		// Add fl to cost of the component - Parent 1		
		i=0;
		for(list<int>::iterator jj = PHI[l].begin(); jj != PHI[l].end(); jj++){	
			v_aux[i]=parent1[*jj];
			i++;
		}
		if (f_comp[l]==-1){
			cost_remaining += Fl[l][binvec2dec(v_aux,k)];
		}
		else{
			comp[ f_comp[l] ].cost_parent1 += Fl[l][binvec2dec(v_aux,k)];	
			// Add fl to cost of the component - Parent 2
			i=0;
			for(list<int>::iterator jj = PHI[l].begin(); jj != PHI[l].end(); jj++){	
				v_aux[i]=parent2[*jj];
				i++;
			}
			comp[ f_comp[l] ].cost_parent2 += Fl[l][binvec2dec(v_aux,k)];	
		}	
	}
	fit_offspring += cost_remaining;
	for (i=0;i<n_comp;i++){
		n_rec_comp++;
		if (comp[i].cost_parent1>comp[i].cost_parent2){
			fit_offspring += comp[i].cost_parent1;
			ii=comp[i].var_id.begin();
			while (ii != comp[i].var_id.end() ){
				offspring[*ii]=parent1[*ii];
				ii++;
			}
		}
		else{
			fit_offspring += comp[i].cost_parent2;
			ii=comp[i].var_id.begin();
			while (ii != comp[i].var_id.end() ){
				offspring[*ii]=parent2[*ii];
				ii++;
			}
		}
	}
		
	fit_offspring/=M;
		
	delete [] comp_id;
	delete [] f_comp;
	delete [] map_var2vert;
	delete [] map_vert2var;
	delete [] v_aux;
	delete [] comp;
	delete G_rec; 

	 	
	return (fit_offspring);
	
}


/******************************************************************************************\
*		Simple epsilon - Partition Crossover	(i.e., without removing common vertices)   *
* See ePX in next function.    															   *
\******************************************************************************************/
double Mk::sePX(int *parent1, int *parent2, int *offspring, double alpha){
	int i, j, c, var_index, i_next, var_index_next, n_comb_fl, n_comp, n_vert=0;
	int *vi_aux, *vx_aux, *vy_aux, *v_aux, *comp_id, *f_comp, *map_var2vert, *map_vert2var, *omegaH, *set_offspring;
	double ml, d_aux, ful, fit_offspring=0.0, cost_remaining1=0.0, cost_remaining2=0.0;
	list<int>::iterator ii;
	Graph *G_rec;	
	recomb_comp *comp;

	comp_id=new int[N];									// indicates the recombining component associated with each variable
	f_comp=new int[M];									// indicates the recombining component associated with each subfunction fl
	omegaH=new int[M];									// indicates that the indice l of the subfunction is in Omega_h (1) or not (0)
	map_var2vert=new int[N];							// vector used to map variables to vertices (-1 indicates that there are no vertices associated to the variable)
	map_vert2var=new int[N];							// vector used to map vertices to variables
	set_offspring=new int[N];							// indicates if element of the offspring was set (1) or not (0)
	v_aux=new int[k];									// (general) auxiliary vector
	vi_aux=new int[k];									// auxiliary vector for the indices of the variables in PHI(l)
	vx_aux=new int[k];									// auxiliary vector for the values of the variables indicated by PHI(l) for parent 1 (x)
	vy_aux=new int[k];									// auxiliary vector for the values of the variables indicated by PHI(l) for parent 2 (y)

	n_rec_comp=0;
	
	// Initializing map_var2vert
	for (i=0;i<N;i++)
		map_var2vert[i]=-1;	
	
	// Step 1: Creating the recombination graph
	for (int l=0;l<M;l++){
		f_comp[l]=-1;									// indicates the component associated to fl	(when equal to -1, no component is associated to fl)
		// Computing ml=max(fl(x),fl(y))
		i=0;
		for(ii = PHI[l].begin(); ii != PHI[l].end(); ii++){	
			vx_aux[i]=parent1[*ii];
			vy_aux[i]=parent2[*ii];
			vi_aux[i]=*ii;
			i++;
		}
		ml=Fl[l][binvec2dec(vx_aux,k)];
		ful=ml;
		d_aux=Fl[l][binvec2dec(vy_aux,k)];
		if (d_aux>ml)
			ml=d_aux;
		else
			ful=d_aux;
		
		// Computing ful=min(fl(x,y)); testing all combinations
		n_comb_fl=(int) pow(2.0,k)-2;   					// Fl is a matrix N x n_comb_fl; remember that 00...0 (for parent 1) and 11...1 (for parent 2) were already examined 
		i=1;
		while ( i<=n_comb_fl && ( fabs(ml-ful)<fabs((1.0-alpha)*ml) ) ){
			dec2binvec(i,v_aux,k);
			for (j=0;j<k;j++){
				if (v_aux[j]==0)
					v_aux[j]=vx_aux[j];	
				else
					v_aux[j]=vy_aux[j];
			}
			d_aux=Fl[l][binvec2dec(v_aux,k)];
			if (d_aux<ful)
				ful=d_aux;	
			i++;			
		}			
	
		// Adding edges to G_rec and updating omegaH
		if ( i>n_comb_fl){
			omegaH[l]=0;
		}
		else {
			omegaH[l]=1;
			// Adding Vertices
			for (i=0;i<k;i++){
				var_index=vi_aux[i];
				if (map_var2vert[var_index]==-1){						
					map_var2vert[var_index]=n_vert;
					map_vert2var[n_vert]=var_index;
					n_vert++;
				}									
				i_next=i+1;
				while (i_next<k){
					var_index_next=vi_aux[i_next];
					if (map_var2vert[var_index_next]==-1){						
						map_var2vert[var_index_next]=n_vert;
						map_vert2var[n_vert]=var_index_next;
						n_vert++;
					}
					i_next++;						
				}
			}
		}			
	}
	G_rec=new Graph(n_vert);									// recombination graph 
	if (n_vert>0){
		for (int l=0;l<M;l++){
			if (omegaH[l]==1){
				i=0;
				for(ii = PHI[l].begin(); ii != PHI[l].end(); ii++){	
					vi_aux[i]=*ii;
					i++;
				}
				// adding edges
				for (i=0;i<k;i++){
					var_index=map_var2vert[vi_aux[i]];
					i_next=i+1;
					while (i_next<k){
						var_index_next=map_var2vert[vi_aux[i_next]];
						if (G_rec->isThereEdge(var_index,var_index_next)==0)
							G_rec->addUndirectedEdge(var_index,var_index_next);						
						i_next++;
					}
				}
			}			
		}
	}
	//G_rec->printGraph();
	
	// Step 2: Finding the connected components of G_rec and computing its properties
	G_rec->connectedComponents(comp_id);
	for (i=0;i<N;i++)
		set_offspring[i]=0;
	n_comp=0;											// number of components
	for (i=0;i<n_vert;i++){
		if (comp_id[i]>n_comp)
			n_comp=comp_id[i];	
	}
	n_comp++;	   										// remember that the first component has label 0
	comp=new recomb_comp[n_comp];						// recombining components		
	for (i=0;i<n_comp;i++){
			comp[i].size=0;
			comp[i].cost_parent1=0.0;
			comp[i].cost_parent2=0.0;
	}
	for (i=0;i<n_vert;i++)
		comp[ comp_id[i] ].size += 1;	
	for (i=0;i<n_vert;i++){
		j=map_vert2var[i];
		c=comp_id[i];
		comp[c].var_id.push_back(j);
		ii=PSI[j].begin();
		while (ii != PSI[j].end()){	
			if(omegaH[*ii]==1)	
				f_comp[*ii]=c;		// indicates that component c is associated to fl
			ii++;
		}
	}	
		
	// Step 3: Generating offspring
	// for homogeneous fl
	for (int l=0;l<M;l++){		
		if (omegaH[l]==1){
			// Add fl to cost of the component - Parent 1		
			i=0;
			for(ii = PHI[l].begin(); ii != PHI[l].end(); ii++){	
				v_aux[i]=parent1[*ii];
				i++;
			}			
			comp[ f_comp[l] ].cost_parent1 += Fl[l][binvec2dec(v_aux,k)];	
			// Add fl to cost of the component - Parent 2
			i=0;
			for(ii = PHI[l].begin(); ii != PHI[l].end(); ii++){	
				v_aux[i]=parent2[*ii];
				i++;
			}
			comp[ f_comp[l] ].cost_parent2 += Fl[l][binvec2dec(v_aux,k)];	
		}
	}
	for (i=0;i<n_comp;i++){
		n_rec_comp++;
		if (comp[i].cost_parent1>comp[i].cost_parent2){
			fit_offspring += comp[i].cost_parent1;
			ii=comp[i].var_id.begin();
			while (ii != comp[i].var_id.end() ){
				offspring[*ii]=parent1[*ii];
				set_offspring[*ii]=1;
				ii++;
			}
		}
		else{
			fit_offspring += comp[i].cost_parent2;
			ii=comp[i].var_id.begin();
			while (ii != comp[i].var_id.end() ){
				offspring[*ii]=parent2[*ii];
				set_offspring[*ii]=1;
				ii++;
			}
		}
	}
	// for heteregenous fl
	for (int l=0;l<M;l++){
		if (omegaH[l]==0){			
			// Add fl to cost of the remaining component - Parent 1		
			i=0;
			for(ii = PHI[l].begin(); ii != PHI[l].end(); ii++){	
				if ( set_offspring[*ii]==0 ) {
					v_aux[i]=parent1[*ii];						
				}
				else
					v_aux[i]=offspring[*ii];
				i++;
			}			
			cost_remaining1 = Fl[l][binvec2dec(v_aux,k)];
			// Add fl to cost of the remaining component - Parent 1		
			i=0;
			for(ii = PHI[l].begin(); ii != PHI[l].end(); ii++){	
				if ( set_offspring[*ii]==0 ) {
					v_aux[i]=parent2[*ii];
				}
				else
					v_aux[i]=offspring[*ii];
				
				i++;
			}			
			cost_remaining2 = Fl[l][binvec2dec(v_aux,k)];
			// Assigning partial solutions (in fl) for best parent
			if (cost_remaining1>=cost_remaining2){
				fit_offspring += cost_remaining1;
				i=0;
				for(ii = PHI[l].begin(); ii != PHI[l].end(); ii++){	
					if ( set_offspring[*ii]==0 ) {
						offspring[*ii]=parent1[*ii];
						set_offspring[*ii]=1;
					}					
					i++;
				}	
			}
			else{
				fit_offspring += cost_remaining2;
				i=0;
				for(ii = PHI[l].begin(); ii != PHI[l].end(); ii++){	
					if ( set_offspring[*ii]==0 ) {
						offspring[*ii]=parent2[*ii];
						set_offspring[*ii]=1;
					}					
					i++;
				}	
			}						
		}
	}
	
	fit_offspring/=M;
	
	delete [] comp_id;
	delete [] f_comp;
	delete [] omegaH;
	delete [] map_var2vert;
	delete [] map_vert2var;
	delete [] set_offspring;
	delete [] v_aux;
	delete [] vi_aux;
	delete [] vx_aux;
	delete [] vy_aux;
	delete [] comp;
	delete G_rec;  
	 	
	return (fit_offspring);
}


/***************************************************************************\
*					epsilon - Partition Crossover				   		   	*
* 	Reference: Tinos, R.; Whitley, D.; Chicano, F. & Ochoa, G. (2021), 		*                    
* 				"Partition Crossover for Continuous Optimization: ePX",		*
*				Proc. of GECCO'2021.										*
\***************************************************************************/
double Mk::ePX(int *parent1, int *parent2, int *offspring, double alpha){
	int i, j, c, var_index, i_next, var_index_next, n_comb_fl, n_comp, n_elements_fl, n_vert=0;
	int *vi_aux, *vx_aux, *vy_aux, *v_aux, *comp_id, *f_comp, *map_var2vert, *map_vert2var, *omegaH, *set_offspring, *parent_equal_vector, *n_elem_equal_in_fl;
	double ml, d_aux, ful, fit_offspring=0.0, cost_remaining1=0.0, cost_remaining2=0.0;
	list<int>::iterator ii;
	Graph *G_rec;	
	recomb_comp *comp;

	comp_id=new int[N];									// indicates the recombining component associated with each variable
	f_comp=new int[M];									// indicates the recombining component associated with each subfunction fl
	omegaH=new int[M];									// indicates that the indice l of the subfunction is in Omega_h (1) or not (0)
	map_var2vert=new int[N];							// vector used to map variables to vertices (-1 indicates that there are no vertices associated to the variable)
	map_vert2var=new int[N];							// vector used to map vertices to variables
	n_elem_equal_in_fl=new int[M];						// indicates the number of elements in fl that are equal in the parents
	parent_equal_vector=new int[N];						// indicates the elements that are equal in the parents
	set_offspring=new int[N];							// indicates if element of the offspring was set (1) or not (0)
	v_aux=new int[k];									// (general) auxiliary vector
	vi_aux=new int[k];									// auxiliary vector for the indices of the variables in PHI(l)
	vx_aux=new int[k];									// auxiliary vector for the values of the variables indicated by PHI(l) for parent 1 (x)
	vy_aux=new int[k];									// auxiliary vector for the values of the variables indicated by PHI(l) for parent 2 (y)

	n_rec_comp=0;
	
	// Initializing  map_var2vert and assigning common variables
	for (i=0;i<N;i++){
		map_var2vert[i]=-1;	
		if (parent1[i]==parent2[i]){
			parent_equal_vector[i]=1;
			offspring[i]=parent1[i];
			set_offspring[i]=1;
		}
		else{
			parent_equal_vector[i]=0;	
			set_offspring[i]=0;
		}
	}

	// Step 1: Creating the recombination graph
	for (int l=0;l<M;l++){
		f_comp[l]=-1;									// indicates the component associated to fl	(when equal to -1, no component is associated to fl)
		// Computing ml=max(fl(x),fl(y))
		n_elem_equal_in_fl[l]=0;							
		i=0;
		for(ii = PHI[l].begin(); ii != PHI[l].end(); ii++){	
			vx_aux[i]=parent1[*ii];
			vy_aux[i]=parent2[*ii];
			if (vx_aux[i]==vy_aux[i])
				n_elem_equal_in_fl[l]++;		
			vi_aux[i]=*ii;
			i++;
		}
		n_elements_fl=i;
		if (n_elem_equal_in_fl[l]<n_elements_fl){
			// partial solutions for fl are different 
			ml=Fl[l][binvec2dec(vx_aux,k)];
			ful=ml;
			d_aux=Fl[l][binvec2dec(vy_aux,k)];
			if (d_aux>ml)
				ml=d_aux;
			else
				ful=d_aux;
		
			// Computing ful=min(fl(x,y)); testing all combinations
			n_comb_fl=(int) pow(2.0,k)-2;   					// Fl is a matrix N x n_comb_fl; remember that 00...0 (for parent 1) and 11...1 (for parent 2) were already examined 
			i=1;
			while ( i<=n_comb_fl && ( fabs(ml-ful)<fabs((1.0-alpha)*ml) ) ){
				dec2binvec(i,v_aux,k);
				for (j=0;j<k;j++){
					if (v_aux[j]==0)
						v_aux[j]=vx_aux[j];	
					else
						v_aux[j]=vy_aux[j];
				}
				d_aux=Fl[l][binvec2dec(v_aux,k)];
				if (d_aux<ful)
					ful=d_aux;	
				i++;			
			}		
		}
			// Adding edges to G_rec and updating omegaH
		if (  n_elem_equal_in_fl[l]==n_elements_fl || i>n_comb_fl ){
			omegaH[l]=0;
		}
		else {
			omegaH[l]=1;
			// Adding Vertices
			for (i=0;i<k;i++){
				var_index=vi_aux[i];
				if (map_var2vert[var_index]==-1 && parent_equal_vector[var_index]==0){						
					map_var2vert[var_index]=n_vert;
					map_vert2var[n_vert]=var_index;
					n_vert++;
				}									
				i_next=i+1;
				while (i_next<k){
					var_index_next=vi_aux[i_next];
					if (map_var2vert[var_index_next]==-1 && parent_equal_vector[var_index_next]==0){						
						map_var2vert[var_index_next]=n_vert;
						map_vert2var[n_vert]=var_index_next;
						n_vert++;
					}
					i_next++;						
				}
			}
		}			
	}
	G_rec=new Graph(n_vert);									// recombination graph 
	if (n_vert>0){
		for (int l=0;l<M;l++){
			if (omegaH[l]==1){
				i=0;
				for(ii = PHI[l].begin(); ii != PHI[l].end(); ii++){	
					vi_aux[i]=*ii;
					i++;
				}
				// adding edges
				for (i=0;i<k;i++){
					var_index=map_var2vert[vi_aux[i]];
					i_next=i+1;
					while (i_next<k){
						var_index_next=map_var2vert[vi_aux[i_next]];
						if ( (parent_equal_vector[vi_aux[i]]==0 && parent_equal_vector[vi_aux[i_next]]==0) && G_rec->isThereEdge(var_index,var_index_next)==0)
							G_rec->addUndirectedEdge(var_index,var_index_next);						
						i_next++;
					}
				}
			}			
		}
	}
	//G_rec->printGraph();
	
	// Step 2: Finding the connected components of G_rec and computing its properties
	G_rec->connectedComponents(comp_id);
	n_comp=0;											// number of components
	for (i=0;i<n_vert;i++)
		if (comp_id[i]>n_comp)
			n_comp=comp_id[i];	
	n_comp++;	   										// remember that the first component has label 0
	comp=new recomb_comp[n_comp];						// recombining components		
	for (i=0;i<n_comp;i++){
			comp[i].size=0;
			comp[i].cost_parent1=0.0;
			comp[i].cost_parent2=0.0;
	}
	for (i=0;i<n_vert;i++)
		comp[ comp_id[i] ].size += 1;
	for (i=0;i<n_vert;i++){
		j=map_vert2var[i];
		c=comp_id[i];
		comp[c].var_id.push_back(j);
		ii=PSI[j].begin();
		while (ii != PSI[j].end()){	
			if(omegaH[*ii]==1)	
				f_comp[*ii]=c;							// indicates that component c is associated to fl
			ii++;
		}
	}

 	// Step 3: Generating offspring
	// for homogeneous fl
	for (int l=0;l<M;l++){		
		if (omegaH[l]==1){
			// Add fl to cost of the component - Parent 1		
			i=0;
			for(ii = PHI[l].begin(); ii != PHI[l].end(); ii++){	
				v_aux[i]=parent1[*ii];
				i++;
			}			
			comp[ f_comp[l] ].cost_parent1 += Fl[l][binvec2dec(v_aux,k)];	
			// Add fl to cost of the component - Parent 2
			i=0;
			for(ii = PHI[l].begin(); ii != PHI[l].end(); ii++){	
				v_aux[i]=parent2[*ii];
				i++;
			}
			comp[ f_comp[l] ].cost_parent2 += Fl[l][binvec2dec(v_aux,k)];	
		}
	}
	for (i=0;i<n_comp;i++){
		n_rec_comp++;
		if (comp[i].cost_parent1>comp[i].cost_parent2){
			fit_offspring += comp[i].cost_parent1;
			ii=comp[i].var_id.begin();
			while (ii != comp[i].var_id.end() ){
				offspring[*ii]=parent1[*ii];
				set_offspring[*ii]=1;
				ii++;
			}
		}
		else{
			fit_offspring += comp[i].cost_parent2;
			ii=comp[i].var_id.begin();
			while (ii != comp[i].var_id.end() ){
				offspring[*ii]=parent2[*ii];
				set_offspring[*ii]=1;
				ii++;
			}
		}
	}
	// for heteregenous fl
	for (int l=0;l<M;l++){
		if (omegaH[l]==0){			
			// Add fl to cost of the remaining component - Parent 1		
			i=0;
			for(ii = PHI[l].begin(); ii != PHI[l].end(); ii++){	
				if ( set_offspring[*ii]==0 ) {
					v_aux[i]=parent1[*ii];						
				}
				else
					v_aux[i]=offspring[*ii];
				i++;
			}			
			cost_remaining1 = Fl[l][binvec2dec(v_aux,k)];
			// Add fl to cost of the remaining component - Parent 1		
			i=0;
			for(ii = PHI[l].begin(); ii != PHI[l].end(); ii++){	
				if ( set_offspring[*ii]==0 ) {
					v_aux[i]=parent2[*ii];
				}
				else
					v_aux[i]=offspring[*ii];
				
				i++;
			}			
			cost_remaining2 = Fl[l][binvec2dec(v_aux,k)];
			// Assigning partial solutions (in fl) for best parent
			if (cost_remaining1>=cost_remaining2){
				fit_offspring += cost_remaining1;
				i=0;
				for(ii = PHI[l].begin(); ii != PHI[l].end(); ii++){	
					if ( set_offspring[*ii]==0 ) {
						offspring[*ii]=parent1[*ii];
						set_offspring[*ii]=1;
					}					
					i++;
				}	
			}
			else{
				fit_offspring += cost_remaining2;
				i=0;
				for(ii = PHI[l].begin(); ii != PHI[l].end(); ii++){	
					if ( set_offspring[*ii]==0 ) {
						offspring[*ii]=parent2[*ii];
						set_offspring[*ii]=1;
					}					
					i++;
				}	
			}						
		}
	}
	
	fit_offspring/=M;

	delete [] comp_id;
	delete [] f_comp;
	delete [] omegaH;
	delete [] map_var2vert;
	delete [] map_vert2var;
	delete [] n_elem_equal_in_fl;
	delete [] parent_equal_vector;
	delete [] set_offspring;
	delete [] v_aux;
	delete [] vi_aux;
	delete [] vx_aux;
	delete [] vy_aux;
	delete [] comp;
	delete G_rec;  
	
	return (fit_offspring);
}


/******************************************************************************\
*								Print Mk information						   *
\******************************************************************************/
void Mk::print(void){
	int aux;
	
	cout<< "Mk Landscape: "<<endl<<"  M="<<M<<endl<<"  k= "<< k<<endl<<"  N= "<< N<<endl;

	aux=(int) pow(2.0,k);
	
   cout<<"Variables: subfunctions "<<endl;
	for (int i=0;i<N;i++){
		cout<<" x_"<<i<<": ";
		for(list<int>::iterator j = PSI[i].begin(); j != PSI[i].end(); j++) 
			cout<<"  f_"<<*j;
		cout<<endl;
	}
	
	cout<<"VIG: variables "<<endl;
	for (int i=0;i<N;i++){
		cout<<" x_"<<i<<": ";
		for(list<int>::iterator j = VIG[i].begin(); j != VIG[i].end(); j++) 
			cout<<"  x_"<<*j;
		cout<<endl;
	}
	
	cout<<"Subfuncions: variables "<<endl;
	for (int l=0;l<M;l++){
		cout<<" f_"<<l<<": ";
		for(list<int>::iterator j = PHI[l].begin(); j != PHI[l].end(); j++) 
			cout<<"  x_"<<*j;
		cout<<endl;
	}
	
	cout<<"Subfuncions: contributions "<<endl;
	for (int l=0;l<M;l++){
		cout<<" f_"<<l<<":";
		for (int j=0;j<aux;j++)
			cout<<"  "<<Fl[l][j];
		cout<<endl;
	}
	
}
