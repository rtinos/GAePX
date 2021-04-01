# GAePX
Genetic Algorithm with Epsilon - Partition Crossover (ePX)

*** Genetic Algorithm with epsilon-Partition Crossover (ePX) for NK Landscapes ***

Description: This is the source code for the the GA with ePX for NK landscapes. 

Reference:  Tinos, R.; Whitley, D.; Chicano, F. & Ochoa, G. (2021), "Partition Crossover for Continuous Optimization: ePX", Proc. of GECCO'2021.	

Contact: Renato Tinos <rtinos@ffclrp.usp.br>


Running the code: ./ga_epx_nkland <N> <K> <model_NK> <crossover_type>

<N>: size of the NK landscapes instance. N is a positive integer.

<K>: controls the epistasis degree of the NK landscapes instance. K is a non-negative integer smaller than or equal to N.

<model_NK>: NK landscapes model. 0-random NK; 1-regular NK (in this case, a table with values of subfuncions f_l for instances solved by dynamic programming is loaded); 2-random NKq with q=2^(K+1)/2-1; ; 3-random NKq with q=2^(K+1)-1.

<crossover_type>: type of crossover. 1-2PointX; 2-UX; 3-PX, 4-simple ePX (without removing common vertices), 5-ePX, 6-Multiple 2PointX, 7-Multiple UX.	


Example for running the code for: random NK landscapes with N=100, K=2, crossover ePX

make

./ga_epx_nkland 100 2 0 5


Observation 1: Class Mk, given in Mk.h, implements Mk landscapes (that is a generalization of NK landscapes). Some functions in Mk.h:

- double Mk::ePX(int *parent1, int *parent2, int *offspring, double alpha);  epsilon-Partition Crossover. Inputs: two parents, pointer to offspring, and alpha=1-epsilon. Outputs: offspring and its fitness.
	
- double Mk::PX(int *parent1, int *parent2, int *offspring);  Partition Crossover. Inputs: two parents, pointer to offspring. Outputs: offspring and its fitness.

- double Mk::LsFi (int *x, double f); First-Improvement Local Search. Inputs: parent x ant its fitness f; Outputs: optimized solution x and its fitness.

- double Mk::compFitness (int *x); Fitness computation for Mk landscapes. Input: solution x; Output: fitness of x.

		
Observation 2: file global.cpp contains the parameters of the GA (examples: number of runs, population size, and crossover rate) and the problem (examples: number of instances and runs).


Observation 3: ga_epx_nkland generates 6 main files
 
- bfi_N%d_K%d_p%d_c%d.dat (example for ./ga_epx_nkland 100 2 0 5: ,bfi_N100_K2_p0_c5.dat): best fitness found in each run
	
- bind_N%d_K%d_p%d_c%d.dat: best individuals found in each run.

- time_N%d_K%d_p%d_c%d.dat: time for each run.

- sucRate_N%d_K%d_p%d_c%d.dat: successful recombination rate for each run.

- impRate_N%d_K%d_p%d_c%d.dat: improvement rate for each run.

- recComp_N%d_K%d_p%d_c%d.dat: saves the number of recombination components for PX, ePX and sePX. 
