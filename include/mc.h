/**************************************************/
/*  This code was made to simulate Lennard-Jones  */
/*  polymers with Monte Carlo simulations.        */
/**************************************************/

#define NDIM 3
#define pi 3.1415926535897932384626433

void initial_pos(void);
void surfinitial_pos(void);
void overlap(int, int*);
void read_in_conf(FILE*);
void gen_suc_monomer(int, int);
void gen_suc_chains(int, int, double);
void gen_suc_linback(int, int);
void dyn_alloc_array(void);
void write_conf(FILE*);
void write_pos(int, FILE*);
void CalcEnergy(int, double*);
void CalcAccProb (double, double, double*);
void Calculations(void);
void NeighborList(void);
void mcmove(void);
void Setup(void);