/**************************************************/
/*  This code was made to simulate Lennard-Jones  */
/*  polymers with Monte Carlo simulations.        */
/**************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <rng.h>
#include <mc.h>

static double box_l = 100.0;
static int N_poly = 1, N_back = 128, N_side = 0, N_chains = 0;
static int poly_N;//number of monomers in a chain
static int N; //total number of monomrs
static int N_bond; //number of bonds
static int N_neighs = 20; //max number of neighbors per atom
static int R_neighs = 2.; //neighbor shell radius
static int t_neighs = 100000; //number of steps before rebuilding the list
static long long i_step, tot_step, eq_step;
static double delta = 0.06; /* monte carlo displacement */
static double **r;
static int *b1, *b2;
static int **neighs;
/* Lennard-Jones parameters */
static double LJ_eps = 1.0, LJ_sig = 1.0, LJ_cut = 1.122462048309373, K_fene = 30.0, R_fene = 1.5;
static double potshift;
/* Surface 9-3 potential parameters */
static const double VDW_epssurf = 5.0, VDW_sigsurf = 1.0; 
/* Log potential parameters */
static double K_logpull = 3.0, R_logpull = 1.0;
/* Simulation parameters */
const int CHAINS = 0; // Determines if chains will be attached.
const int SURF = 0; // Determines if there is a surface.
const int BONDS = 1; // Determines if whether there are FENE bonds.
const int LJ = 1;  // Determines whether to use the Lennard-Jones potential.
const int PARAB = 0; // Parabolic potential.
const int LOGPULL = 0; // Logarithmic potential between the ends of two chains
static double acptcnt, cnt;
static double rrmax, rrmin, rravg;
static double bondMax, bondMin, bondTens, chainR2;

static time_t wallClockTime0;

FILE *tapePtr, *binPtr;

/********************************************************/
/* Algegraic Functions */
double dsqr(double x) { return x * x; }
double min(double a, double b) { return (a < b ? a : b);};
double max(double a, double b) { return (a > b ? a : b);};

/********************************************************/
/* Initial Configuration Function */
void initial_pos(void)
{
    int h,i,j,k,i_back,i_side,acc_flag,bond_counter; 
    FILE *posPtr;
    char pos_file[100];
    /*open data file*/
    sprintf(pos_file, "initialpos.dat");
    if ((posPtr =fopen (pos_file, "w")) == NULL) {
        fprintf(stderr, "file %s could not be opened!", pos_file);
        exit(1);
    }
    
    /* Repeat the initial particle generation */
    bond_counter = 0;
    for(h = 0; h < N_poly; h++){    
        /* Generate the first particle in each polymer */
        do{
            acc_flag = 0;
            for(k = 0; k < 3; k++) {
                r[h * poly_N][k] = (Randomr()-0.5)*box_l;
            }
            if(h != 0) overlap(h * poly_N, &acc_flag);
        } while(acc_flag == 1);
        
        /* print the initial positions */
        //printf("%5d%15.6le%15.6le%15.6le\n", h*poly_N, r[h * poly_N][0], r[h * poly_N][1], r[h * poly_N][2]);  
        fprintf(posPtr, "%15.6le%15.6le%15.6le\n", r[h * poly_N][0], r[h * poly_N][1], r[h * poly_N][2]);
        /* generating the remaining N-1 particles in a backbone*/
        for(i = 1; i < N_back; i++){     
            i_back = h * poly_N + i;
            gen_suc_monomer(i_back, i_back-1);
            fprintf(posPtr, "%15.6le%15.6le%15.6le\n", r[i_back][0], r[i_back][1], r[i_back][2]);
            //printf("%5d%15.6le%15.6le%15.6le\n",i_back, r[i_back][0], r[i_back][1], r[i_back][2]);
            b1[bond_counter] = i_back - 1;
            b2[bond_counter] = i_back;
            bond_counter++;
        }
        if (CHAINS) {
            /* Generate brush chains */    
            for(i=0; i < N_back; i++){   
                i_back = h * poly_N + i;
                /* Generate the first monomer on the side chain */
                i_side = h * poly_N + N_back + i * N_side;
                gen_suc_monomer(i_side, i_back);
                fprintf(posPtr, "%15.6le%15.6le%15.6le\n", r[i_side][0], r[i_side][1], r[i_side][2]);
                //printf("%5d%5d%15.6le%15.6le%15.6le\n",i_back,i_side,r[i_side][0], r[i_side][1], r[i_side][2]);
                b1[bond_counter] = i_back;
                b2[bond_counter] = i_side;
                bond_counter++;
                /* Generate succesive monomers on the side chain */
                for(j=1; j < N_side; j++){
                    i_side = h * poly_N + N_back + i * N_side + j;    
                    gen_suc_monomer(i_side, i_side-1);
                    fprintf(posPtr, "%15.6le%15.6le%15.6le\n", r[i_side][0], r[i_side][1], r[i_side][2]);
                    //printf("%5d%15.6le%15.6le%15.6le\n",i_side,r[i_side][0], r[i_side][1], r[i_side][2]);
                    b1[bond_counter] = i_side - 1;
                    b2[bond_counter] = i_side;
                    bond_counter++;
                }
            }
        }
    }
    fclose(posPtr);
    if (0) {
        printf("---- Bonds ----\n");
        for(i=0;i<N_bond;i++) {printf("%5d%5d\n",b1[i],b2[i]);}
    }
}

/************************************************************/
void gen_suc_monomer(int i, int j)
{
    int acc_flag;
    double theta, phi;
    
    do{
        acc_flag = 0;
        theta = Randomr()*2.*pi;
        phi = Randomr()*pi;
        r[i][0] = r[j][0] + (0.97 * cos(phi) * cos(theta));
        r[i][1] = r[j][1] + (0.97 * cos(phi) * sin(theta));
        r[i][2] = r[j][2] + (0.97 * sin(phi));
        /* To test whether the ith particles is overlapping with all
         the previously generated particles */ 
        overlap(i, &acc_flag);
    } while(acc_flag == 1);  
}

/*************************************************************/
/* Initial Configuration Function Surface */
void surf_initial_pos(void)
{
    int h,i,j,k,i_back,i_side,bond_counter;
    double theta; 
    FILE *posPtr;
    char pos_file[100];
    /*open data file*/
    sprintf(pos_file, "surfinitialpos.dat");
    if ((posPtr =fopen (pos_file, "w")) == NULL) {
        fprintf(stderr, "file %s could not be opened!", pos_file);
        exit(1);
    }
    
    /* Repeat the initial particle generation */
    bond_counter = 0;
    for(h = 0; h < N_poly; h++){    
        /* Generate the first particle in each polymer */
        for(k = 0; k < 3; k++) {
            r[h * poly_N][0] = 0;
            r[h * poly_N][1] = 0;
            r[h * poly_N][2] = (N_side * 1.07) + 1;
        }
        
        /* print the initial positions */
        //printf("%5d%15.6le%15.6le%15.6le\n", h*poly_N, r[h * poly_N][0], r[h * poly_N][1], r[h * poly_N][2]);  
        fprintf(posPtr, "%15.6le%15.6le%15.6le\n", r[h * poly_N][0], r[h * poly_N][1], r[h * poly_N][2]);
        /* generating the remaining N-1 particles in a backbone*/
        for(i = 1; i < N_back; i++){     
            i_back = h * poly_N + i;
            gen_suc_linback(i_back, i_back-1);
            fprintf(posPtr, "%15.6le%15.6le%15.6le\n", r[i_back][0], r[i_back][1], r[i_back][2]);
            //printf("%5d%15.6le%15.6le%15.6le\n",i_back, r[i_back][0], r[i_back][1], r[i_back][2]);
            b1[bond_counter] = i_back - 1;
            b2[bond_counter] = i_back;
            bond_counter++;
        }
        if (CHAINS) {
            /* Generate brush chains */
            for(i=0; i < N_back; i++){   
                for(k=0; k < N_chains; k++){ 
                    i_back = h * poly_N + i;
                    theta = ((2*pi) / N_chains) * k;
                    /* Generate the first monomer on the side chain */
                    i_side = (h * poly_N) + N_back + (N_chains *  N_side * i) + (k * N_side);
                    
                    r[i_side][0] = r[i_back][0];
                    r[i_side][1] = r[i_back][1] + (1.07 * cos(theta));
                    r[i_side][2] = r[i_back][2] + (1.07 * sin(theta));
                    
                    fprintf(posPtr, "%15.6le%15.6le%15.6le\n", r[i_side][0], r[i_side][1], r[i_side][2]);
                    //printf("%5d%5d%15.6le%15.6le%15.6le\n",i_back,i_side,r[i_side][0], r[i_side][1], r[i_side][2]);
                    b1[bond_counter] = i_back;
                    b2[bond_counter] = i_side;
                    bond_counter++;
                    /* Generate succesive monomers on the side chain */
                    for(j=1; j < N_side; j++){
                        i_side = (h * poly_N) + N_back + (N_chains * N_side * i) + (k * N_side) + j;   
                        r[i_side][0] = r[i_side-1][0];
                        r[i_side][1] = r[i_side-1][1] + (1.07 * cos(theta));
                        r[i_side][2] = r[i_side-1][2] + (1.07 * sin(theta));
                        fprintf(posPtr, "%15.6le%15.6le%15.6le\n", r[i_side][0], r[i_side][1], r[i_side][2]);
                        //printf("%5d%15.6le%15.6le%15.6le\n",i_side,r[i_side][0], r[i_side][1], r[i_side][2]);
                        b1[bond_counter] = i_side - 1;
                        b2[bond_counter] = i_side;
                        bond_counter++;
                    }
                }
            }
        }
    }
    fclose(posPtr);
    //for(i=0;i<N_bond;i++) printf("%5d%5d\n",b1[i],b2[i]); 
}

/**********************************************************************/
void gen_suc_linback(int i, int j)
{
    r[i][0] = r[j][0] + 1.07;
    r[i][1] = r[j][1];
    r[i][2] = r[j][2];
}

/***********************************************************/
void gen_suc_chains(int i, int j, double theta)
{
    r[i][0] = r[j][0] + (1.07 * cos(theta));
    r[i][1] = r[j][1] + (1.07 * sin(theta));
    r[i][2] = r[j][2] + (1.07); 
}

/************************************************************/
/* judge the overlap of the particle being generated */
void overlap(int mon_n, int *overlap_flag)
{
    int j,k;
    double d, r2;
    
    for(j = 0; j < mon_n; j ++){
        r2 = 0.0;
        for(k = 0; k < 3; k++){
            d = r[mon_n][k] - r[j][k];
            //d -= floor((d/box_l)+0.5)*box_l;
            //correct for particles that "loop" back around          
            r2 += d*d;
        }
        if (r2 < 0.92) *overlap_flag = 1;
    }
}

/*************************************************************/
/* dynamically allocate the arrays for particles*/
void dyn_alloc_array(void)
{
    int i;
    
    r= malloc(N*sizeof(double*));
    for (i=0; i<N; i++) r[i] = malloc(NDIM*sizeof(double));
    b1= malloc(N_bond*sizeof(int));
    b2= malloc(N_bond*sizeof(int));
    neighs = (int **) malloc(N*sizeof(int *));
    for (i=0;i<N;i++) {
        neighs[i] = (int *) malloc(N_neighs*sizeof(int));
    }
}

/*****************************************************************/
/* Read in the saved configurations */
void read_in_conf(FILE *inPtr)
{
    int i,j;
    
    fscanf(inPtr,"%lld\n",&i_step);
    for(i = 0; i < N; i++){
        for(j = 0; j <3; j++)
            fscanf(inPtr, "%le",&r[i][j]);;
        fscanf(inPtr,"\n");
    }
}
/******************************************************************/
void NeighborList(void)
{
    int i,j,k,cnt[N];
    double d[3],r2,r2neighs;

    r2neighs = R_neighs * R_neighs;
    for (i=0; i<N; i++) {
        cnt[i]=0;
        for (j=0; j<N_neighs; j++) {
            neighs[i][j] = -1;
        }
    }
    for (i=0; i<N-1; i++) {
        for (j=i+1; j<N; j++) {
            r2=0.;
            for (k=0; k<3; k++) {
                d[k] = (r[i][k] - r[j][k]);
                r2 += d[k] * d[k];
            }
            if (r2 <= r2neighs) {
                neighs[i][cnt[i]] = j;
                cnt[i]++;
                neighs[j][cnt[j]] = i;
                cnt[j]++;
            }
            if ( (cnt[i] > N_neighs) || (cnt[j] > N_neighs) ) {
                printf("Error: Number of allowed neighbors is too small!\n");
            }
        }
    }
}

/******************************************************************/
void CalcEnergy(int M, double *UTOT)
{
    int i, k;
    double d[3], r1, r2, r3, r6;
    double ljcut2, ljsig2, ljconst;
    double rfene2, forceConst;
    //double u1,u2;
    
    ljcut2 = LJ_cut * LJ_cut;
    ljsig2 = LJ_sig * LJ_sig;
    *UTOT = 0.;
    
    if (PARAB) {
        for (i = 0; i < N; i++) {
            r2 = 0.0;
            if (i != M) {
                for (k = 0; k < 3; k++) {
                    d[k] = r[M][k] - r[i][k];
                    //if (k < 2) d[k] -= floor((d[k] / box_l) + 0.5) * box_l;
                    //d[k] -= floor((d[k] / box_l) + 0.5) * box_l;
                    r2 += d[k] * d[k];
                }
                //x = sqrtl(r2);
                *UTOT += 50 * r2;
            }
        }
    }
    
    /* Lennard-Jones Potential */
    if (LJ) {
        ljconst = ljcut2 * ljsig2;
        forceConst = 4.0 * LJ_eps;
        //u1 = u2 = 0.;
        for (i = 0; i < N_neighs; i++) {
            if (neighs[M][i] != -1) {
                r2 = 0.;
                for (k = 0; k < 3; k++) {
                    d[k] = r[M][k] - r[neighs[M][i]][k];
                    r2 += d[k] * d[k];
                }
                if (r2 <= ljconst) {
                    r2 = ljsig2 / r2;
                    r6 = r2 * r2 * r2;
                    //u1 += (forceConst * r6 * (r6 - 1.));
                    //u1 -= potshift;
                    *UTOT += (forceConst * r6 * (r6 - 1.));
                    *UTOT -= potshift;
                }
            } else {
                break;
            }
        }

        /* TEST */
        /*
        for (i = 0; i < N; i++) {
            r2 = 0.;
            if (i != M) {
                for (k = 0; k < 3; k++) {
                    d[k] = r[M][k] - r[i][k];
                    r2 += d[k] * d[k];
                }
                if (r2 <= ljconst) {
                    r2 = ljsig2 / r2;
                    r6 = r2 * r2 * r2;
                    u2 += (forceConst * r6 * (r6 - 1.));
                    u2 -= potshift;
                }
            }
        }
        if (u1-u2 != 0. ) printf("%f\n", u1-u2); */
    }
    
    if (BONDS) {
        /* FENE Potential */
        /* Only for linear chain! */
        rfene2 = R_fene * R_fene;
        r2 = 0.;
        if (M != N_back-1) {
            for (k = 0; k < 3; k++) {
                d[k] = r[b1[M]][k] - r[b2[M]][k];
                //d[k] -= floor((d[k] / box_l) + 0.5) * box_l;
                r2 += d[k] * d[k];
            }
            r2 /= ljsig2;
            if (rfene2 < r2) printf("%lld. Broken bond rr = %g (%d:%d)\n", i_step, r2, b1[M], b2[M]);
            if (rfene2 <= r2) {
                *UTOT += 10000.;
            }
            else {
                *UTOT -= 0.5 * K_fene * rfene2 * log1p(-r2/rfene2);
            }
        }
        r2 = 0.;
        if (M != 0) {
            for (k = 0; k < 3; k++) {
                d[k] = r[b1[M-1]][k] - r[b2[M-1]][k];
                //d[k] -= floor((d[k] / box_l) + 0.5) * box_l;
                r2 += dsqr(d[k]);
            }
            r2 /= ljsig2;
            if (rfene2 < r2) printf("%lld. Broken bond rr = %g (%d:%d)\n", i_step, r2, b2[M-1], b1[M-1]);
            if (rfene2 <= r2) {
                *UTOT += 1000.;
            }
            else {
                *UTOT -= 0.5 * K_fene * rfene2 * log1p(-r2/rfene2);
            }
        }
    }
    
    if (SURF) {
        r3 = r[M][2] * r[M][2] * r[M][2];
        r6 = r3 * r3;
        *UTOT += 1.5 * VDW_epssurf * (VDW_sigsurf / r3) * ((VDW_sigsurf / 3.) * (1. / r6) - 1.);
    }
    
    if (LOGPULL) {
        if (M == N_back-1) {
            r2 = 0.;
            for (k = 0; k < 3; k++) {
                d[k] = r[M][k] - r[(M)+1][k];
                //d[k] -= floor((d[k] / box_l) + 0.5) * box_l;
                r2 += d[k] * d[k];
            }
            r1 = sqrt(r2);
            r1 /= R_logpull;
            if (r1 <= 1.) {
                *UTOT += K_logpull * (r1 - 1.);
            } else {
                *UTOT += K_logpull * log(r1);
            }
        }
        if (M == N_back) {
            r2 = 0.;
            for (k = 0; k < 3; k++) {
                d[k] = r[M][k] - r[M-1][k];
                //d[k] -= floor((d[k] / box_l) + 0.5) * box_l;
                r2 += d[k] * d[k];
            }
            r1 = sqrt(r2);
            r1 /= R_logpull;
            if (r1 <= 1.) {
                *UTOT += K_logpull * (r1 - 1.);
            } else {
                *UTOT += K_logpull * log(r1);
            }
        }
    }
}

void VRand (double *p)
{
    double s, x, y;
    
    s = 2.;
    while (s > 1.) {
        x = 2. * Randomr() - 1.;
        y = 2. * Randomr() - 1.;
        s = dsqr(x) + dsqr(y);
    }
    p[2] = 1. - 2. * s;
    s = 2. * sqrt(1. - s);
    p[0] = s * x;
    p[1] = s * y;
}

void CalcAccProb (double oldU, double newU, double *acc)
{
    double B, C;
    
    B = -(newU - oldU);
    C = exp(B);
    *acc = min(1,C);
}

/*****************************************************************/
void mcmove(void)
/* Metropolis Algorithm */
{
    int m, k;
    double acc;
    double p[3], r0[3];
    double U1, U2;
    
    m = (int)(Randomr() * N);
    CalcEnergy(m,&U1);
    VRand(p);
    for (k = 0; k < 3; k++) {
        r0[k] = r[m][k];
        r[m][k] += p[k] * delta;
    }
    CalcEnergy(m,&U2);
    
    //printf("%g\n",U2-U1);
    CalcAccProb(U1,U2,&acc);
    if (Randomr() < acc) {
        acptcnt += 1.;
    } else {
        for (k = 0; k < 3; k++) r[m][k] = r0[k];
    }
}

/****************************************************************/
void Setup(void)
{
    char bin_file[100], tape_file[100];
    
    if (CHAINS) poly_N = N_back * ((N_side * N_chains) + 1);
    else poly_N = N_back;
    N = N_poly * poly_N;  // total number of monomers
    if (BONDS || PARAB) N_bond = N_poly * ((N_back - 1) + (N_chains * N_side * N_back));
    else N_bond = 0;
    printf("Total Number of Monomers = %d\t Total Number of Bonds = %d\n",N,N_bond);
    potshift = 4.0 * LJ_eps * (pow((LJ_sig/LJ_cut),12) - pow((LJ_sig/LJ_cut),6));
    //RandSeed((unsigned)time(NULL));
    RandSeed(13131);

    /* open data files */
    sprintf(tape_file, "Pos.dat");
    if((tapePtr = fopen(tape_file, "w")) == NULL){
        fprintf(stderr, "file %s could not be opened !", tape_file);
        exit(1);
    }
    sprintf(bin_file, "binPos.dat");
    if((binPtr = fopen(bin_file, "w")) == NULL){
        fprintf(stderr, "file %s could not be opened !", bin_file);
        exit(1);
    }
    
    /* Initilization */
    printf("---- Started ----\n");
    dyn_alloc_array();
    if (SURF) surf_initial_pos();
    else initial_pos();
    printf("---- Initial configuration generated ----\n");
    write_pos(2, binPtr); 
    
}

/****************************************************************/
/* calculate and write out the final results for each system  */
void write_conf(FILE *outPtr)
{
    int i,j;
    fprintf(outPtr,"%lld\n",i_step);
    for(i = 0; i < N; i++){
        for(j = 0; j <3; j++)
            fprintf(outPtr, "%15.6le",r[i][j]);
        fprintf(outPtr, "\n");
    }
}
/****************************************************************/
void write_pos(int rec, FILE *outPtr)
{
    int i, k;
    int N3;
    short tb[2];
    double tr[3], time;
    char code;
    
    time = i_step;
    
    if (rec) {
        code = rec;
        fwrite(&code,sizeof(code),1,outPtr);
        if (rec == 1) {
            N3 = N;
            fwrite(&N3,sizeof(N3),1,outPtr);
            fwrite(&time,sizeof(time),1,outPtr);
            for(i = 0; i < N; i++) {
                for (k = 0; k < 3; k++) tr[k] = r[i][k];
                fwrite(&tr[0],sizeof(tr[0]),3,outPtr);
            }
        } else
            if (rec == 2) {
                N3 = N_bond;
                fwrite(&N3,sizeof(N3),1,outPtr);
                fwrite(&time,sizeof(time),1,outPtr);    
                for(i = 0; i < N_bond; i++) {
                    tb[0] = b1[i] + 1;
                    tb[1] = b2[i] + 1; 
                    fwrite(&tb[0],sizeof(tb[0]),2,outPtr);
                }
            }
    }
}

void Calculations(void)
{
    int k;
    double r1,r2,r3,r7,r13,ftens;
    double Ro = dsqr(R_fene);

    r2=0.;
    for (k = 0; k < 3; k++) {
        r2 += dsqr(r[b1[0]][k] - r[b2[0]][k]);
    }

    r1 = sqrt(r2);
    rravg += r2;
    rrmax = max(r2,rrmax);
    rrmin = min(r2,rrmin);
    
    r3 = r2 * r1;
    r7 = r3 * r3 * r1;
    r13 = r7 * r3 * r3;
    ftens = -(K_fene * Ro * r1) / (r2 - Ro);
    ftens -= (24. * (LJ_eps) * (2. * (1. / r13) - (1. / r7)));
    bondTens += ftens;
    
    r2=0.;
    for (k = 0; k < 3; k++) {
        r2 += dsqr(r[0][k] - r[N_back-1][k]);
    }
    chainR2 += r2;
}
/*****************************************************************/

int main(int argc, char **argv)
{
    //int n_tape = 100;
    int p_tape = 10000000;
    //int samp_time = 5;
    eq_step = 1e+05;
    tot_step = 1e+09;
    
    time (&wallClockTime0);
    Setup();
    
    /* Equilibration */
    printf("----- Equilibration Step -----\n\t");
    for(i_step=0; i_step < eq_step; i_step++) {
        if((i_step%t_neighs==0)) NeighborList();
        mcmove();
        if((i_step%p_tape==0) && (i_step > 0)) {
            printf("%lld  ",i_step);
        }
    }
    
    acptcnt = cnt = 0.;
    rrmax = rravg = chainR2 = bondMax = bondMin = bondTens = 0.;
    rrmin = 10.;
    /* main loop */
    printf("\n----- Step Number -----\n\t");
    for(i_step=0; i_step < tot_step; i_step++) {
        if((i_step%t_neighs==0)) NeighborList();
        mcmove();
        //if(i_step%n_tape==0) {
            //write_conf(tapePtr);
            //write_pos(1, binPtr);
        //}
        //if(i_step%samp_time==0) Calculations();
        Calculations();
        if((i_step%p_tape==0) && (i_step > 0)) {
            bondMax += rrmax;
            bondMin += rrmin;
            cnt += 1.;
            printf("%lld  ",i_step);
            fflush(stdout);
            rrmax = 0.;
            rrmin = 10.;
        }
        if((i_step%(tot_step/10)==0) && (i_step > 0)) {
            printf("\n\t");
            fflush(stdout);
        }    
    }
    printf("\n");
    printf("----------- Average Squared End-to-End Distance = %f----------\n", chainR2 / tot_step);
    printf("----------- Average Bond Distance = %f----------\n", rravg / tot_step);
    printf("----------- Average Tension = %f ----------\n", bondTens / tot_step);
    printf("----------- Average Bond Max = %f----------\n", bondMax / cnt);
    printf("----------- Average Bond Min = %f----------\n", bondMin / cnt);
    printf("----------- Acceptance Ratio = %f ----------\n", acptcnt / tot_step);
    printf ("----------- Run time %ld seconds -----------\n", time (NULL) - wallClockTime0);
    printf("----------- Finished Normally -----------\n");
    
    return 0;
}
