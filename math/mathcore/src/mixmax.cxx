/*
 *  mixmax.c
 *  A Pseudo-Random Number Generator
 *
 *  Created by Konstantin Savvidy on Sun Feb 22 2004.
 *  As of version 0.99 and later, the code is being released under GNU Lesser General Public License v3
 *
 *	Generator described in 
 *	N.Z.Akopov, G.K.Savvidy and N.G.Ter-Arutyunian, Matrix Generator of Pseudorandom Numbers, 
 *	J.Comput.Phys. 97, 573 (1991); 
 *	Preprint EPI-867(18)-86, Yerevan Jun.1986;
 */


#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#define __MIXMAX_C  // do NOT define it in your own program, just include mixmax.h

#include "mixmax.h"

int iterate(rng_state_t* X){
	X->counter = 0;
	X->sumtot = iterate_raw_vec(X->V, X->sumtot);
	return 0;
}
 
inline myuint iterate_raw_vec(myuint* Y, uint64_t sumtotOld){
	// operates with a raw vector, uses known sum of elements of Y
	int i;
	myuint  tempP, tempV;
#if (SPECIAL != 0)
	myuint temp2 = Y[1];
#endif
	Y[0] = (tempV = MOD_MERSENNE(Y[0] + sumtotOld));
	__uint128_t sumtot = 0; // will keep a running sum of all new elements (except Y[0])
	tempP = 0;              // will keep a partial sum of all old elements (except Y[0])
	for (i=1; i<N; i++){
		tempP = MOD_MERSENNE(tempP + Y[i]);
		Y[i] = (tempV = MOD_MERSENNE(tempV + tempP) );   
		sumtot += tempV;
	}
#if (SPECIAL != 0)
	temp2 = MOD_MULSPEC(temp2);
	Y[2] = MOD_MERSENNE( Y[2] + temp2 );
	sumtot += temp2;
#endif
	return mod128(sumtot);
}


void iterate_and_fill_array(rng_state_t* X, int stride, double *array){
	myuint* Y=X->V;
	int i;
	myuint  tempP, tempV;
#if (SPECIAL != 0)
	myuint temp2 = Y[1];
#endif
	Y[0] = (tempV = MOD_MERSENNE(Y[0] + X->sumtot));
	array[0+stride*N] = (double)tempV * (double)(INV_MERSBASE);
	__uint128_t sumtot = 0; // will keep a running sum of all new elements (except Y[0])
	tempP = 0;             // will keep a partial sum of all old elements (except Y[0])
	for (i=1; i<N; i++){
		tempP = MOD_MERSENNE(tempP + Y[i]);
		Y[i] = (tempV = MOD_MERSENNE(tempV + tempP) );   
		sumtot += tempV;
		array[i+stride*N] = (double)tempV * (double)(INV_MERSBASE);
	}
#if (SPECIAL != 0)
	temp2 = MOD_MULSPEC(temp2);
	Y[2] = MOD_MERSENNE( Y[2] + temp2 );
	array[2+stride*N] = (double)Y[2] * (double)(INV_MERSBASE);
	sumtot += temp2;
#endif
	X->sumtot = mod128(sumtot);
	X->counter = 0;
}




myuint get_next(rng_state_t* X){	
	int i;
	i=X->counter;
	if (i>=(N) ){
		i=0;
		iterate(X);
		X->counter=1;
		return X->V[0];
	}
	else{
	X->counter++;
  }
	return X->V[i];
}



double get_next_float(rng_state_t* X){
	return ( ( get_next(X) ) * (INV_MERSBASE));
}

void fill_array(rng_state_t* X, unsigned int n, double *array)
{
	// Return an array of n random numbers uniformly distributed in 0,1]
	unsigned int i=0,j;
	while ( i<(n/N) ){
		iterate_and_fill_array(X, i, array);
		i++;
	}
	if ((n % N)) {
	iterate(X);
	for (j=0; j< (n % N); j++){
		array[N*i+j] = X->V[j] * (double)(INV_MERSBASE);
	}
	}
//	X->counter = j; // set this if you want to continue with single fetches from the exact spot, but it is not usually necessary 
}


rng_state_t* rng_alloc() 
{
/* allocate the state */
	rng_state_t  *p = (rng_state_t*)malloc(sizeof(rng_state_t)); 
	p->fh=stdout; // by default, set the output file handle to stdout  
	return p;
}

int rng_free(rng_state_t* X) /* free the memory occupied by the state */
{
	free(X);
	return 0;
}

rng_state_t*  rng_copy(myuint *Y)
{
	/* copy the vector stored at Y, and return pointer to the newly allocated and initialized state.  
	 It is the user's responsibility  to make sure that Y is properly allocated with rng_alloc, 
	 then pass Y->V or it can also be an array -- such as myuint Y[N+1] and Y[1]...Y[N] have been set to legal values [0 .. MERSBASE-1]
	 Partial sums on this new state are recalculated, and counter set to zero, so that when get_next is called, 
	 it will output the initial vector before any new numbers are produced, call iterate(X) if you want to advance right away */
	rng_state_t* X = rng_alloc();
	__uint128_t sumtmp;

	X->counter = 0;
	sumtmp = 0;
	X->V[0] = Y[0];
	for (int i=1; i < N; i++){
		X->V[i] = Y[i]; 
		sumtmp +=  X->V[(i)] ; 
	}
	X->sumtot = mod128(sumtmp);
	return X;
}

void seed_vielbein(rng_state_t* X, unsigned int index)	
{
int i;
	if (index<N){
		for (i=0; i < N; i++){
			X->V[i] = ZERO; 
		}
		X->V[index] = SEEDVAL;
	}else{
		fprintf(stderr, "Out of bounds index, is not ( 0 <= index < N  )\n"); exit(-1);
	}
	X->counter = N;  // set the counter to N if iteration should happen right away
	//precalc(X);
	X->sumtot = (index ? 1:0);
	if (X->fh==NULL){X->fh=stdout;}	
}

void seed_lcg(rng_state_t* X, myuint seed)	/* an LCG(MULT, 2^61-1) generator of Wu is used to seed */
{
#define MULT 1073217536
	int i;
	myuint temp;
	__uint128_t sumtmp;
	sumtmp = 0; 
	if (seed == 0){
		fprintf(stderr, " try seeding with nonzero seed next time!");
		seed = SEEDVAL;
	}

	temp =  (myuint)seed;
	X->V[0] = temp & MERSBASE;
	if (X->fh==NULL){X->fh=stdout;} // if the filehandle is not yet set, make it stdout
	fprintf(X->fh, "Seed vector of the RNG was: {%llu",X->V[0]);
	for (i=1; i < N; i++){
		X->V[i] = modmulM61(X->V[i-1], MULT);
		fprintf(X->fh, ", %llu",X->V[i]);
		sumtmp += (X->V[i]);
	}
	fprintf(X->fh," }\n");
	X->counter = 0;  // set the counter to N if we want  iteration to happen right away
	X->sumtot= mod128(sumtmp);
}

void seed_spbox(rng_state_t* X, uint64_t seed)
{ // a 64-bit LCG from Knuth line 26, in combination with a bit swap is used to seed
	const uint64_t MULT64=6364136223846793005ULL; 
	int i;
	__uint128_t sumtmp;
	if (seed == 0){
		fprintf(stderr, " try seeding with nonzero seed next time!");
		exit(-1);
	}
	
	uint64_t l = seed;
	sumtmp = 0; 

	X->V[0] = l & MERSBASE;
	if (X->fh==NULL){X->fh=stdout;} // if the filehandle is not yet set, make it stdout
	for (i=1; i < N; i++){
		l*=MULT64; l = (l << 32) ^ (l>>32);
		X->V[i] = l & MERSBASE;
		sumtmp += (X->V[i]);
	}
	X->counter = 0;  // set the counter to N if iteration should happen right away
	X->sumtot= mod128(sumtmp);
}

int precalc(rng_state_t* X){
	int i;
	myuint temp;
	temp = ZERO; 
	for (i=1; i < N; i++){
		temp = MOD_MERSENNE(temp + X->V[i]);
	}	
	X->sumtot = temp; 
	return 0;
}


int rng_get_N(void){return N;}

#define MASK32 0xFFFFFFFFULL

inline uint64_t mod128(__uint128_t s){
	uint64_t s1;
	s1 = ( (  ((uint64_t)s)&MERSBASE )    + (  ((uint64_t)(s>>64)) * 8 )  + ( ((uint64_t)s) >>BITS) );
	return	MOD_MERSENNE(s1);
}

inline uint64_t modmulM61(uint64_t a, uint64_t b){
	// my best modmul so far
	__uint128_t temp;
	temp = (__uint128_t)a*(__uint128_t)b;
	return mod128(temp);
}


inline uint64_t fmodmulM61(uint64_t cum, uint64_t a, uint64_t b){
	__uint128_t temp;
	temp = (__uint128_t)a*(__uint128_t)b + cum;
	return mod128(temp);
}

void print_state(rng_state_t* X){
    int j;
	for (j=0; (j< (rng_get_N()-1) ); j++) {
		fprintf(X->fh, "%llu, ", X->V[j] );
	}
	fprintf(X->fh, "%llu", X->V[rng_get_N()-1] );
}

#define FUSEDMODMULVEC \
{ for (i =0; i<N; i++){         \
cum[i] = MOD_MERSENNE(cum[i] + modmulM61( coeff ,  Y[i] )) ; \
} }


#define SKIPISON 1

#if ( ( (N==88)||(N==256) ) && BITS==61 && SKIPISON!=0)
void seed_uniquestream( rng_state_t* Xin, myID_t clusterID, myID_t machineID, myID_t runID, myID_t  streamID ){
	FILE* fin; 
	if(  ( fin = fopen("motherseed.conf", "r") ) ){
		char l=0;
		while ( l != 0x7B ) {
			l=fgetc(fin); // proceed until hitting opening bracket
		}
		l=0;
		for(int i = 0; i < rng_get_N(); i++){
			fscanf(fin, "%llu, ", &Xin->V[i]);
		}
		fscanf(fin, "%llu", &Xin->V[rng_get_N()-1]);
		//printf("read in successfully\n\n"); print_state(X);
		if (clusterID || machineID){printf("Please dont try to reapply skipping with non-zero clusterID/machineID to your mother seed in motherseed.conf\n");}
		Xin->sumtot = apply_bigskip(Xin->V, Xin->V,  clusterID,  machineID,  runID,   streamID );
	}else{
		//printf("could not open the file for reading in the mother seed, we apply it ourself\n" );
		seed_vielbein(Xin,0);
		Xin->sumtot = apply_bigskip(Xin->V, Xin->V,  clusterID,  machineID,  runID,   streamID );
	}
}

void branch_inplace( rng_state_t* Xin, myID_t* IDvec ){
	Xin->sumtot = apply_bigskip(Xin->V, Xin->V,  IDvec[3],  IDvec[2],  IDvec[1],   IDvec[0] );
}

uint64_t apply_bigskip(myuint* Vout, myuint* Vin, myID_t clusterID, myID_t machineID, myID_t runID, myID_t  streamID ){
	/*
	 makes a derived state vector, Vout, from the mother state vector Vin
	 by skipping a large number of steps, determined by the given seeding ID's
	 
	 it is mathematically guaranteed that the substreams derived in this way from the SAME (!!!) Vin will not collide provided
	 1) at least one bit of ID is different
	 2) less than 10^100 numbers are drawn from the stream 
	 (this is good enough : a single CPU will not exceed this in the lifetime of the universe, 10^19 sec, 
	 even if it had a clock cycle of Planch time, 10^44 Hz )
	 
	 Caution: never apply this to a derived vector, just choose some mother vector Vin, for example the unit vector by seed_vielbein(X,0), 
	 and use it in all your runs, just change runID to get completely nonoverlapping streams of random numbers on a different day.
	 
	 clusterID and machineID are provided for the benefit of large organizations who wish to ensure that a simulation 
	 which is running in parallel on a large number of  clusters and machines will have non-colliding source of random numbers.
	 
	 did i repeat it enough times? the non-collision guarantee is absolute, not probabilistic
	 
	 */
	
	
	const	uint64_t skipMat[128][N] = 
	
#if (N==88) 
#include "mixmax_skip_N88.c"  // to make this file, delete all except some chosen 128 rows of the coefficients table
#elif (N==256) 
#include "mixmax_skip_N256.c"
//#include "mixmax_skip_N256.dev.c"
#endif
	;
	
	myID_t IDvec[4] = {streamID, runID, machineID, clusterID};
	int r,i,j,  IDindex;
	myID_t id;
	myuint Y[N], cum[N];
	myuint coeff;
	myuint* rowPtr;
	__uint128_t sumtot = 0;
	
	for (i=0; i<N; i++) { Y[i] = Vin[i]; sumtot += Vin[i]; } ; sumtot -= Vin[0]; sumtot = mod128(sumtot) ;
	for (IDindex=0; IDindex<4; IDindex++) { // go from lower order to higher order ID
		id=IDvec[IDindex];
		//printf("now doing ID at level %d, with ID = %d\n", IDindex, id);     
		r = 0;
		while (id){
			if (id & 1) { 
				rowPtr = (uint64_t*)skipMat[r + IDindex*8*sizeof(myID_t)];
				//printf("free coeff for row %d is %llu\n", r, rowPtr[0]);
				for (i=0; i<N; i++){ cum[i] = 0; }    
				for (j=0; j<N; j++){              // j is lag, enumerates terms of the poly
					// for zero lag Y is already given
					coeff = rowPtr[j]; // same coeff for all i
					FUSEDMODMULVEC;
					sumtot = iterate_raw_vec(Y, sumtot); 
				}
				sumtot=0;
				for (i=0; i<N; i++){ Y[i] = cum[i]; sumtot += cum[i]; } sumtot -= Y[0]; sumtot = mod128(sumtot) ;	
			}
		id = (id >> 1); r++; // bring up the r-th bit in the ID		
		}		
	}
	sumtot=0;
	for (i=0; i<N; i++){ Vout[i] = Y[i]; sumtot += Y[i]; } ; sumtot -= Y[0]; // returns sumtot, and copy the vector over to Vout 
	return mod128(sumtot) ;
}


#endif

void seed_CTR_1x64(rng_state_t* X, uint64_t a){
	// X needs to be allocated before calling this, we take care of the rest	
	seed_spbox(X, a);
	doRounds(X);
}

void seed_CTR_1x128(rng_state_t* X, __uint128_t a){
	// X needs to be allocated before calling this, we take care of the rest
	__uint128_t l=a ;
	__uint128_t sumtmp=0;

	const __uint128_t MULT128 = ((__uint128_t)6364136223846793005<<64) + ((__uint128_t)1745191708<<32) + ((__uint128_t)1546275797); 
	// MULT has been checked to be a bi-jection on Z/2^128*Z
	X->V[0]=l&MERSBASE;
	int index=1;
	while (index<N){ // Fill up
		l*=MULT128; l = (l << 64) ^ (l>>64); // we do need something nonlinear so long as it does not create collisions
		X->V[index] = ((uint64_t)l) & MERSBASE;
		sumtmp += (X->V[index]);
		index++;
	}
	X->sumtot = mod128(sumtmp);
	doRounds(X);
	
}


void doRounds(rng_state_t* X) {iterate(X);}

