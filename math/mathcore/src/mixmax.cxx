/*
 *  mixmax.c
 *  A Pseudo-Random Number Generator
 *
 *  Created by Kosta on Sun Feb 22 2004.
 *  As of version 0.99, the code is being released under GNU Lesser General Public License v3
 *
 *	Generator described in 
 *	N.Z.Akopov,G.K.Savvidy and N.G.Ter-Arutyunian, Matrix Generator of Pseudorandom Numbers, 
 *	J.Comput.Phys. 97, 573 (1991); 
 *	Preprint EPI-867(18)-86, Yerevan Jun.1986;
 */


#include <stdio.h>
#include <stdlib.h>

#include "mixmax.h"  // stdint.h is included there

void iterate(rng_state_t *X){
	int i;
	myuint temp2, temp3, temp4;
	
	//myuint temp1 = X->V[1];
	temp2 = X->V[2]; temp3 = X->V[3]; temp4 = X->V[4];
	X->counter = 0;
	X->V[1] = MOD_MERSENNE(X->V[1] + X->partial[N]);
	X->V[2] = MOD_MERSENNE(X->V[1] + X->partial[2]);  // order matters here, first calculate new X, then update the partial sum, next X uses its own partial
	X->V[3] = MOD_MERSENNE(X->V[2] + X->partial[3]);  
	X->V[4] = MOD_MERSENNE(X->V[3] + X->partial[4]);  
	X->V[5] = MOD_MERSENNE(X->V[4] + X->partial[5]);  
	X->V[6] = MOD_MERSENNE(X->V[5] + X->partial[6]);  
#if (N==60 && BITS==61)                                  /* the special entries of the matrix are taken care of here! */
	X->V[3] = MOD_DBL(X->V[3] + 4*temp2 );              	/* The maximal period matrix for N=60 has 737 pattern instead of 3's in these places, 4=7-3 */
	X->V[5] = MOD_DBL(X->V[5] + 4*temp4 );
#else
#ifdef SPECIAL
	X->V[3] = MOD_MERSENNE(X->V[3] + MOD_MULSPEC(temp2) );        /* all the newly recommended N have just one special entry, here we add SPECIAL*V[2] */
#else
	X->V[3] = MOD_MERSENNE(X->V[3] + MERSBASE - temp2);          /* the old matrix had a 2 instead of 3 there, substract V[2] */
#endif
#endif
	X->partial[1] = ZERO;	
	X->partial[2] = MOD_MERSENNE(X->partial[1] + X->V[2]);  
	X->partial[3] = MOD_MERSENNE(X->partial[2] + X->V[3]);  
	X->partial[4] = MOD_MERSENNE(X->partial[3] + X->V[4]);  
	X->partial[5] = MOD_MERSENNE(X->partial[4] + X->V[5]);  
	X->partial[6] = MOD_MERSENNE(X->partial[5] + X->V[6]);  
	for (i=7; i<=N; i++){
		X->V[i] = MOD_MERSENNE(X->V[i-1] + X->partial[i]);   
		X->partial[i] = MOD_MERSENNE(X->partial[i-1] + X->V[i]);   // precalculate for next iteration
	}
}

#define WPRECALC  
// 10% faster overall, and no waiting for a full iteration in between spitting out 

#ifdef WPRECALC
myuint get_next(rng_state_t *X){
	// (void *vstate){  
	// rng_state_t *X = (rng_state_t *) vstate;
	
	myuint temp;
	int i;

	i=X->counter;
	if (i>=N) {i=1;}else{i++;};
	temp = X->V[i];                                  // hold the old X
	if (i <= 6){
		if (i==1){
			X->V[1] = MOD_MERSENNE(X->V[1] + X->partial[N]);
			X->temp = X->V[1];
			X->partial[1] = ZERO;	
		} else {
			// precalculate for next iteration	
			X->V[i] = MOD_MERSENNE(X->temp + X->partial[i]);  // use uncorected X[i-1], is held in X->temp , very important
			X->temp = X->V[i];                               // hold it 
#if (N==60 && BITS==61)                              /* the special entries of the matrix are taken care of here! */
			if (i==3 || i==5) { X->V[i] = MOD_DBL(X->V[i] + 4*X->PrevValue );}        /* The maximal period matrix for N=60 has 7 instead of 3 in these places, 4=7-3 */
#elif ( (N==3150 || N==1260 || N==508 || N==256 || N==240 ) && BITS==61)
			if (i==3) { X->V[3] = MOD_MERSENNE(X->V[3] + MOD_MULSPEC(X->PrevValue)  );} 
			// the charpoly is irreducible for these N and has maximal period for N=508, 256, half period for 1260, and 1/12 period for 3150
#else
			if (i==3) { X->V[3] = MOD_MERSENNE(X->V[3] + MERSBASE - X->PrevValue );}     /* the old matrix has a 2 instead of 3 there, substract old X[2] */
#endif	
			X->PrevValue = temp;  // store it for later
			X->partial[i] = MOD_MERSENNE(X->partial[i-1] + X->V[i]); // use actual new value for partial sum     	
		}
	}else{
		X->V[i] = MOD_MERSENNE(X->temp + X->partial[i]);
		X->temp = X->V[i];
		X->partial[i] = MOD_MERSENNE(X->partial[i-1] + X->V[i]);   // precalculate for next iteration
	}
	X->counter = i;
	return temp ;
}
#else
myuint get_next(rng_state_t *X){
	// (void *vstate){  
	// rng_state_t *X = (rng_state_t *) vstate;
	
	int i;
	i=X->counter;
	if (i>=N){
		i=1;
		IT(X);
	} else {
		i++;
	}
	X->counter = i;       // beware, i points to the X which is already output, below
	return X->V[i];
}
#endif

double get_next_float(rng_state_t *X){
	return ( ( (double)get_next(X) ) * (double)(INV_MERSBASE));
}

rng_state_t* rng_alloc() 
{
/* allocate the state */
	rng_state_t  *p = (rng_state_t*)malloc(sizeof(rng_state_t)); // 1 for counter, N for the vector 1..N, 2 more for temps, and N for partial sums 1..N
	p->fh=stdout; // by defualt, set the output file handle to stdout  
	return p;
}

int rng_free(rng_state_t *X) /* free memory occupied by the state */
{
	free(X);
	return 0;
}

rng_state_t*  rng_copy(myuint *Y)
{
	int i;
	/* copy the state Y, returns pointer to the newly allocated and initialized state,  
	 it is users responsibility  to make sure Y is properly allocated with rng_alloc, then pass Y->V or is an array -- myuint Y[N+1] and Y[1]...Y[N] are set to legal values */
	/* partial sums on this new state are recalculated, and counter set to zero, if get_next is used will output all values in it until... */
	rng_state_t *X = rng_alloc();
	X->counter = 0;
	X->partial[1] = ZERO;
	X->V[1] = Y[1];
	for (i=2; i <= N; i++){
		X->V[i] = Y[i]; 
		X->partial[i] = MOD_MERSENNE(X->partial[i-1] + X->V[(i)]); 
	}
	if (X->fh==NULL){X->fh=stdout;}
	return X;
}

void seedvielbein(rng_state_t *X, unsigned int seed)	
{
int i;
	if (seed<=N && seed>=1){
		for (i=0; i <= N; i++){
			X->V[i] = ZERO; 
		}
		X->V[seed] = SEEDVAL;
	}else{
		fprintf(stderr, "Out of bounds seed index, is not ( <=n and >=1 )\n"); exit(-1);
	}
	X->temp=ZERO; 
	X->PrevValue=ZERO;
	X->counter = N;  // set the counter to N if iteration should happen right away
	precalc(X);
	if (X->fh==NULL){X->fh=stdout;}	
}

int seed_lcg(rng_state_t *X, myuint seed)	/* an LCG(MULT, 2^61-1) generator of Wu is used to seed */
{
#define MULT 1073217536
	int i;
	myuint temp;
	if (seed == 0){
		seed = SEEDVAL;
	}

	temp =  (myuint)seed;
	X->V[1] = temp & MERSBASE;
	if (X->fh==NULL){X->fh=stdout;} // if the filehandle is not set, make it stdout
	fprintf(X->fh, "Seed vector was: %llu",X->V[1]);
	for (i=2; i <= N; i++){
		X->V[i] = modmulM61(X->V[i-1], MULT);
		fprintf(X->fh, ",  %llu",X->V[i]);
	}
	fprintf(X->fh,"\n");
	precalc(X);
	X->temp=ZERO; 
	X->PrevValue=ZERO;
	X->counter = N;  // set the counter to N so that iteration happens right away
	return 0;
}

int precalc(rng_state_t *X){
	int i;
	X->partial[0] = ZERO;
	X->partial[1] = ZERO;
	for (i=2; i <= N; i++){
		X->partial[i] = MOD_MERSENNE(X->partial[i-1] + X->V[(i)]);  //partial sums, excluding X[1]
	}
	return 0;
}


#if __LP64__
#if ((N==60) && BITS==61)
rng_state_t *skip_1M(rng_state_t *X){
	/*
	 * Skips 2^20 = 1048576 steps, by applying a precomputed matrix to X
	 */
#if (N==60)
	__uint128_t skipMat[N][N] = 
#include "skip1M_N60.c"
	;
#endif
	int i,j;
	__uint128_t temp, Y[N+1];
	//__uint64_t *ptr, t1,t2;
	Y[0]=0;
	for (i=0; i<N; i++){
		Y[i+1] = 0;
		for (j=0; j<N; j++){
			temp = Y[i+1] + ( skipMat[j][i] * (__uint128_t)X->V[j+1]) ;
			Y[i+1] = temp - (temp>>BITS)*MERSBASE;
		}
	}
	for (j=0; j<=N; j++){ X->V[j] = (myuint)Y[j];}
	precalc(X);
	return X;
}

rng_state_t *skip_1G(rng_state_t *X){
	/*
	 * Skips 2^30 = 1073741824 steps, by applying a precomputed matrix to X
	 */
#if (N==60) 
__uint128_t skipMat[N][N] = 
#include "skip1G_N60.c"
	;
#endif
	int i,j;
	__uint128_t temp, Y[N+1];
	//__uint64_t *ptr, t1,t2;
	Y[0]=0;
	for (i=0; i<N; i++){
		Y[i+1] = 0;
		for (j=0; j<N; j++){
			/*
			 temp = ( skipMat[j][i] * (__uint128_t)X->V[j+1]);
			temp = Y[i+1] + temp;
			ptr = &temp; // the only place where endianness matters
			t1 = *ptr;
			ptr++;
			t2 = *ptr;
			//			fprintf(stderr, "Lower and upper parts are %llu , %llu\n", t1, t2);
			//			fprintf(stderr, "Will sub 2^61-1 x n, n=%llu", 8*t2 + (t1>>61));
			temp = temp - ((1<<(64-BITS))*t2 + (t1>>BITS))*MERSBASE;
			Y[i+1] = (temp); 
			 */
			temp = Y[i+1] + ( skipMat[j][i] * (__uint128_t)X->V[j+1]) ;
			Y[i+1] = temp - (temp>>BITS)*MERSBASE;
		}
	}
	for (j=0; j<=N; j++){ X->V[j] = (myuint)Y[j];}
	precalc(X);
	return X;
}
#endif
#endif

int rng_get_N(void){return N;}

void printbitssimple(uint32_t n) {
	uint8_t i = 32 ; // (sizeof(n) * 8 - 1);
	register uint32_t tmp=n;
	do {
		i--;
		putchar(0x30 + (tmp&1) ); /* 0 is ASCII 0x30 and 1 is ASCII 0x31  */
		tmp=tmp>>1;
	} while (i > 0);
	putchar(0x0A);//putchar(0x0C);// newline 0A, carriage return 0D, form feed 0C
}

inline uint64_t modmulM61(uint64_t s, uint64_t a)
{
#define MASK32 0xFFFFFFFFULL
	register uint64_t o,ph,pl,ah,al,tmp;
	o=(s)*a;
	ph = ((s)>>32);
	pl = (s) & MASK32;
	ah = a>>32;
	al = a & MASK32;
	o = (o & M61) + ((ph*ah)<<3) + ((ah*pl+al*ph + ((al*pl)>>32))>>29) ;
	if (o >= M61){ tmp = o-M61; o = tmp; } // "c" putchar(0x63); 	
    return o;
}

