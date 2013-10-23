#include <stdint.h>

#ifndef MIXMAX_H_
#define MIXMAX_H_

#ifdef __cplusplus
extern "C" {
#endif
	
#ifndef _N
#define N 256 
/* The currently recommended N are 3150, 1260, 508, 256, 240, 88
   Since the algorithm is linear in N, the cost per number is almost independent of N.
 */
#else
#define N _N
#endif

#ifndef __LP64__
typedef uint64_t myuint;
//#warning typedefining 'myuint' as 'uint64_t'
#else
typedef unsigned long long int myuint;
//#warning typedefining 'myuint' as 'unsigned long long int'
#endif

typedef struct
{
    myuint V[N];
	uint64_t sumtot;
	myuint tempP;
	myuint tempV;
    int counter;
	FILE* fh;
}
rng_state_t;


rng_state_t  *rng_alloc();                 /* allocate the state */
int           rng_free(rng_state_t* X);    /* free memory occupied by the state */
rng_state_t  *rng_copy(myuint *Y);         /* init from vector, takes the vector Y, returns pointer to the newly allocated and initialized state */
int           rng_get_N(void);


//   FUNCTIONS FOR SEEDING

typedef uint32_t myID_t;

void seed_uniquestream(rng_state_t* X, myID_t clusterID, myID_t machineID, myID_t runID, myID_t  streamID );
/*
 best choice: will make a state vector from which you can get at least 10^100 numbers 
 guaranteed mathematically to be non-colliding with any other stream prepared from another set of IDs, 
 so long as it is different by at least one bit in at least one of the four IDs
			-- useful if you are running a parallel simulation with many clusters, many CPUs each
 */

void seed_spbox(rng_state_t* X, uint64_t seed);    // non-linear method, makes certified unique vectors, but only probabilistically non-colliding

void seed_lcg(rng_state_t* X, myuint seed);        // seeds with a LCG(MULT, 2^61-1), only ok for non-parallel simulations

void seed_vielbein(rng_state_t* X, unsigned int i); // seeds with the i-th unit vector, i = 0..N-1

void seed_CTR_1x64(rng_state_t* X, uint64_t a);     // experimental counter mode, see driver_CTRmode.c for example usage
void seed_CTR_1x128(rng_state_t* X, __uint128_t a);

	int iterate(rng_state_t* X);
	myuint iterate_raw_vec(myuint* Y, myuint sumtotOld);


//   FUNCTIONS FOR GETTING RANDOM NUMBERS


#ifdef __MIXMAX_C
	myuint get_next(rng_state_t* X);         // returns 64-bit int, which is between 1 and 2^61-1 inclusive
	double get_next_float(rng_state_t* X);   // returns double precision floating point number in (0,1]

#else
#define get_next(X) GET_BY_MACRO(X)
	myuint GET_BY_MACRO(rng_state_t* X) {
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
		//return ( ( (X->counter)<(N) ) || iterate(X) ) + X->V[X->counter++ ] 	;
	}	
	
#define get_next_float(X) get_next_float_BY_MACRO(X)
	
	double get_next_float_BY_MACRO(rng_state_t* X){
		return ( ( GET_BY_MACRO(X) ) * ((0x1p-61)));
	}
	
#endif  //__MIXMAX_C
	

void fill_array(rng_state_t* X, unsigned int n, double *array);  // fastest method: set n to a multiple of N (e.g. n=256)

void iterate_and_fill_array(rng_state_t* X, int stride, double *array); // fills the array with N numbers starting at array+N*stride

void print_state(rng_state_t* X);


int precalc(rng_state_t* X); 
/* needed if the state has been changed by something other than  iterate, but no worries, seeding functions call this for you if necessary */
void iterate_backward(rng_state_t* X); // iterates one step backwards
myuint apply_bigskip(myuint* Vout, myuint* Vin, myID_t clusterID, myID_t machineID, myID_t runID, myID_t  streamID );
// applies a skip of some number of steps calculated from the four IDs
void branch_inplace( rng_state_t* Xin, myID_t* ID ); // almost the same as apply_bigskip, but in-place and from a vector of IDs

uint64_t mod128(__uint128_t s);
uint64_t modmulM61(uint64_t s, uint64_t a);
uint64_t fmodmulM61(uint64_t cum, uint64_t s, uint64_t a);

void doRounds(rng_state_t* X);

//void printbitssimple(uint32_t n) ;

#define IT(a) iterate(a)

#ifndef _BITS
#define BITS  61
#endif

/* Mersenne Numbers */
#define M61   2305843009213693951

#define SUFF(i) i##ULL
#define xSUFF(i) SUFF(i)
#define ZERO SUFF(0)
#define ONE SUFF(1)
#define SEEDVAL SUFF(1)

#define ROCL( value, shift) ( ((value) << (shift) ) | ((value) >> (sizeof(value) * 8 - (shift) )) )
#define ROCR( value, shift) ( ((value) >> (shift) ) | ((value) << (sizeof(value) * 8 - (shift) )) )

#define MERSBASE xSUFF(M61) // 2^61-1 = 2305843009213693951
#define MOD_KOSTAS(k) ((k) - (((k)+1) >> BITS)*MERSBASE ) 
#define MOD_PAYNE(k) ((((k)) & MERSBASE) + (((k)) >> BITS) )  // slightly faster than mine, ok for addition
#define MOD_ROCL(k) ((((k)) & MERSBASE) + (ROCL( (k), (64 - BITS))&7) ) 
#define MOD_IF(k)    ( ((k)<MERSBASE ) ? (k): ((k)-MERSBASE) )
#define MOD_REM(k) ((k) % MERSBASE )  // latest Intel CPU is supposed to do this in one CPU cycle, but on my machines it seems to be 20% slower than my trick
#define MOD_MUL(m,k) ((m)*(k) - ( ((m)*(k)) >> BITS)*MERSBASE ) // m must be less or equal to 8
#define NEG(k) ((k) ? (MERSBASE-(k)): ZERO)
#define MOD_DBL(k) MOD_REM( (k) ) // must work in all cases
//#define MOD_MERSENNE(k) MOD_ROCL(k)
#define MOD_MERSENNE(k) MOD_PAYNE(k)

//#define INV_MERSBASE 0.433680868994201773791060216479542685926876E-18L // that is 1/(2^61-1), this is for numbers in [0,1), if numbers in [0,1] are desired, use 1/(2^61-2)
#define INV_MERSBASE (0x1p-61)



#if (BITS==61)
// the charpoly is irreducible for these combinations of N and SPECIAL and has maximal period for N=508, 256, half period for 1260, and 1/12 period for 3150

#if (N==1260)
#define SPECIAL 15
#define MOD_MULSPEC(k)  ( (3*( ( 5*(k) ) % MERSBASE)) % MERSBASE) 

#elif (N==3150)
#define SPECIAL -11
#define MOD_MULSPEC(k) (( ( ( 6*(MERSBASE-(k) )) % MERSBASE) + ( ( 5*(MERSBASE-(k) )) % MERSBASE) ) % MERSBASE )

#elif (N==1000)
#define SPECIAL 0
#define MOD_MULSPEC(k) (0)

#elif (N==720)
#define SPECIAL 1
#define MOD_MULSPEC(k) (k)

#elif (N==508)
#define SPECIAL 5
#define MOD_MULSPEC(k) (( SPECIAL*(k)) % MERSBASE)

#elif (N==256)
#define SPECIAL -1
#define MOD_MULSPEC(k) (MERSBASE-(k))

#elif (N==88)
#define SPECIAL 1
#define MOD_MULSPEC(k) (k)

#elif (N==64)
#define SPECIAL 6
#define MOD_MULSPEC(k) (( SPECIAL*(k)) % MERSBASE)

#elif (N==44)
#define SPECIAL 0
#define MOD_MULSPEC(k) 0

#elif (N==40)
#define SPECIAL 1
#define MOD_MULSPEC(k) (k)

#elif (N==30)
#define SPECIAL 3
#define MOD_MULSPEC(k) (( SPECIAL*(k)) % MERSBASE)

#elif (N==16)
#define SPECIAL 6 
#define MOD_MULSPEC(k) (( SPECIAL*(k)) % MERSBASE)

#elif (N==10)
#define SPECIAL -1 
#define MOD_MULSPEC(k)  (MERSBASE-(k) )

#elif (N==4)
#define SPECIAL 1073217533  // not like the others, uses a large special, but is still a mixmax with determinant equal to 1 
#define MOD_MULSPEC(k) ( modmulM61(SPECIAL, (k)) )
//#define MOD_MULSPEC(k) (( SPECIAL*(k)) % MERSBASE)

#else
#define SPECIAL -1
#define MOD_MULSPEC(k)  (MERSBASE-(k) )

#endif // list of interesting N for 61-BITS ends here

#else
#define SPECIAL -1
#define MOD_MULSPEC(k)  (MERSBASE-(k) )

#endif


//#define FLOAT_OUTPUT 1
#ifndef FLOAT_OUTPUT
#define OUT(x) printf("%llu,\t", x);
#else
#define OUT(x) printf("%1.20LF\t", (long double)x* INV_MERSBASE); // 20 digits of which 18 or 19 are good
#endif

#ifdef __cplusplus
}
#endif
		
#endif // closing MIXMAX_H_