#include <stdint.h>

#ifndef MIXMAX_H_
#define MIXMAX_H_

#ifndef BITS
#define BITS  61
#endif

/* Mersenne Numbers */
#define M7    127
#define M31   2147483647
#define M61   2305843009213693951
#define M(B)  M##B
#define xM(B) M(B)
#define BASE  xM(BITS)

#define SUFF(i) i##ULL
#define xSUFF(i) SUFF(i)
#define ZERO SUFF(0)
#define ONE SUFF(1)
#define SEEDVAL SUFF(1)
#define xBEXP(i) 0x1p-##i
#define BEXP(i) xBEXP(i) 


#define MERSBASE xSUFF(BASE) // 2^61-1 = 2305843009213693951
#define MOD_KOSTAS(k) ((k) - (((k)+1) >> BITS)*MERSBASE ) 
#define MOD_PAYNE(k) (((k) & MERSBASE) + ((k) >> BITS) )
#define MOD_IF(k) ((k)<MERSBASE ) ? (k): ((k)-MERSBASE) )
#define MOD_REM(k) ((k) % MERSBASE )  // latest Intel CPU supposed to do this in one CPU cycle, but on my machines it seems to be 20% slower than my trick
#define MOD_MUL(m,k) ((m)*(k) - ((m)*((k)) >> BITS)*MERSBASE ) // m must be less or equal to 8
#define MOD_DBL(k) MOD_MERSENNE( MOD_PAYNE(k) ) // do it recursively, twice
#define NEG(k) ((k) ? (MERSBASE-(k)): ZERO)
#define MOD_MERSENNE(k) MOD_KOSTAS(k)

/*inline myuint MOD_MERSENNE(myuint k){
 //if ((k)<MERSBASE){ return (k); }else{	return ((k)-MERSBASE);}
 return ((k) - (((k)+1) >> BITS)*MERSBASE );
 }
 */

//#define INV_MERSBASE 0.433680868994201773791060216479542685926876E-18L // that is 1/(2^61-1), for numbers in [0,1), if numbers in [0,1] are desired, use 1/(2^61-2), but remember, this only matters only once in a quintillion.
//#define INV_MERSBASE (1.0L/(0x1p61L-1.0L))
#define INV_MERSBASE (BEXP(BITS)*(1.0+BEXP(BITS) + BEXP(BITS)*BEXP(BITS)))

#ifndef _N
#define N 256 //508 //3150 //1260 //508 //256 //240 //60
/* The currently recommended N are 3150, 1260, 508, 256, 240, 60
   Will generate (and print) a billion numbers in under 4min, on Intel Xeon 2.8G . 
   Printing, I think is the more expensive part! 
   Since the algorithm is linear in N, the cost per number is roughly independent of N.
   */
#else
#define N _N
#endif

#ifndef SPECIAL

#if (BITS==61)

#if (N==1260)
#define SPECIAL 15
#define MOD_MULSPEC(k)  ( (3*( ( 5*(k) ) % MERSBASE)) % MERSBASE) 

#elif (N==3150)
#define SPECIAL -11
#define MOD_MULSPEC(k) (( ( ( 6*(MERSBASE-(k) )) % MERSBASE) + ( ( 5*(MERSBASE-(k) )) % MERSBASE) ) % MERSBASE )

#elif (N==256)
#define SPECIAL -1
#define MOD_MULSPEC(k) NEG(k)

#elif (N==240)
#define SPECIAL 24
#define MOD_MULSPEC(k) (4*( ( 6*(k) ) % MERSBASE) % MERSBASE)

#elif (N==508)
#define SPECIAL 5
#define MOD_MULSPEC(k) (( SPECIAL*(k)) % MERSBASE)

#elif (N==60)
#define SPECIAL 747
//#undef SPECIAL

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
#define SPECIAL 6 //but beware this is different from the SPECIAL=-1 for which the skip matrix was made!
#define MOD_MULSPEC(k) (( SPECIAL*(k)) % MERSBASE)

#elif (N==10)
#define SPECIAL -1 
#define MOD_MULSPEC(k)  (MERSBASE-(k) )

#elif (N==4)
#define SPECIAL 2 
#define MOD_MULSPEC(k) (( SPECIAL*(k)) % MERSBASE)

#else
#define SPECIAL -1
#define MOD_MULSPEC(k)  (MERSBASE-(k) )

#endif // list of interesting N for 61-BITS ends here

#elif (BITS==31)

#if (N==4)
#define SPECIAL 2
#define MOD_MULSPEC(k) (((k)<<1) % MERSBASE)

#elif (N==10)
#define SPECIAL 3 
#define MOD_MULSPEC(k) (( SPECIAL*(k)) % MERSBASE)

#elif (N==16)
#define SPECIAL 0
#define MOD_MULSPEC(k) 0

#elif (N==24)
#define SPECIAL 6
#define MOD_MULSPEC(k) (( SPECIAL*(k)) % MERSBASE)


#else
#define SPECIAL -1
#define MOD_MULSPEC(k)  (MERSBASE-(k) )

#endif // list of interesting N for 31-BITS ends here

#else
// for BITS=7, there are N=16, N=18 and N=26,N=44,N=106 for all of which good SPECIAL is -1
#define SPECIAL -1
#define MOD_MULSPEC(k)  (MERSBASE-(k) )

#endif
#endif


#if (BITS==61)
#ifndef __LP64__
typedef uint64_t myuint;
#else
typedef unsigned long long int myuint;
#endif
#elif (BITS==31)
typedef unsigned long myuint;
#elif (BITS==7)
typedef unsigned char myuint;
#else
#error "unsupported Mersenne base"
#endif

typedef struct
{
    int counter;
    myuint V[N+1];
    myuint temp;
    myuint PrevValue;
    myuint partial[N+1];
	FILE* fh;
}
rng_state_t;


rng_state_t *rng_alloc();        /* allocate the state */
int rng_free(rng_state_t *X);    /* free memory occupied by the state */
rng_state_t *rng_copy(myuint *Y); /* init from vector, takes the vector Y, returns pointer to the newly allocated and initialized state */

void seedvielbein(rng_state_t *X, unsigned int seed);
int seed_lcg(rng_state_t *X, myuint seed);

void iterate(rng_state_t *X);
int precalc(rng_state_t *X); /* needed if the state has been changed by something other that iterate, and/or we are using get_next with precalculation */
void iterate_backward(rng_state_t *X);
myuint get_next(rng_state_t *X);
double get_next_float(rng_state_t *X);
int rng_get_N(void);

rng_state_t *skip_1M(rng_state_t *X);
rng_state_t *skip_1G(rng_state_t *X);

uint64_t modmulM61(uint64_t s, uint64_t a);
void printbitssimple(uint32_t n) ;

//#define IT(a) iterate_backward(a)
#define IT(a) iterate(a)



//#define FLOAT_OUTPUT 1
#ifndef FLOAT_OUTPUT
#define OUT(x) printf("%llu,\t", x);
#else
#define OUT(x) printf("%1.20LF\t", (long double)x* INV_MERSBASE); // 20 digits of which 18 or 19 are good
#endif


/* // comletely untested, gsl prototype declaration
 #include <gsl/gsl_rng.h>
 static const gsl_rng_type mixmax_type =
 {"mixmax",                       // name 
 MERSBASE,                    //  RAND_MAX =2^61-1
 0,                           // RAND_MIN 
 sizeof (rng_state_t),
 &seed_lcg,
 &get_next,
 &get_next_float};
 */

#endif
