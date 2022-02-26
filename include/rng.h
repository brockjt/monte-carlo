#ifndef RNG_H_
#define RNG_H_
/*
   This macro determines which random number generator will be used:
   RANMAR or MT19937 or internal (~1.5 times quicker).
   When none of USE_RANMAR or MT19937 is defined, under 32bit WATCOM C/C++ will be used
   LC RNG from Delphi RTL, under all other systems and compilers RNG from the
   standard library will be used.
*/
//#define USE_RANMAR
#define USE_MT19937

#if defined(__WATCOMC__) || defined(_MSC_VER)
#define INLINE __inline
#elif defined(__GNUC__)
#define INLINE static __inline
#else
#define INLINE static
#endif


#ifdef USE_RANMAR
extern const long __ranmax;

/* RANMAR initializer */
void rmarin(unsigned ij, unsigned kl);
long ranmar(void);
extern void RandSeed(unsigned long seed);

#if defined(__WATCOMC__) && defined(__386__)
unsigned long _ranmar_range(unsigned long range, long value);
#pragma aux _ranmar_range = \
        ".386"              \
        "mul  edx"          \
        "shld edx,eax,8"    \
        parm [edx] [eax]    \
        modify [eax]        \
        value [edx];

#elif defined(__GNUC__) && defined(i386)
INLINE unsigned long _ranmar_range(unsigned long range, long value)
{
    unsigned long __res, ignore;    
    __asm__ (	"mull %%edx\n"
		"shld $8, %%eax, %%edx\n"
	    : "=d" (__res), "=a" (ignore)
	    : "d" (range), "a" (value) );
    return __res;
}

#elif defined(not__GNUC__) || defined(not__INTEL_COMPILER) 
INLINE unsigned long _ranmar_range(unsigned long range, long value)
{ /* disabled, sometimes fails */
  long long i = ((long long)range * (long long)value) >> 24;
  return i;
}
#else 
INLINE unsigned long _ranmar_range(unsigned long range, long value)
{
  unsigned long l1, l2, h1, h2;
  l1 = (unsigned short)range;
  l2 = (unsigned short)value;
  h1 = (unsigned short)(range >> 16);
  h2 = (unsigned short)(value >> 16);
  return ((((l1 * l2) >> 16) + ((h1 * l2) + (h2 * l1))) >> 8) +
         ((h1 * h2) << 8);
/*
  unsigned short l1, l2, h1, h2;
  l1 = (unsigned short)range;
  l2 = (unsigned short)value;
  h1 = (unsigned short)(range >> 16);
  h2 = (unsigned short)((unsigned long)value >> 16);
  return (((unsigned long)l1 * (unsigned long)l2) >> 24) +
         ((((unsigned long)h1 * (unsigned long)l2) + ((unsigned long)h2 * (unsigned long)l1)) >> 8) +
         (((unsigned long)h1 * (unsigned long)h2) << 8);
*/
}
#endif /* WATCOM, GNUC, etc */

#define _MINUS24 (-24)

#if defined(__WATCOMC__) && defined(__386__)
extern double _minus24;
double _ranmar_double(long value);
#pragma aux _ranmar_double = \
        ".386"                            \
        "FLD     _minus24"                \
        "PUSH    0"                       \
        "PUSH    EAX"                     \
        "FILD    qword ptr [ESP]"         \
        "ADD     ESP,8"                   \
        "FSCALE"                          \
        "FSTP    ST(1)"                   \
         parm [eax]                       \
         modify [eax 8087]                \
         value [8087];
#else

INLINE double _ranmar_double(long value)
{
  return (double)value / __ranmax;
}

#endif /* ranmar_double, WATCOM, GNUC, etc */

#define RAND_LIMIT (__ranmax - 1)

#elif defined(USE_MT19937)

typedef unsigned long _mt_uint32;

extern _mt_uint32   *_mt_next;
extern int      _mt_left;

extern _mt_uint32 reloadMT(void);
extern void RandSeed(unsigned long seed);

INLINE _mt_uint32 randomMT(void)
 {
    _mt_uint32 y;

    if(--_mt_left < 0)
        return(reloadMT());

    y  = *_mt_next++;
    y ^= (y >> 11);
    y ^= (y <<  7) & 0x9D2C5680U;
    y ^= (y << 15) & 0xEFC60000U;
    return(y ^ (y >> 18));
 }


#if defined(__GNUC__) && defined(i386)
INLINE unsigned long _mt_range(unsigned long range,  _mt_uint32 value)
{
    unsigned long __res, ignore;    
    __asm__ (	"mull %%edx\n"
	    : "=d" (__res), "=a" (ignore)
	    : "d" (range), "a" (value) );
    return __res;
}

#else 
INLINE unsigned long _mt_range(unsigned long range, _mt_uint32 value)
{
  unsigned long l1, l2, h1, h2;
  l1 = (unsigned short)range;
  l2 = (unsigned short)value;
  h1 = (unsigned short)(range >> 16);
  h2 = (unsigned short)(value >> 16);
  return ((((l1 * l2) >> 16) + ((h1 * l2) + (h2 * l1))) >> 16) + (h1 * h2);
}
#endif /* mt_range, GNUC, etc */

#define RAND_LIMIT (0xFFFFFFFFU)

static const double _mt_inv = 2.3283064365387e-10;
INLINE double _mt_double(_mt_uint32 value)
{
//  static const double __inv = 0.5 / 0x80000000U;
  return _mt_inv * value;
}

#else /* use fast LG or builtin algorithm */

#if defined(__WATCOMC__) && defined(__386__)

 /* This random number generator is taken form DELPHI run-time library */
#define _DELPHI_RANDOM_GEN_

extern unsigned long _randseed;
extern double _minus32;

unsigned long _randint(unsigned long range);
#pragma aux _randint = \
        ".386"                            \
        "imul    edx,_randseed,08088405h" \
        "inc     edx"                     \
        "mov     _randseed,edx"           \
        "mul     edx"                     \
        "mov     eax,edx"                 \
        parm [eax]                        \
        modify [edx]                      \
        value [eax];

double _randdouble(void);
#pragma aux _randdouble = \
        ".386"                            \
        "IMUL    EDX,_randseed,08088405H" \
        "INC     EDX"                     \
        "MOV     _randseed,EDX"           \
        "FLD     _minus32"                \
        "PUSH    0"                       \
        "PUSH    EDX"                     \
        "FILD    qword ptr [ESP]"         \
        "ADD     ESP,8"                   \
        "FSCALE"                          \
        "FSTP    ST(1)"                   \
         modify [edx 8087]                \
         value [8087];

INLINE void RandSeed(unsigned long seed)
{
  _randseed = seed;
}

#define RAND_LIMIT ULONG_MAX

#else

INLINE void RandSeed(unsigned long seed)
{
  srand((unsigned)seed);
}

#define RAND_LIMIT RAND_MAX

#endif

#endif /* USE_RANMAR or USE_MT19937 */

/* Wrapper functions for random number generator */

INLINE unsigned Random(unsigned Range)
{
#if defined (USE_RANMAR)
  return _ranmar_range(Range, ranmar());
#elif defined(USE_MT19937)
  return _mt_range(Range, randomMT());
#else
#ifdef _DELPHI_RANDOM_GEN_
  return _randint(Range);
#else
//  return (unsigned)(((unsigned long)rand()*(unsigned long)Range)/((unsigned)RAND_MAX+1));
  return (unsigned)(((unsigned long long)rand()*(unsigned long long)Range)/((unsigned long long )RAND_MAX+1));
#endif
#endif
}

double INLINE Randomr(void)
{
#ifdef USE_RANMAR
  return _ranmar_double(ranmar());
#elif defined (USE_MT19937)
  return _mt_double(randomMT());
#elif defined(_DELPHI_RANDOM_GEN_)
  return _randdouble();
#else
  return ((double)rand())/((unsigned long)RAND_MAX+1);
#endif
}

extern void Randomize(void);
extern double NormRand(void);
#undef INLINE
#endif
