#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <rng.h>

#ifdef _DELPHI_RANDOM_GEN_
dword _randseed = 0;
double _minus32 = -32.0;
#endif

void Randomize(void)
#ifdef DEBUG
{ RandSeed( 6142 ); }
#else
{ RandSeed( time(NULL) ); }
#endif

double NormRand1(void)
/* Based on central limiting theorem */
{
  int i;
  double r, r2, sum = 0;
  for (i = 0; i < 12; i++)
    sum += Randomr();
  r = ( sum - 6.0 ) / 4.0;
  r2 = r * r;
  return (((( 0.029899776 * r2 + 0.008355968 ) * r2 + 0.076542912 )
                          * r2 + 0.252408784 ) * r2 + 3.949846138 ) * r;
}

static double __normRandStor;
static int __hasnormRandStor = 0;
double NormRand(void)
/* Based on Box-Muller transformation */
{
    double x1, x2, w;
 
    if (__hasnormRandStor) {
	__hasnormRandStor = 0;
	return __normRandStor;
    } else {
	do {
    	    x1 = 2.0 * Randomr() - 1.0;
    	    x2 = 2.0 * Randomr() - 1.0;
    	    w = x1 * x1 + x2 * x2;
	} while ( w >= 1.0 );

	w = sqrt( (-2.0 * log( w ) ) / w );
	__normRandStor = x1 * w;
	__hasnormRandStor = 1;
	return x2 * w;
    }
}

long ran0(long idum)
/* "Minimal" random number generator of Park and Miller. Returns a uniform
   random deviate between 0.0 and 1.0. Set or reset idum to any integer
   value to initialize the sequence; idum must not be altered between calls
   for successive deviates in a sequence.
--D.G.-- Source for this RNG was Numerical Recipes {C}. Here it is used
   for simple tasks - like supplying ranmar with initial values.
*/
{
  long k;
  k = idum / 127773;

/* Compute idum = (IA * idum) % IM without overflows by Schrage's method. */
  idum = 16807 * (idum - k * 127773) - 2836 * k;
  if (idum < 0) idum += 2147483647;
  return idum;
}


#ifdef USE_RANMAR
/*===========================================================================
 This is the initialization routine for the random number generator RANMAR()
 NOTE: The seed variables can have values between:    0 <= IJ <= 31328
                                                      0 <= KL <= 30081
 The random number sequences created by these two seeds are of sufficient
 length to complete an entire calculation with. For example, if several
 different groups are working on different parts of the same calculation,
 each group could be assigned its own IJ seed. This would leave each group
 with 30000 choices for the second seed. That is to say, this random
 number generator can create 900 million different subsequences -- with
 each subsequence having a length of approximately 10^30.

 Use IJ = 1802 & KL = 9373 to test the random number generator. The
 subroutine RANMAR should be used to generate 20000 random numbers.
 Then display the next six random numbers generated multiplied by 4096*4096
 If the random number generator is working properly, the random numbers
 should be:
           6533892.0  14220222.0  7275067.0
           6172232.0  8354498.0   10633180.0
==========================================================================*/

const unsigned long __ranparm1 = 31328;
const unsigned long __ranparm2 = 30081;
const unsigned long __ranparm3 = 362436;
const unsigned long __ranparm4 = 7654321;
const unsigned long __ranparm5 = 16777213;
const long __ranmax = 16777216;

long ranmar_params_u[97];
long ranmar_params_c;
int ranmar_params_i97;
int ranmar_params_j97;
double _minus24 = _MINUS24;


/* RANMAR - Random number array generator. Translated from FORTRAN source */
/* This is integer implementation - returns numbers in the range 0 - 2^24 */

long ranmar(void)
{
  long uni, l;
  long i = ranmar_params_i97;
  long j = ranmar_params_j97;

  uni = (ranmar_params_u[i] - ranmar_params_u[j]) & (__ranmax - 1);
  ranmar_params_u[i] = uni;
  if (--i < 0) i = 96;
  if (--j < 0) j = 96;
  ranmar_params_i97 = i;
  ranmar_params_j97 = j;
  l = ranmar_params_c - __ranparm4;
  if (l < 0) l += __ranparm5;
  ranmar_params_c = l;
  uni -= l;
  uni &= __ranmax - 1;
  return uni;
}

void RandSeed(unsigned long seed)
{
  unsigned ij, kl;
  seed = ran0(seed);
  ij = seed % __ranparm1;
  kl = ran0(seed) % __ranparm2;
  rmarin(ij, kl);
}

void rmarin(unsigned ij, unsigned kl)
{
  unsigned i, j, k, l, m, ii, jj;
  unsigned long s, t;

  if (ij > __ranparm1) ij = __ranparm1;
  if (kl > __ranparm2) kl = __ranparm2;

  i = ((ij / 177) % 177) + 2;
  j = (ij % 177) + 2;
  k = ((kl / 169) % 178) + 1;
  l = (kl % 169);

  for (ii = 0; ii < 97; ii++)
    {
      s = 0;
      t = __ranmax / 2;
      for (jj = 0; jj < 24; jj++)
        {
          m = (((i * j) % 179) * k) % 179;
          i = j;
          j = k;
          k = m;
          l = ((53 * l) + 1) % 169;
          if (((l * m) % 64) >= 32) s = s + t;
          t = t / 2;
        }
      ranmar_params_u[ii] = s;
    }

  ranmar_params_c = __ranparm3;

  ranmar_params_i97 = 96;/* 97; */
  ranmar_params_j97 = 32;/* 33; */
}

#endif

#ifdef USE_MT19937
// This is the ``Mersenne Twister'' random number generator MT19937, which
// generates pseudorandom integers uniformly distributed in 0..(2^32 - 1)
// starting from any odd seed in 0..(2^32 - 1).  This version is a recode
// by Shawn Cokus (Cokus@math.washington.edu) on March 8, 1998 of a version by
// Takuji Nishimura (who had suggestions from Topher Cooper and Marc Rieffel in
// July-August 1997).
//
// Effectiveness of the recoding (on Goedel2.math.washington.edu, a DEC Alpha
// running OSF/1) using GCC -O3 as a compiler: before recoding: 51.6 sec. to
// generate 300 million random numbers; after recoding: 24.0 sec. for the same
// (i.e., 46.5% of original time), so speed is now about 12.5 million random
// number generations per second on this machine.
//
// According to the URL <http://www.math.keio.ac.jp/~matumoto/emt.html>
// (and paraphrasing a bit in places), the Mersenne Twister is ``designed
// with consideration of the flaws of various existing generators,'' has
// a period of 2^19937 - 1, gives a sequence that is 623-dimensionally
// equidistributed, and ``has passed many stringent tests, including the
// die-hard test of G. Marsaglia and the load test of P. Hellekalek and
// S. Wegenkittl.''  It is efficient in memory usage (typically using 2506
// to 5012 bytes of static data, depending on data type sizes, and the code
// is quite short as well).  It generates random numbers in batches of 624
// at a time, so the caching and pipelining of modern systems is exploited.
// It is also divide- and mod-free.
//
// This library is free software; you can redistribute it and/or modify it
// under the terms of the GNU Library General Public License as published by
// the Free Software Foundation (either version 2 of the License or, at your
// option, any later version).  This library is distributed in the hope that
// it will be useful, but WITHOUT ANY WARRANTY, without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See
// the GNU Library General Public License for more details.  You should have
// received a copy of the GNU Library General Public License along with this
// library; if not, write to the Free Software Foundation, Inc., 59 Temple
// Place, Suite 330, Boston, MA 02111-1307, USA.
//
// The code as Shawn received it included the following notice:
//
//   Copyright (C) 1997 Makoto Matsumoto and Takuji Nishimura.  When
//   you use this, send an e-mail to <matumoto@math.keio.ac.jp> with
//   an appropriate reference to your work.
//
// It would be nice to CC: <Cokus@math.washington.edu> when you write.
//

//
// uint32 must be an unsigned integer type capable of holding at least 32
// bits; exactly 32 should be fastest, but 64 is better on an Alpha with
// GCC at -O3 optimization so try your options and see what's best for you
//

#define _mt_N              (624)                 // length of state vector
#define _mt_M              (397)                 // a period parameter
#define _mt_K              (0x9908B0DFU)         // a magic constant
#define _mt_hiBit(u)       ((u) & 0x80000000U)   // mask all but highest   bit of u
#define _mt_loBit(u)       ((u) & 0x00000001U)   // mask all but lowest    bit of u
#define _mt_loBits(u)      ((u) & 0x7FFFFFFFU)   // mask     the highest   bit of u
#define _mt_mixBits(u, v)  (_mt_hiBit(u)|_mt_loBits(v))  // move hi bit of u to hi bit of v

static _mt_uint32   _mt_state[_mt_N+1];     // state vector + 1 extra to not violate ANSI C
_mt_uint32   *_mt_next;          // next random value is computed from here
int      _mt_left = -1;      // can *next++ this many times before reloading

//void seedMT(_mt_uint32 seed)
void RandSeed(_mt_uint32 seed)
 {
    //
    // We initialize state[0..(N-1)] via the generator
    //
    //   x_new = (69069 * x_old) mod 2^32
    //
    // from Line 15 of Table 1, p. 106, Sec. 3.3.4 of Knuth's
    // _The Art of Computer Programming_, Volume 2, 3rd ed.
    //
    // Notes (SJC): I do not know what the initial state requirements
    // of the Mersenne Twister are, but it seems this seeding generator
    // could be better.  It achieves the maximum period for its modulus
    // (2^30) iff x_initial is odd (p. 20-21, Sec. 3.2.1.2, Knuth); if
    // x_initial can be even, you have sequences like 0, 0, 0, ...;
    // 2^31, 2^31, 2^31, ...; 2^30, 2^30, 2^30, ...; 2^29, 2^29 + 2^31,
    // 2^29, 2^29 + 2^31, ..., etc. so I force seed to be odd below.
    //
    // Even if x_initial is odd, if x_initial is 1 mod 4 then
    //
    //   the          lowest bit of x is always 1,
    //   the  next-to-lowest bit of x is always 0,
    //   the 2nd-from-lowest bit of x alternates      ... 0 1 0 1 0 1 0 1 ... ,
    //   the 3rd-from-lowest bit of x 4-cycles        ... 0 1 1 0 0 1 1 0 ... ,
    //   the 4th-from-lowest bit of x has the 8-cycle ... 0 0 0 1 1 1 1 0 ... ,
    //    ...
    //
    // and if x_initial is 3 mod 4 then
    //
    //   the          lowest bit of x is always 1,
    //   the  next-to-lowest bit of x is always 1,
    //   the 2nd-from-lowest bit of x alternates      ... 0 1 0 1 0 1 0 1 ... ,
    //   the 3rd-from-lowest bit of x 4-cycles        ... 0 0 1 1 0 0 1 1 ... ,
    //   the 4th-from-lowest bit of x has the 8-cycle ... 0 0 1 1 1 1 0 0 ... ,
    //    ...
    //
    // The generator's potency (min. s>=0 with (69069-1)^s = 0 mod 2^32) is
    // 16, which seems to be alright by p. 25, Sec. 3.2.1.3 of Knuth.  It
    // also does well in the dimension 2..5 spectral tests, but it could be
    // better in dimension 6 (Line 15, Table 1, p. 106, Sec. 3.3.4, Knuth).
    //
    // Note that the random number user does not see the values generated
    // here directly since reloadMT() will always munge them first, so maybe
    // none of all of this matters.  In fact, the seed values made here could
    // even be extra-special desirable if the Mersenne Twister theory says
    // so-- that's why the only change I made is to restrict to odd seeds.
    //

    register _mt_uint32 x = (seed | 1U) & 0xFFFFFFFFU, *s = _mt_state;
    register int    j;

    for(_mt_left=0, *s++=x, j=_mt_N; --j;
        *s++ = (x*=69069U) & 0xFFFFFFFFU);
 }

_mt_uint32 reloadMT(void)
 {
    _mt_uint32 *p0=_mt_state, *p2=_mt_state+2, *pM=_mt_state+_mt_M, s0, s1;
    int    j;

//    if(_mt_left < -1) RandSeed(4357U);

    _mt_left=_mt_N-1, _mt_next=_mt_state+1;

    for(s0=_mt_state[0], s1=_mt_state[1], j=_mt_N-_mt_M+1; --j; s0=s1, s1=*p2++)
        *p0++ = *pM++ ^ (_mt_mixBits(s0, s1) >> 1) ^ (_mt_loBit(s1) ? _mt_K : 0U);

    for(pM=_mt_state, j=_mt_M; --j; s0=s1, s1=*p2++)
        *p0++ = *pM++ ^ (_mt_mixBits(s0, s1) >> 1) ^ (_mt_loBit(s1) ? _mt_K : 0U);

    s1=_mt_state[0], *p0 = *pM ^ (_mt_mixBits(s0, s1) >> 1) ^ (_mt_loBit(s1) ? _mt_K : 0U);
    s1 ^= (s1 >> 11);
    s1 ^= (s1 <<  7) & 0x9D2C5680U;
    s1 ^= (s1 << 15) & 0xEFC60000U;
    return(s1 ^ (s1 >> 18));
 }

#endif
