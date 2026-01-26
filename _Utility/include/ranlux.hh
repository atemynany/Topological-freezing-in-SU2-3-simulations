// ********************



// ranlux.hh

// Author: Marc Wagner
//   (this is essentially a compilation of functions written by Martin Luescher
//    for random number generation [ranlux-3.3])
// Date: January 2010



// ********************



#ifndef __RANLUX_HH__

#define __RANLUX_HH__



// ********************



// performs the tests from Luescher's "testlx.c" program
void testlx();


// initializes the ranlxd random number generator at highest quality with seed seed
void InitializeRand(int seed);


// returns a uniformly chosen integer i with min <= i <= max
int IRand(int min, int max);


// returns a uniformly chosen double d in [0.0,1.0)
double DRand();



// ********************



#endif



// ********************
