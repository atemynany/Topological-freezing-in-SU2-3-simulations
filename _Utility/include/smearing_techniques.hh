// ********************



// smearing_techniques.hh



// ********************



#ifndef __SMEARING_TECHNIQUES_HH__

#define __SMEARING_TECHNIQUES_HH__



// ********************



#include "linear_algebra.hh"



// ********************



// Computes fat links in time direction.

void Fat_Time_Links(double *gauge_field, double *smeared_gauge_field, int T, int L, double time_link_epsilon);


// Computes HYP links in time direction.

void HYP_Time_Links(double *gauge_field, double *smeared_gauge_field, int T, int L, double time_link_alpha1, double time_link_alpha2, double time_link_alpha3);


// Performs an APE smearing step.

void APE_Smearing_Step(double *smeared_gauge_field, int T, int L, double APE_smearing_alpha);
void APE_Smearing_all(double *smeared_gauge_field, int T, int L, double APE_smearing_alpha);


// Performs an APE smearing step on a given timeslice.

void APE_Smearing_Step_Timeslice(double *smeared_gauge_field, int T, int L, double APE_smearing_alpha, int timeslice);

// ********************



#endif



// ********************
