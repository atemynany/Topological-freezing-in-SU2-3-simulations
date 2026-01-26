// ********************



// fields.hh



// ********************



#ifndef __FIELDS_HH__

#define __FIELDS_HH__



// ********************



#include <stdio.h>



// ********************



// Allocates and frees memory for a gauge field (lattice size T * L^3).

void Gauge_Field_Alloc(double **gauge_field, int T, int L);

void Gauge_Field_Alloc_Mu_Fixed(double **gauge_field, int T, int L);

void Gauge_Field_Alloc_Timeslice(double **gauge_field, int L);

void Gauge_Field_Free(double **gauge_field);



// Copies a gauge field (lattice size T * L^3).

void Gauge_Field_Copy(double *gauge_field_dst, double *gauge_field_src, int T, int L);

// Copies a timeslice of a gauge field (lattice size T * L^3).

void Gauge_Field_Copy_Timeslice(double *timeslice_gauge_field_dst, double *gauge_field_src, int T, int L, int timeslice);



// Generates a unit gauge field.

void Gauge_Field_Unity(double *gauge_field, int T, int L);

// Generates a random gauge field.

void Gauge_Field_Random(double *gauge_field, int T, int L);



// ********************



#endif



// ********************
