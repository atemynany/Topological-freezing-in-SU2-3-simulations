// ********************



// fields.cc



// ********************



#include "fields.hh"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "geometry.hh"
#include "linear_algebra.hh"
#include "ranlux.hh"



// ********************



// Allocates and frees memory for a gauge field (lattice size T * L^3).

void Gauge_Field_Alloc(double **gauge_field, int T, int L)
{
  //fprintf(stderr,
	 // "void Gauge_Field_Alloc(...   -->   Trying to allocate %ld M ...",
	//  T*L*L*L * 4 * 8 * sizeof(double) / 1000000);

  if((*gauge_field = (double *)malloc(T*L*L*L * 4 * 8 * sizeof(double))) ==
     NULL)
    {
      fprintf(stderr, "\nError: void Gauge_Field_Alloc(...\n");
      exit(EXIT_FAILURE);
    }

  // fprintf(stderr, " o.k.\n");
}

void Gauge_Field_Alloc_Mu_Fixed(double **gauge_field, int T, int L)
{
  //fprintf(stderr,
   // "void Mu_Fixed_Gauge_Field_Alloc(...   -->   Trying to allocate %ld M ...",
   // T*L*L*L * 8 * sizeof(double) / 1000000);

  if((*gauge_field = (double *)malloc(T*L*L*L * 8 * sizeof(double))) ==
     NULL)
    {
      fprintf(stderr, "\nError: void Mu_Fixed_Gauge_Field_Alloc(...\n");
      exit(EXIT_FAILURE);
    }

  //fprintf(stderr, " o.k.\n");
}

void Gauge_Field_Alloc_Timeslice(double **gauge_field, int L)
{
 // fprintf(stderr,
  //  "void Timeslice_Gauge_Field_Alloc(...   -->   Trying to allocate %ld M ...",
   // L*L*L * 4 * 8 * sizeof(double) / 1000000);

  if((*gauge_field = (double *)malloc(L*L*L * 4 * 8 * sizeof(double))) ==
     NULL)
    {
      fprintf(stderr, "\nError: void Gauge_Field_Alloc(...\n");
      exit(EXIT_FAILURE);
    }

  //fprintf(stderr, " o.k.\n");
}

void Gauge_Field_Free(double **gauge_field)
{
  if(*gauge_field == NULL)
    {
      fprintf(stderr, "Error: void Gauge_Field_Free(...\n");
      exit(EXIT_FAILURE);
    }

  free(*gauge_field);
  *gauge_field = NULL;
}



// ********************



// Copies a gauge field (lattice size T * L^3).

void Gauge_Field_Copy(double *gauge_field_dst, double *gauge_field_src, int T, int L)
{
  memcpy(gauge_field_dst, gauge_field_src, T*L*L*L * 4 * 8 * sizeof(double));
}

// Copies a timeslice of a gauge field (lattice size T * L^3).

void Gauge_Field_Copy_Timeslice(double *timeslice_gauge_field_dst, double *gauge_field_src, int T, int L, int timeslice)
{
  memcpy(timeslice_gauge_field_dst,
	 gauge_field_src + ggi(get_index(timeslice, 0, 0, 0, T, L), 0),
	 L*L*L * 4 * 8 * sizeof(double));
}



// ********************



// Generates a unit gauge field.

void Gauge_Field_Unity(double *gauge_field, int T, int L)
{
  int i1;
  int it, ix, iy, iz;

  for(it = 0; it < T; it++)
    {
      for(ix = 0; ix < L; ix++)
	{
	  for(iy = 0; iy < L; iy++)
	    {
	      for(iz = 0; iz < L; iz++)
		{
		  for(i1 = 0; i1 < 4; i1++)
		    {
		      int index = ggi(get_index(it, ix, iy, iz, T, L), i1);

		      cm_eq_id(gauge_field + index);
		    }
		}
	    }
	}
    }  
}

// Generates a random gauge field.

void Gauge_Field_Random(double *gauge_field, int T, int L)
{
  int i1;
  double SU2_1[8];

  int it, ix, iy, iz;

  for(it = 0; it < T; it++)
    {
      for(ix = 0; ix < L; ix++)
	{
	  for(iy = 0; iy < L; iy++)
	    {
	      for(iz = 0; iz < L; iz++)
		{
		  for(i1 = 0; i1 < 4; i1++)
		    {
		      double h[4];

		      double norm;

		      while(1)
			{
			  h[0] = 2.0 * DRand() - 1.0;
			  h[1] = 2.0 * DRand() - 1.0;
			  h[2] = 2.0 * DRand() - 1.0;
			  h[3] = 2.0 * DRand() - 1.0;

			  norm = pow(h[0], 2.0) + pow(h[1], 2.0) + pow(h[2], 2.0) + pow(h[3], 2.0);

			  if(norm >= 0.0000000001 && norm <= 1.0)
			    break;
			}

		      norm = 1.0 / sqrt(norm);

		      h[0] *= norm;
		      h[1] *= norm;
		      h[2] *= norm;
		      h[3] *= norm;

		      cm_from_h(SU2_1, h);

		      int index = ggi(get_index(it, ix, iy, iz, T, L), i1);

		      cm_eq_cm(gauge_field + index, SU2_1);
		    }
		}
	    }
	}
    }  
}



// ********************
