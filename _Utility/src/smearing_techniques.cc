// ********************


// APE_smearing.cc


// erweitern für timelinks, um die topologische Ladung zu berechnen



// ********************



#include "smearing_techniques.hh"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fields.hh"
#include "geometry.hh"



// ********************



// Performs an APE smearing step.

void APE_Smearing_all(double *smeared_gauge_field, int T, int L, double APE_smearing_alpha)
{
  int it, ix, iy, iz;
  double M1[18], M2[18];


  double *smeared_gauge_field_old;
  Gauge_Field_Alloc(&smeared_gauge_field_old, T, L);
  Gauge_Field_Copy(smeared_gauge_field_old, smeared_gauge_field, T, L);


  for(it = 0; it < T; it++)
    {
      for(ix = 0; ix < L; ix++)
	{
	  for(iy = 0; iy < L; iy++)
	    {
	      for(iz = 0; iz < L; iz++)
		{
		  int index;

		  int index_mt_1, index_mt_2, index_mt_3;
		  int index_pt_1, index_pt_2, index_pt_3;
		  int index_mx_1, index_mx_2, index_mx_3;
		  int index_px_1, index_px_2, index_px_3;
		  int index_my_1, index_my_2, index_my_3;
		  int index_py_1, index_py_2, index_py_3;
		  int index_mz_1, index_mz_2, index_mz_3;
		  int index_pz_1, index_pz_2, index_pz_3;
		  
		
		  double *U;


  // **********
  // Links in x-direction.
  // **********


  index = ggi(get_index(it, ix, iy, iz, T, L), 1);

  index_mt_1 = ggi(get_index(it-1, ix, iy, iz, T, L), 0);
  index_mt_2 = ggi(get_index(it-1, ix, iy, iz, T, L), 1);
  index_mt_3 = ggi(get_index(it-1, ix+1, iy, iz, T, L), 0);

  index_pt_1 = ggi(get_index(it, ix, iy, iz, T, L), 0);
  index_pt_2 = ggi(get_index(it+1, ix, iy, iz, T, L), 1);
  index_pt_3 = ggi(get_index(it, ix+1, iy, iz, T, L), 0);

  index_my_1 = ggi(get_index(it, ix, iy-1, iz, T, L), 2);
  index_my_2 = ggi(get_index(it, ix, iy-1, iz, T, L), 1);
  index_my_3 = ggi(get_index(it, ix+1, iy-1, iz, T, L), 2);

  index_py_1 = ggi(get_index(it, ix, iy, iz, T, L), 2);
  index_py_2 = ggi(get_index(it, ix, iy+1, iz, T, L), 1);
  index_py_3 = ggi(get_index(it, ix+1, iy, iz, T, L), 2);

  index_mz_1 = ggi(get_index(it, ix, iy, iz-1, T, L), 3);
  index_mz_2 = ggi(get_index(it, ix, iy, iz-1, T, L), 1);
  index_mz_3 = ggi(get_index(it, ix+1, iy, iz-1, T, L), 3);

  index_pz_1 = ggi(get_index(it, ix, iy, iz, T, L), 3);
  index_pz_2 = ggi(get_index(it, ix, iy, iz+1, T, L), 1);
  index_pz_3 = ggi(get_index(it, ix+1, iy, iz, T, L), 3);

if (index_mt_1>=0 && index_mt_2>=0 && index_mt_3 && index_pt_2>=0){

  U = smeared_gauge_field + index;	
  cm_eq_zero(U);

// negative t-direction


cm_eq_cm_ti_cm(M1, smeared_gauge_field_old + index_mt_2,
	smeared_gauge_field_old + index_mt_3);

cm_eq_cm_dag_ti_cm(M2, smeared_gauge_field_old + index_mt_1, M1);

cm_pl_eq_cm(U, M2);


// positive t-direction

cm_eq_cm_ti_cm_dag(M1, smeared_gauge_field_old + index_pt_2,
	smeared_gauge_field_old + index_pt_3);

cm_eq_cm_ti_cm(M2, smeared_gauge_field_old + index_pt_1, M1);

cm_pl_eq_cm(U, M2);


  // negative y-direction

  cm_eq_cm_ti_cm(M1, smeared_gauge_field_old + index_my_2,
	     smeared_gauge_field_old + index_my_3);

  cm_eq_cm_dag_ti_cm(M2, smeared_gauge_field_old + index_my_1, M1);

  cm_pl_eq_cm(U, M2);


  // positive y-direction

  cm_eq_cm_ti_cm_dag(M1, smeared_gauge_field_old + index_py_2,
	     smeared_gauge_field_old + index_py_3);

  cm_eq_cm_ti_cm(M2, smeared_gauge_field_old + index_py_1, M1);

  cm_pl_eq_cm(U, M2);


  // negative z-direction

  cm_eq_cm_ti_cm(M1, smeared_gauge_field_old + index_mz_2,
	     smeared_gauge_field_old + index_mz_3);

  cm_eq_cm_dag_ti_cm(M2, smeared_gauge_field_old + index_mz_1, M1);

  cm_pl_eq_cm(U, M2);


  // positive z-direction

  cm_eq_cm_ti_cm_dag(M1, smeared_gauge_field_old + index_pz_2,
	     smeared_gauge_field_old + index_pz_3);

  cm_eq_cm_ti_cm(M2, smeared_gauge_field_old + index_pz_1, M1);

  cm_pl_eq_cm(U, M2);


  cm_ti_eq_re(U, APE_smearing_alpha);


  // center

  cm_pl_eq_cm(U, smeared_gauge_field_old + index); // U+\alpha C


  // Projection to SU(2).
  cm_proj(U);
}


  // **********
  // Links in y-direction.
  // **********


  index = ggi(get_index(it, ix, iy, iz, T, L), 2);

  index_mt_1 = ggi(get_index(it-1, ix, iy, iz, T, L), 0);
  index_mt_2 = ggi(get_index(it-1, ix, iy, iz, T, L), 2);
  index_mt_3 = ggi(get_index(it-1, ix, iy+1, iz, T, L), 0);

  index_pt_1 = ggi(get_index(it, ix, iy, iz, T, L), 0);
  index_pt_2 = ggi(get_index(it+1, ix, iy, iz, T, L), 2);
  index_pt_3 = ggi(get_index(it, ix, iy+1, iz, T, L), 0);

  index_mx_1 = ggi(get_index(it, ix-1, iy, iz, T, L), 1);
  index_mx_2 = ggi(get_index(it, ix-1, iy, iz, T, L), 2);
  index_mx_3 = ggi(get_index(it, ix-1, iy+1, iz, T, L), 1);

  index_px_1 = ggi(get_index(it, ix, iy, iz, T, L), 1);
  index_px_2 = ggi(get_index(it, ix+1, iy, iz, T, L), 2);
  index_px_3 = ggi(get_index(it, ix, iy+1, iz, T, L), 1);

  index_mz_1 = ggi(get_index(it, ix, iy, iz-1, T, L), 3);
  index_mz_2 = ggi(get_index(it, ix, iy, iz-1, T, L), 2);
  index_mz_3 = ggi(get_index(it, ix, iy+1, iz-1, T, L), 3);

  index_pz_1 = ggi(get_index(it, ix, iy, iz, T, L), 3);
  index_pz_2 = ggi(get_index(it, ix, iy, iz+1, T, L), 2);
  index_pz_3 = ggi(get_index(it, ix, iy+1, iz, T, L), 3);

if (index_mt_1>=0 && index_mt_2>=0 && index_mt_3 && index_pt_2>=0){

  U = smeared_gauge_field + index;
  cm_eq_zero(U);

// negative t-direction


cm_eq_cm_ti_cm(M1, smeared_gauge_field_old + index_mt_2,
	smeared_gauge_field_old + index_mt_3);

cm_eq_cm_dag_ti_cm(M2, smeared_gauge_field_old + index_mt_1, M1);

cm_pl_eq_cm(U, M2);


// positive t-direction

cm_eq_cm_ti_cm_dag(M1, smeared_gauge_field_old + index_pt_2,
	smeared_gauge_field_old + index_pt_3);

cm_eq_cm_ti_cm(M2, smeared_gauge_field_old + index_pt_1, M1);

cm_pl_eq_cm(U, M2);


  // negative x-direction

  cm_eq_cm_ti_cm(M1, smeared_gauge_field_old + index_mx_2,
	     smeared_gauge_field_old + index_mx_3);

  cm_eq_cm_dag_ti_cm(M2, smeared_gauge_field_old + index_mx_1, M1);

  cm_pl_eq_cm(U, M2);


  // positive x-direction

  cm_eq_cm_ti_cm_dag(M1, smeared_gauge_field_old + index_px_2,
	     smeared_gauge_field_old + index_px_3);

  cm_eq_cm_ti_cm(M2, smeared_gauge_field_old + index_px_1, M1);

  cm_pl_eq_cm(U, M2);


  // negative z-direction

  cm_eq_cm_ti_cm(M1, smeared_gauge_field_old + index_mz_2,
	     smeared_gauge_field_old + index_mz_3);

  cm_eq_cm_dag_ti_cm(M2, smeared_gauge_field_old + index_mz_1, M1);

  cm_pl_eq_cm(U, M2);


  // positive z-direction

  cm_eq_cm_ti_cm_dag(M1, smeared_gauge_field_old + index_pz_2,
	     smeared_gauge_field_old + index_pz_3);

  cm_eq_cm_ti_cm(M2, smeared_gauge_field_old + index_pz_1, M1);

  cm_pl_eq_cm(U, M2);


  cm_ti_eq_re(U, APE_smearing_alpha);


  // center - Frage: was macht das?

  cm_pl_eq_cm(U, smeared_gauge_field_old + index);


  // Projection to SU(3). Frage: brauche ich das? (Ich nutze SU(2)??)
  cm_proj(U);
}


  // **********
  // Links in z-direction.
  // **********


  index = ggi(get_index(it, ix, iy, iz, T, L), 3);

  index_mt_1 = ggi(get_index(it-1, ix, iy, iz, T, L), 0);
  index_mt_2 = ggi(get_index(it-1, ix, iy, iz, T, L), 3);
  index_mt_3 = ggi(get_index(it-1, ix, iy, iz+1, T, L), 0);

  index_pt_1 = ggi(get_index(it, ix, iy, iz, T, L), 0);
  index_pt_2 = ggi(get_index(it+1, ix, iy, iz, T, L), 3);
  index_pt_3 = ggi(get_index(it, ix, iy, iz+1, T, L), 0);

  index_mx_1 = ggi(get_index(it, ix-1, iy, iz, T, L), 1);
  index_mx_2 = ggi(get_index(it, ix-1, iy, iz, T, L), 3);
  index_mx_3 = ggi(get_index(it, ix-1, iy, iz+1, T, L), 1);

  index_px_1 = ggi(get_index(it, ix, iy, iz, T, L), 1);
  index_px_2 = ggi(get_index(it, ix+1, iy, iz, T, L), 3);
  index_px_3 = ggi(get_index(it, ix, iy, iz+1, T, L), 1);

  index_my_1 = ggi(get_index(it, ix, iy-1, iz, T, L), 2);
  index_my_2 = ggi(get_index(it, ix, iy-1, iz, T, L), 3);
  index_my_3 = ggi(get_index(it, ix, iy-1, iz+1, T, L), 2);

  index_py_1 = ggi(get_index(it, ix, iy, iz, T, L), 2);
  index_py_2 = ggi(get_index(it, ix, iy+1, iz, T, L), 3);
  index_py_3 = ggi(get_index(it, ix, iy, iz+1, T, L), 2);

if (index_mt_1>=0 && index_mt_2>=0 && index_mt_3 && index_pt_2>=0){

  U = smeared_gauge_field + index; 
  cm_eq_zero(U);

// negative t-direction


cm_eq_cm_ti_cm(M1, smeared_gauge_field_old + index_mt_2,
	smeared_gauge_field_old + index_mt_3);

cm_eq_cm_dag_ti_cm(M2, smeared_gauge_field_old + index_mt_1, M1);

cm_pl_eq_cm(U, M2);


// positive t-direction

cm_eq_cm_ti_cm_dag(M1, smeared_gauge_field_old + index_pt_2,
	smeared_gauge_field_old + index_pt_3);

cm_eq_cm_ti_cm(M2, smeared_gauge_field_old + index_pt_1, M1);

cm_pl_eq_cm(U, M2);


  // negative x-direction

  cm_eq_cm_ti_cm(M1, smeared_gauge_field_old + index_mx_2,
	     smeared_gauge_field_old + index_mx_3);

  cm_eq_cm_dag_ti_cm(M2, smeared_gauge_field_old + index_mx_1, M1);

  cm_pl_eq_cm(U, M2);


  // positive x-direction

  cm_eq_cm_ti_cm_dag(M1, smeared_gauge_field_old + index_px_2,
	     smeared_gauge_field_old + index_px_3);

  cm_eq_cm_ti_cm(M2, smeared_gauge_field_old + index_px_1, M1);

  cm_pl_eq_cm(U, M2);


  // negative y-direction

  cm_eq_cm_ti_cm(M1, smeared_gauge_field_old + index_my_2,
	     smeared_gauge_field_old + index_my_3);

  cm_eq_cm_dag_ti_cm(M2, smeared_gauge_field_old + index_my_1, M1);

  cm_pl_eq_cm(U, M2);


  // positive y-direction

  cm_eq_cm_ti_cm_dag(M1, smeared_gauge_field_old + index_py_2,
	     smeared_gauge_field_old + index_py_3);

  cm_eq_cm_ti_cm(M2, smeared_gauge_field_old + index_py_1, M1);

  cm_pl_eq_cm(U, M2);


  cm_ti_eq_re(U, APE_smearing_alpha);


  // center

  cm_pl_eq_cm(U, smeared_gauge_field_old + index);


  // Projection to SU(3).
  cm_proj(U);
 }

 
  // **********
  // Links in t-direction.
  // **********


  index = ggi(get_index(it, ix, iy, iz, T, L), 0);

  index_mx_1 = ggi(get_index(it, ix-1, iy, iz, T, L), 1);
  index_mx_2 = ggi(get_index(it, ix-1, iy, iz, T, L), 0);
  index_mx_3 = ggi(get_index(it+1, ix-1, iy, iz, T, L), 1);

  index_px_1 = ggi(get_index(it, ix, iy, iz, T, L), 1);
  index_px_2 = ggi(get_index(it, ix+1, iy, iz, T, L), 0);
  index_px_3 = ggi(get_index(it+1, ix, iy, iz, T, L), 1);

  index_my_1 = ggi(get_index(it, ix, iy-1, iz, T, L), 2);
  index_my_2 = ggi(get_index(it, ix, iy-1, iz, T, L), 0);
  index_my_3 = ggi(get_index(it+1, ix, iy-1, iz, T, L), 2);

  index_py_1 = ggi(get_index(it, ix, iy, iz, T, L), 2);
  index_py_2 = ggi(get_index(it, ix, iy+1, iz, T, L), 0);
  index_py_3 = ggi(get_index(it+1, ix, iy, iz, T, L), 2);

  index_mz_1 = ggi(get_index(it, ix, iy, iz-1, T, L), 3);
  index_mz_2 = ggi(get_index(it, ix, iy, iz-1, T, L), 0);
  index_mz_3 = ggi(get_index(it+1, ix, iy, iz-1, T, L), 3);

  index_pz_1 = ggi(get_index(it, ix, iy, iz, T, L), 3);
  index_pz_2 = ggi(get_index(it, ix, iy, iz+1, T, L), 0);
  index_pz_3 = ggi(get_index(it+1, ix, iy, iz, T, L), 3);


  if (index >= 0 && index_mx_3>=0 && index_px_3>=0 && index_my_3>=0 && index_py_3>=0 && index_mz_3>=0 && index_pz_3>=0){

  U = smeared_gauge_field + index; 
  cm_eq_zero(U);


  // negative x-direction


  cm_eq_cm_ti_cm(M1, smeared_gauge_field_old + index_mx_2,
	     smeared_gauge_field_old + index_mx_3);

  cm_eq_cm_dag_ti_cm(M2, smeared_gauge_field_old + index_mx_1, M1);

  cm_pl_eq_cm(U, M2);


  // positive x-direction

  cm_eq_cm_ti_cm_dag(M1, smeared_gauge_field_old + index_px_2,
	     smeared_gauge_field_old + index_px_3);

  cm_eq_cm_ti_cm(M2, smeared_gauge_field_old + index_px_1, M1);

  cm_pl_eq_cm(U, M2);


  // negative y-direction

  cm_eq_cm_ti_cm(M1, smeared_gauge_field_old + index_my_2,
	     smeared_gauge_field_old + index_my_3);

  cm_eq_cm_dag_ti_cm(M2, smeared_gauge_field_old + index_my_1, M1);

  cm_pl_eq_cm(U, M2);


  // positive y-direction

  cm_eq_cm_ti_cm_dag(M1, smeared_gauge_field_old + index_py_2,
	     smeared_gauge_field_old + index_py_3);

  cm_eq_cm_ti_cm(M2, smeared_gauge_field_old + index_py_1, M1);

  cm_pl_eq_cm(U, M2);


  cm_ti_eq_re(U, APE_smearing_alpha);


  // negative z-direction

cm_eq_cm_ti_cm(M1, smeared_gauge_field_old + index_mz_2,
	smeared_gauge_field_old + index_mz_3);

cm_eq_cm_dag_ti_cm(M2, smeared_gauge_field_old + index_mz_1, M1);

cm_pl_eq_cm(U, M2);


// positive z-direction

cm_eq_cm_ti_cm_dag(M1, smeared_gauge_field_old + index_pz_2,
	smeared_gauge_field_old + index_pz_3);

cm_eq_cm_ti_cm(M2, smeared_gauge_field_old + index_pz_1, M1);

cm_pl_eq_cm(U, M2);


  // center

  cm_pl_eq_cm(U, smeared_gauge_field_old + index);


  // Projection to SU(3).
  cm_proj(U);
}


		}
	    }
	}
    }


  Gauge_Field_Free(&smeared_gauge_field_old);
}



// ********************
