#ifndef INSTANTON_HH
#define INSTANTON_HH

#include "../_Utility/include/linear_algebra.hh"
#include "../_Utility/include/geometry.hh"
#include "../_Utility/include/fields.hh"

void insert_instanton(double *gauge_field, double rho, double a, int L, int T, bool is_anti);

double instanton_topological_charge(double rho, double a, int L, int T,
                                    int smearing_steps, double smearing_alpha, bool is_anti);

void instanton_Q_vs_smearing(double rho, double a, int L, int T,
                             int max_smear_steps, double smearing_alpha, bool is_anti);

#endif
