#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// --- USER CONFIGURATION SECTION ---

#define DIM 4
#define L_T 16               // Lattice Time extent
#define L_S 16               // Lattice Spatial extent
#define VOL (L_T*L_S*L_S*L_S)
#define SUN 2                // Set to 3 for SU(3)
#define PI 3.14159265358979323846

// ----------------------------------

typedef struct { double r, i; } complex;
typedef struct { complex e[SUN][SUN]; } sun_mat;

// Global Gauge Field
sun_mat *pu[VOL][DIM]; 
int neib[VOL][2*DIM];