// ********************



// io.cc



// ********************



#include "io.hh"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "geometry.hh"
#include "linear_algebra.hh"



// ********************



// Reads/writes a gauge field configuration.

// for all it ...
//   for all ix ...
//     for all iy ...
//       for all iz ...
//         for all mu ...

// Data starts after "BEGIN".

// h[0], h[1], h[2], h[3] for each link.

void read_gauge_field(double *gauge_field, const char *filename, const int T, const int L)
{
  char string1[1000];


  FILE *fd;

  if((fd = fopen(filename, "r")) == NULL)
    {
      fprintf(stderr, "Error: void read_gauge_field(...!\n");
      fprintf(stderr, "  Filename does not exist.\n");
      exit(EXIT_FAILURE);
    }

  while(1)
    {
      if(fgets(string1, 995, fd) == NULL)
	{
	  fprintf(stderr, "Error: void read_gauge_field(...!\n");
	  fprintf(stderr, "  Wrong file format.\n");
	  exit(EXIT_FAILURE);
	}

      if(strcmp(string1, "BEGIN\n") == 0)
	break;
    }


  int it, ix, iy, iz, mu;

  for(it = 0; it < T; it++)
    {
      for(ix = 0; ix < L; ix++)
	{
	  for(iy = 0; iy < L; iy++)
	    {
	      for(iz = 0; iz < L; iz++)
		{
		  for(mu = 0; mu < 4; mu++)
		    {
  double h[4];

  if(fread(h, sizeof(double), 4, fd) != 4)
    {
      fprintf(stderr, "Error: void read_gauge_field(...!\n");
      fprintf(stderr, "  File too short!\n");
      exit(EXIT_FAILURE);
    }

  if(fabs(h[0]*h[0] + h[1]*h[1] + h[2]*h[2] + h[3]*h[3] - 1.0) > 0.0000000001)
    {
      fprintf(stderr, "Error: void read_gauge_field(...!\n");
      fprintf(stderr, "  it = %2d, ix = %2d, iy = %2d, iz = %2d, mu = %2d   -->   h^2 = %.6lf.\n", it, ix, iy, iz, mu, h[0]*h[0] + h[1]*h[1] + h[2]*h[2] + h[3]*h[3]);
      exit(EXIT_FAILURE);
    }

  int index = ggi(get_index(it, ix, iy, iz, T, L), mu);
  cm_from_h(gauge_field + index, h);
		    }
		}
	    }
	}
    }

  double h[4];

  if(fread(h, sizeof(double), 4, fd) != 0)
    {
      fprintf(stderr, "Error: void read_gauge_field(...!\n");
      fprintf(stderr, "  File too long!\n");
      exit(EXIT_FAILURE);
    }


  fclose(fd);
}


void read_gauge_field_header(double *gauge_field, const char *filename, const int T, const int L, char *header)
{
    char string1[1000];
    FILE *fd;

    if ((fd = fopen(filename, "r")) == NULL) {
        fprintf(stderr, "Error: void read_gauge_field(...)\n");
        fprintf(stderr, "  Filename does not exist.\n");
        exit(EXIT_FAILURE);
    }

    // Header speichern
    header[0] = '\0';  // Sicherstellen, dass der Header leer ist
    while (1) {
        if (fgets(string1, 995, fd) == NULL) {
            fprintf(stderr, "Error: void read_gauge_field(...)\n");
            fprintf(stderr, "  Wrong file format.\n");
            exit(EXIT_FAILURE);
        }

        // Header speichern, solange es nicht "BEGIN\n" ist
        if (strcmp(string1, "BEGIN\n") == 0) {
            strcat(header, "BEGIN\n");
            break;
        } else {
            strcat(header, string1);
        }
    }

    int it, ix, iy, iz, mu;

    for (it = 0; it < T; it++) {
        for (ix = 0; ix < L; ix++) {
            for (iy = 0; iy < L; iy++) {
                for (iz = 0; iz < L; iz++) {
                    for (mu = 0; mu < 4; mu++) {
                        double h[4];

                        if (fread(h, sizeof(double), 4, fd) != 4) {
                            fprintf(stderr, "Error: void read_gauge_field(...)\n");
                            fprintf(stderr, "  File too short!\n");
                            exit(EXIT_FAILURE);
                        }

                        if (fabs(h[0] * h[0] + h[1] * h[1] + h[2] * h[2] + h[3] * h[3] - 1.0) > 1e-10) {
                            fprintf(stderr, "Error: void read_gauge_field(...)\n");
                            fprintf(stderr, "  it = %2d, ix = %2d, iy = %2d, iz = %2d, mu = %2d   -->   h^2 = %.6lf.\n", it, ix, iy, iz, mu, h[0] * h[0] + h[1] * h[1] + h[2] * h[2] + h[3] * h[3]);
                            exit(EXIT_FAILURE);
                        }

                        int index = ggi(get_index(it, ix, iy, iz, T, L), mu);
                        cm_from_h(gauge_field + index, h);
                    }
                }
            }
        }
    }

    double h[4];

    if (fread(h, sizeof(double), 4, fd) != 0) {
        fprintf(stderr, "Error: void read_gauge_field(...)\n");
        fprintf(stderr, "  File too long!\n");
        exit(EXIT_FAILURE);
    }

    fclose(fd);
}



void write_gauge_field(double *gauge_field, const char *filename, const int T, const int L, const char *header)
{
  char string1[1000];


  FILE *fd;

  if((fd = fopen(filename, "w")) == NULL)
    {
      fprintf(stderr, "Error: void write_gauge_field(...!\n");
      exit(EXIT_FAILURE);
    }


  fprintf(fd, "# %s\n", header);
  fprintf(fd, "BEGIN\n");


  int it, ix, iy, iz, mu;

  for(it = 0; it < T; it++)
    {
      for(ix = 0; ix < L; ix++)
	{
	  for(iy = 0; iy < L; iy++)
	    {
	      for(iz = 0; iz < L; iz++)
		{
		  for(mu = 0; mu < 4; mu++)
		    {
  double h[4];
  h_from_cm(h, gauge_field + ggi(get_index(it, ix, iy, iz, T, L), mu));

  if(fabs(h[0]*h[0] + h[1]*h[1] + h[2]*h[2] + h[3]*h[3] - 1.0) > 0.0000000001)
    {
      fprintf(stderr, "Error: void write_gauge_field(...!\n");
      fprintf(stderr, "  it = %2d, ix = %2d, iy = %2d, iz = %2d, mu = %2d   -->   h^2 = %.6lf.\n", it, ix, iy, iz, mu, h[0]*h[0] + h[1]*h[1] + h[2]*h[2] + h[3]*h[3]);
      exit(EXIT_FAILURE);
    }

  if(fwrite(h, sizeof(double), 4, fd) != 4)
    {
      fprintf(stderr, "Error: void write_gauge_field(...!\n");
      exit(EXIT_FAILURE);
    }
		    }
		}
	    }
	}
    }


  fclose(fd);
}



// ********************
