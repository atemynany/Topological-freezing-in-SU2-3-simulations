//
//  topologicalcharge.c
//  glueballs
//
//  Created by Carolin Riehl on 14.05.20.
//  Copyright © 2020 Carolin Riehl. All rights reserved.
//

#define TOPCHARGE_C

#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include"gauge.h"
#include"headers.h"
#include"modules.h"





static void U_plpl(sun_mat *u, int ii,int jj, int kk)
{
    int n;
    sun_mat *un;
    
    
    un=malloc(3*sizeof(sun_mat));
    
    // Compute Plaquette U_mu nu.
    un[0]=*pu[ii][jj];
    n=neib[ii][jj];
    un[1]=*pu[n][kk];
    sun_mul(un[2],un[0],un[1]);
    n=neib[ii][kk];
    sun_dag(un[0],*pu[n][jj]);
    sun_mul(un[1],un[2],un[0]);
    sun_dag(un[2],*pu[ii][kk]);
    sun_mul(un[0],un[1],un[2]);
    
    u[0]=un[0];
    free(un);
}

static void U_mipl(sun_mat *u, int ii,int jj, int kk)
{
    int n;
    sun_mat *un;
    
    
    un=malloc(3*sizeof(sun_mat));
    
    // Compute Plaquette U_mu nu.
    n=neib[ii][jj+DIM];
    sun_dag(un[0],*pu[n][jj]);
    un[1]=*pu[n][kk];
    sun_mul(un[2],un[0],un[1]);
    
    n=neib[n][kk];
    
    sun_mul(un[1],un[2],*pu[n][jj]);
    sun_dag(un[2],*pu[ii][kk]);
    sun_mul(un[0],un[1],un[2]);
    
    u[0]=un[0];
    free(un);
}

static void U_mimi(sun_mat *u, int ii,int jj, int kk)
{
    int n;
    sun_mat *un;
    
    
    un=malloc(3*sizeof(sun_mat));
    
    // Compute Plaquette U_ -mu -nu.
    n=neib[ii][jj+DIM];
    sun_dag(un[0],*pu[n][jj]);
    n=neib[n][kk+DIM];
    sun_mul_dag(un[2],un[0],*pu[n][kk]);
    sun_mul(un[1],un[2],*pu[n][jj]);
    n=neib[ii][kk+DIM];
    sun_mul(un[0],un[1],*pu[n][kk]);
    
    u[0]=un[0];
    free(un);
}

static void U_plmi(sun_mat *u, int ii,int jj, int kk)
{
    int n;
    sun_mat *un;
    
    
    un=malloc(3*sizeof(sun_mat));
    
    // Compute Plaquette U_ +mu -nu.
    un[0]=*pu[ii][jj];
    n=neib[ii][jj];
    n=neib[n][kk+DIM];
    sun_mul_dag(un[2],un[0],*pu[n][kk]);
    n=neib[ii][kk+DIM];
    sun_dag(un[0],*pu[n][jj]);
    sun_mul(un[1],un[2],un[0]);
    sun_mul(un[0],un[1],*pu[n][kk]);
    
    u[0]=un[0];
    free(un);
}


void clover(sun_mat *clov, int x, int d1, int d2)
{
    double tr;
     sun_mat *un, *p;
     
     //checkpoint("compute_clover");
     
     un=malloc(3*sizeof(sun_mat));
     p=malloc(2*sizeof(sun_mat));
     sun_zero(p[0]);
     sun_zero(p[1]);
     

      
      
      
      U_plpl(un, x, d1, d2);
      sun_add(p[1],un[0],p[0]);
     
      U_mipl(un, x, d2, d1);
      sun_add(p[0],un[0],p[1]);
      U_mimi(un, x, d1, d2);
      sun_add(p[1],un[0],p[0]);
      U_plmi(un, x, d2, d1);
      sun_add(p[0],un[0],p[1]);

      un[0]=p[0];
     
    // Im(C_{mu nu}).
    sun_dag(un[1], un[0]);
    sun_sub(p[0], un[0], un[1]);
    
    
    sun_dble_div(p[0], 8.);
    *clov=p[0];
     

}


double compute_clover_products(int x,int d1,int d2, int d3, int d4)
{
    
   double tr;
    sun_mat *un, *unn, *p;
    
    //checkpoint("compute_plaquette");
    
    un=malloc(3*sizeof(sun_mat));
    unn=malloc(3*sizeof(sun_mat));
    p=malloc(2*sizeof(sun_mat));
    sun_zero(p[0]);
    sun_zero(p[1]);
    

     
     
     
     U_plpl(un, x, d1, d2);
    
     sun_add(p[1],un[0],p[0]);
    
     U_mipl(un, x, d2, d1);
     sun_add(p[0],un[0],p[1]);
     U_mimi(un, x, d1, d2);
     sun_add(p[1],un[0],p[0]);
     U_plmi(un, x, d2, d1);
     sun_add(p[0],un[0],p[1]);

     un[0]=p[0];
     sun_dble_div(un[0], 4.);
    
    

     
     sun_zero(p[0]);
     sun_zero(p[1]);
         
      U_plpl(unn, x, d3, d4);
    
      sun_add(p[1],unn[0],p[0]);
      U_mipl(unn, x, d4, d3);
      sun_add(p[0],unn[0],p[1]);
      U_mimi(unn, x, d3, d4);
      sun_add(p[1],unn[0],p[0]);
      U_plmi(unn, x, d4, d3);
      sun_add(p[0],unn[0],p[1]);
     unn[0]=p[0];
     sun_dble_div(unn[0], 4.);
    
    
    // compute trace of product of imaginary part of plaquettes C_{mu nu} C_{rho sigma}.

    // compute Im (C) * Im (C).
    sun_zero(p[0]);
    sun_zero(p[1]);
    sun_zero(un[1]);
    sun_zero(unn[1]);
    sun_zero(un[2]);
    
    // Im(C_{mu nu}).
    sun_dag(un[1], un[0]);
    sun_sub(p[0], un[0], un[1]);
    

    // Im(C_{rho sigma}).
    sun_dag(unn[1], unn[0]);
    sun_sub(p[1], unn[0], unn[1]);
    
    sun_dble_div(p[0], 2.);
    sun_dble_div(p[1], 2.);
    // ImC*ImC.
    sun_mul(un[2], p[0], p[1]);
    sun_trace(tr, un[2]);
    

    
    
    free(un);
    free(unn);
    free(p);
    

   return tr;
}

void meas_topologicalcharge(double *q)
{
    int mu,nu,rho,sig;
    int epsilon;
    int n=0;
    double tr;
    double tstart1,tend1;
    

    checkpoint("meas_topologicalcharge");
    
    //********************* Explcit formula (slower version) *******************************************
    /*
    tstart1=get_time();

    *q=0.;
    
int count=0;
        for (mu=0; mu < DIM; mu++){
            for (nu=0; nu<DIM; nu++) {
                for (rho=0; rho< DIM; rho++) {
                    for (sig=0; sig<DIM; sig++) {
                        
                        epsilon= ( (nu-mu)*(rho-mu)*(sig-mu)*(rho-nu)*(sig-nu)*(sig-rho) )/12;
                        
                        
                        if (epsilon !=0) {
                            //printf("(mu,nu,rho,sigma) = (%d,%d,%d,%d) -> eps = %d\n", mu,nu,rho,sig, epsilon);
                            count++;
                        for(n=0; n< VOL; n++){
                            tr=compute_clover_products(n, mu, nu, rho, sig);
                            
                            //printf("Tr[ Im P_{mu nu} Im P_{rho sigma} ] = %g\n", tr);
                            *q+= epsilon*tr;
                            
                            //printf("... Q =  %e\n", *q);
                        }
                        
                    }
                }
            }
        }
    }
    *q/=(32.*M_PI*M_PI);
    printf("Q =  %e\n", *q);
    tend1=get_time();
    printf("Zeit 1 = %e count %d\n", tend1-tstart1, count);
    logging("Q1 %d %d %e\n", inlmeas.isms, 10000, *q);
    */
    //*********************************************************************************************
    
  
    tstart1=get_time();
    double q2=0;
    
 
    //********************* Reduced formula (faster version)  *******************************************
    
    for(n=0; n< VOL; n++){
        q2+= -compute_clover_products(n, 0, 1, 3, 2) -compute_clover_products(n, 1, 2, 3, 0) -compute_clover_products(n, 2, 0, 3, 1);
    }
 
    
    
    tend1=get_time();
    //printf("Time for computation of topological charge = %e\n", tend1-tstart1);
    
    
    *q=q2 /(4.*M_PI*M_PI );
    
    
    
}

