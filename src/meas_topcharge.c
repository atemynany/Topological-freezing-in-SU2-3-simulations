//
//  main_topcharge.c
///// Originally /////
//  Created by Carolin Riehl on 15.05.20.
//  Copyright © 2020 Carolin Riehl. All rights reserved.
//
//// Edited by Alexander de Barros Noll on 2024-06-10. ////

#define MAIN_C

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include"../../Utility/include/ranlxd.h"
#include"headers.h"
#include"modules.h"

static int sconf,fconf, binsize;
char config_dir[NAME_SIZE];
char base[NAME_SIZE];
static int init=0;



static void read_inp_meas(int *iseed,char *dir,int argc,char *argv[])
{
    FILE *inf=NULL,*ftest=NULL;
    int i,ii,ifail,id;
    
    int ismt,isms;
    double beta,ismpt,ismps;
    
    ftest=fopen("SETUP_ERRORS","r");
    if(ftest!=NULL)
    {
        fclose(ftest);
        remove("SETUP_ERRORS");
    }
    sprintf(LOG_FILE,"SETUP_ERRORS");
    sprintf(OUT_FILE, "OUT");
    
    for(i=1,ii=0;i<argc;i++)
        if (strcmp(argv[i],"-i")==0)
            ii=i+1;
    error(ii==0,"read_inp_meas [meas.c]",
          "Syntax: meas -i <input file> [-c <input config>]");
    inf=fopen(argv[ii],"r");
    error(inf==NULL,"read_inp_meas [meas.c]",
          "Unable to open input file!");
    
    ifail=fscanf(inf,"id      %d\n",&id);
    ifail=fscanf(inf,"odir    %[^'\n']\n",dir);
    ifail=fscanf(inf,"beta    %lf\n",&beta);
    ifail=fscanf(inf,"seed    %d\n",iseed);
    ifail=fscanf(inf,"sconf   %d\n",&sconf);
    ifail=fscanf(inf,"fconf   %d\n",&fconf);
    ifail=fscanf(inf,"smear_t %d %lf\n",&ismt,&ismpt);
    ifail=fscanf(inf,"smear_s %d %lf\n",&isms,&ismps);

    
    runp.id=id;
    runp.beta=beta;
    inlmeas.isms=isms;
    inlmeas.ismt=ismt;
    inlmeas.smpars=ismps;
    inlmeas.smpart=ismpt;

    if(DIM==2)
        sprintf(base,"%dx%d_SU%d_b%.3f_id%d",LENGT,LENGS1,
                SUN,beta,id);
    else if(DIM==3)
        sprintf(base,"%dx%dx%d_SU%d_b%.3f_id%d",LENGT,LENGS1,LENGS2,
                SUN,beta,id);
    else
        sprintf(base,"%dx%dx%dx%d_SU%d_b%.1f_id%d",LENGT,LENGS1,LENGS2,LENGS3,
                SUN,beta,id);
    
    
    sprintf(LOG_FILE,"%s_topQ.test",base);
    sprintf(CNFG_FILE,"%s/%s_n%d",dir,base,sconf);
    sprintf(config_dir, "%s", dir);
    
    // Check FILES.
    FILE *fset=NULL;
    char check[NAME_SIZE];
      
    
    ftest=fopen(LOG_FILE,"w");
    error(ftest==NULL,"read_inp_meas [meas.c]",
          "Unable to create .meas file!");
    fclose(ftest);

}



 
int main(int argc,char *argv[])
{
    int n,ir,it,nr,nt,il;
    int seed;
    char out_dir[NAME_SIZE];
    
    double t1,t2;
    double *qtop;
    
    
    read_inp_meas(&seed,out_dir,argc,argv);
    rlxd_init(1,seed);
    neib_init();
    init_gauge(0);

    
    
    
    qtop= malloc(sizeof(double));
    
    for (n=sconf; n <=fconf; n+=500) {
        
        t1=get_time();
        printf("Configuration : %04d\n", n);
        sprintf(CNFG_FILE, "%s/conf.%d", config_dir, n);

        
        for(inlmeas.isms=0; inlmeas.isms<=200; inlmeas.isms+=20){

            
            read_config_fromSU2code(pu, CNFG_FILE);
        // read_config(CNFG_FILE);
            inlmeas.ismt=inlmeas.isms;
            
           
        // Smearing.
        //smearing_APE_spatial(inlmeas.isms,inlmeas.smpars,pu);
        //smearing_APE_temporal(inlmeas.ismt, inlmeas.smpart, pu);
            smearing_APE_all(inlmeas.isms, inlmeas.smpars, pu);
        
        // Measurement.
        meas_topologicalcharge(qtop);
        
        // Printing.
        logging("Q %d %d %e\n",inlmeas.isms, n, *qtop);
            printf("Q %d %d %e\n",inlmeas.isms, n, *qtop);
        }
        
        
    }

    finish_gauge();
    


    return 0;
}
