#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
 
#define pi 3.14159265
 
#define data_num  2250
 
#define P 3            /* the number of input units */
#define hd1 30         /* 15 the number of first hidden layer units */
#define OUT_UNIT 1     /* the number of output units */
static int shortcut=0;

#define total_iteration 5000000
#define stepscale  10000
#define WARM           1
#define scale          1 
#define N              1         /* population size */

static double lowE=0.0, maxEE=6000.0, maxE, range=0.0;
static double tau=0.6, rho=1.0, delta,stepsize, hightem=10.0; 
static double **solution,*min,*refden,**data_mat, **data_y, VARIANCE;
static double accept_loc, total_loc,accept_dir, total_dir;
static int dim, sze, Best=5;

#include "nrutil.h"
#include "nrutil.c"
#include "lib.c"
#include "cost.c"
#include "metmove2.c"
#include "metrop.c"
#include "metrop2.c"

main()
{
int stime;
long ltime;
int i,j,k,k0,k1,k2,it,iter,grid;
double **x, **hist, *fvalue, sum, max, var;
double fret, ave, tem;
int cg=1, Repeat,batchsize;
FILE *ins;

ltime=time(NULL);
stime=(unsigned int)ltime/2;
stime=745708026;
srand(stime);

ins=fopen("ee.log","a");
fprintf(ins, "random seed=%d  N=%d tau=%g rho=%g  scale=%d hd1=%d data_num=%d\n", 
stime,N,tau,rho,scale,hd1,data_num);
fclose(ins);

ins=fopen("ee.sum","a");
fprintf(ins, "random seed=%d  N=%d tau=%g rho=%g  scale=%d hd1=%d data_num=%d\n",
stime,N,tau,rho,scale,hd1,data_num);
fclose(ins);


grid=ceil((maxEE+2*range-lowE)*scale);
if(shortcut==0) dim=(P+1)*hd1+(1+hd1)*OUT_UNIT;
  else dim=(P+1)*hd1+(1+hd1+P)*OUT_UNIT;

x=dmatrix(1,N,1,dim);  /* store the current samples */
fvalue=dvector(1,N);
hist=dmatrix(0,grid,1,3);
refden=dvector(0,grid);
solution=dmatrix(0,Best,0,dim);
min=dvector(0, Best);
data_mat=dmatrix(1,data_num,1,P);
data_y=dmatrix(1,data_num,1,OUT_UNIT);


ins=fopen("toy2.dat","r");
if(ins==NULL){ printf("can't open the data file\n"); return 1; }
for(i=1; i<=data_num; i++){
    for(j=1; j<=P; j++) fscanf(ins, " %lf", &data_mat[i][j]);
    for(j=1; j<=OUT_UNIT; j++) fscanf(ins, " %lf", &data_y[i][j]);
   }
fclose(ins);


for(Repeat=1; Repeat<=1; Repeat++){
                                                                                                                            
   /* Initialization */
   for(sum=0.0, i=0; i<=grid; i++){ refden[i]=exp(-0.001*i); sum+=refden[i]; }
   for(i=0; i<=grid; i++){ refden[i]/=sum; printf("i=%d %g\n",i,refden[i]); }
   for(i=0; i<=grid; i++){
       hist[i][1]=lowE+i*1.0/scale;
       hist[i][2]=0.0;
       hist[i][3]=0.0;
     }

   for(i=1; i<=N; i++){
     for(j=1; j<=dim; j++) x[i][j]=gasdev()*0.01;
     fvalue[i]=cost(x[i]);
     printf("favlue[%d]=%g\n", i, fvalue[i]);
    /*
     while(fvalue[i]>=maxEE+range){
         for(j=1; j<=1000; j++) metropolis(x[i],&fvalue[i]);
         printf("fvalue[%d]=%g\n", i,fvalue[i]);
        }
     */
    }
  /*
   ins=fopen("ee.sol","r");
   for(i=1; i<=N; i++){
       fscanf(ins, " %lf", &fvalue[i]);
       for(j=1; j<=dim; j++) fscanf(ins, " %lf", &x[i][j]);
       fvalue[i]=cost(x[i]);
      }
   fclose(ins);
  */

   for(i=1; i<=Best; i++){
      min[i]=fvalue[1];
      for(j=1; j<=dim; j++) solution[i][j]=x[1][j];
    }
   maxE=max_vector(fvalue,N);


   accept_loc=accept_dir=total_loc=total_dir=0.1;
   for(iter=1; iter<=total_iteration; iter++){

       if(iter<=WARM*stepscale) delta=rho;
           else delta=rho*exp(-tau*log(1.0*(iter-(WARM-1)*stepscale)/stepscale));

       if(iter<=100) tem=hightem+0.2;
           else tem=hightem/sqrt(1.0*iter/100.0)+0.2;

       /*
       if(delta<0.1) stepsize=2.0;
         else if(delta<0.05) stepsize=1.0;
            else  stepsize=0.5;
       */

       if(iter<500000) stepsize=0.5;
          else if(iter<1000000) stepsize=1.0;
             else stepsize=2.0;
       
       MetmoveSA(x,fvalue,hist,tem);
       
       
       if(iter%10000==0){
         printf("tem=%g delta=%g iter=%d stepsize=%g min=%g max=%g\n",tem,delta,iter,stepsize,min[Best],maxE);
         for(i=1; i<=N; i++) printf(" %g",fvalue[i]);
         printf("\n");
        }
     
       
       if(min[Best]<0.21) goto ABC;
     }
 
ABC:

    for(sum=0.0,k0=0,i=0; i<=sze; i++)
       if(hist[i][3]<=0.0){ sum+=refden[i]; k0++; }
    if(k0>0) ave=sum/k0;
       else ave=0.0;
    for(i=0; i<=sze; i++) hist[i][2]=hist[i][2]+log(refden[i]+ave);
    max=hist[0][2];
    for(i=1; i<=sze; i++)
       if(hist[i][2]>max) max=hist[i][2];
    for(sum=0.0, i=0; i<=sze; i++){ hist[i][2]-=max; sum+=exp(hist[i][2]); }
    for(i=0; i<=sze; i++) hist[i][2]=hist[i][2]-log(sum)+log(100.0);

   ins=fopen("ee.log","a");
   fprintf(ins, "minimum energy=%12.6f\n",min[Best]);
   for(j=1; j<=dim; j++){
       fprintf(ins, " %10.6f",solution[Best][j]);
       if(j%8==0) fprintf(ins,"\n");
      }
   fprintf(ins,"\n");
   fprintf(ins, "rates=%g  %g  min=%g\n",accept_loc/total_loc,accept_dir/total_dir,min[Best]);
   fclose(ins);
                                                                                                        
   ins=fopen("ee.est", "w");
   fprintf(ins, "delta=%g \n", delta);
   if(ins==NULL){ printf("Can't write to file\n"); return 1; }
   for(i=0; i<=sze; i++){
       fprintf(ins, "%5d  %10.6f  %10.6f  %10.6f  %g\n",i,hist[i][1],exp(hist[i][2]),
               hist[i][3],hist[i][2]);
       hist[i][3]=0.0;
     }
   fclose(ins);
                                                                                                        

  ins=fopen("ee.sol","a");
  for(i=1; i<=Best; i++){
      fprintf(ins, " %g\n", min[i]);
      for(j=1; j<=dim; j++){
          fprintf(ins, " %10.6f", solution[i][j]);
          if(j%8==0) fprintf(ins,"\n");
        }
      fprintf(ins, "\n\n");
     }
  fclose(ins);


  fret=min[Best];
  /*
  metropolis_post(solution[Best],&fret);
  */
  ins=fopen("ee.sum","a");
  fprintf(ins, "%d  %10d  %10.4f\n",Repeat,iter-1,fret);
  fclose(ins);
} /* end repeat */

free_dmatrix(x,1,N,1,dim);  
free_dvector(fvalue,1,N);
free_dmatrix(hist,0,grid,1,3);
free_dvector(refden,0,grid);
free_dmatrix(solution,0,Best,0,dim);
free_dvector(min, 0, Best);
free_dmatrix(data_mat,1,data_num,1,P);
free_dmatrix(data_y,1,data_num,1,OUT_UNIT);


return 0;
}