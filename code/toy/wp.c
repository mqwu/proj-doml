#include "comm.h"
#include "nrutil.h"
#include "nrutil.c"
#include "lib.c"
#include "cost.c"

main()
{
int stime;
long ltime;
int i, j, k, L;
double *x, *zin, *zout, *ave, energy;
FILE *ins, *ins0, *outs;
char buf[255];


ltime=time(NULL);
stime=(unsigned int)ltime/2;
srand(stime);

if(shortcut==0) dim=(P+1)*hd1+(1+hd1)*OUT_UNIT;
  else dim=(P+1)*hd1+(1+hd1+P)*OUT_UNIT;
x=dvector(1,dim);
zin=dvector(1,P);
zout=dvector(1,OUT_UNIT);
ave=dvector(1,40000);

for(i=1; i<=40000; i++) ave[i]=0.0;

ins0=fopen("neural.ww","r");
fscanf(ins0, " %lf",&energy);
L=0;
while(!feof(ins0)){

  for(i=1; i<=dim; i++) fscanf(ins0, " %lf", &x[i]);

  ins=fopen("plot.dat","r");
  fscanf(ins, " %lf", &zin[1]); 
  k=0;
  while(!feof(ins)){
     for(i=2; i<=P; i++) fscanf(ins, " %lf", &zin[i]);
     fden(zin, x, zout);

     k++;
     /*
     if(zout[1]>0.5) ave[k]+=1;
     */
     ave[k]+=zout[1];

     fscanf(ins, " %lf", &zin[1]);
  }
  fclose(ins);
  
  fscanf(ins0, " %lf", &energy);
  L++;

printf("L=%d\n", L);
 } 
  fclose(ins0);

for(i=1; i<=k; i++) ave[i]=ave[i]/L;

ins=fopen("plot.dat","r");
outs=fopen("plot.out","w");
fscanf(ins, " %lf", &zin[1]);
k=0;
while(!feof(ins)){
    for(i=2; i<=P; i++) fscanf(ins, " %lf", &zin[i]);
    k++;
    for(i=1; i<=P; i++) fprintf(outs, " %g", zin[i]);
    if(ave[k]>0.5) fprintf(outs,  " %g   %d\n",ave[k], 1);
      else fprintf(outs, "  %g   %d\n", ave[k],0);
   
   fscanf(ins, " %lf", &zin[1]);
  }
fclose(outs);
fclose(ins);

free_dvector(x,1,dim);
free_dvector(zin,1,P);
free_dvector(zout,1,OUT_UNIT);
free_dvector(ave,1,40000);

return 0;
}
