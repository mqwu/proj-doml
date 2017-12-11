#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
 
#define pi 3.14159265
 
#include "nrutil.h"
#include "nrutil.c"
#include "lib2.c"

int main(int argc, char *argv[])
{
int stime;
long ltime;
int i,j,l,k,v, xdim, ydim, datanum, iternum, miniter, maxiter;
int NI, NR, ND, NES;
double *x, *y;
int mode, sample, iter, succ, Okay;
FILE *ins, *ins1, *ins2, *outs;
char frname[15];

ins=fopen(argv[1], "r");
fscanf(ins, " %d %d %d %d", &NI, &NR, &ND, &NES);
fclose(ins);

xdim=NI+NR;
ydim=ND; 
datanum=NES;

x=dvector(1,xdim);
y=dvector(1,ydim);
 
miniter=maxiter=1;

ins1=fopen(argv[2],"r");
ins2=fopen(argv[3],"r");
 
fscanf(ins2, " %d", &mode);
while(!feof(ins2)){

   fscanf(ins2, " %d %d %d", &sample, &iter, &succ);

   if(iter<miniter) miniter=iter;
      else if(iter>maxiter) maxiter=iter; 
   
   if(iter<10){
     frname[0]='B'; frname[1]='i'; frname[2]='m'; frname[3]='G'; frname[4]=(char)(iter+48);
     frname[5]='.'; frname[6]='d'; frname[7]='a'; frname[8]='t';  frname[9]='\0';  }
    else if(iter<100){
      frname[0]='B'; frname[1]='i'; frname[2]='m';  frname[3]='G';
      frname[4]=(char)(iter/10+48); frname[5]=(char)(iter%10+48);
      frname[6]='.'; frname[7]='d'; frname[8]='a'; frname[9]='t'; frname[10]='\0'; }
     else{ printf("iter is too large\n"); return 1; }

     outs=fopen(frname, "w");
     for(k=1; k<=datanum; k++){

        for(i=1; i<=ydim; i++) fscanf(ins2, " %lf", &y[i]);
        Okay=1; i=0;   
        while(Okay==1 && i<ydim){
              i++;
              if(y[i]< -999.0) Okay=0;
             }
  
        fscanf(ins1, " %d %d %d", &mode, &sample, &iter); 
        for(i=1; i<=xdim; i++) fscanf(ins1, " %lf", &x[i]); 

        if(succ==0 || Okay==0){ printf("Failed simulation: mode=%d sample=%d iter=%d\n", mode, sample, iter); }
           else{
             
             fprintf(outs, " %3d %3d %3d", mode, sample, iter);
             for(i=1; i<=xdim; i++) fprintf(outs, " %g", x[i]);
             for(i=1; i<=ydim; i++) fprintf(outs, " %g", y[i]);
             fprintf(outs, "\n");
            }

         fscanf(ins2, " %d %d %d %d", &mode, &sample, &iter, &succ);
        }
      fclose(outs);

      fscanf(ins2, " %d", &mode); 
    } 
       
  fclose(ins1); 
  fclose(ins2);

  ins=fopen("Dim.par","w");
  fprintf(ins, " %d %d %d %d %d\n", xdim, ydim, miniter, maxiter,NES);
  fclose(ins); 

free_dvector(x,1,xdim);  
free_dvector(y,1,ydim);

return 0;
}
