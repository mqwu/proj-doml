#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
 
#define pi 3.14159265
 
#define maxrow 100000
#define maxcol 10000

#include "nrutil.h"
#include "nrutil.c"

int main(int argc, char *argv[])
{
int stime;
long ltime;
double **data_mat, **data_y;
int data_num, DATA_NUM, SELECTION, P, Q, Ycolumn;
int i,j,l,k,k0,k1,k2,m,it,iter,grid, ok;
double **dist, aaa, max;
int *sel, *rem, **aux, S, R, choice;
FILE *ins;

ltime=time(NULL);
stime=(unsigned int)ltime/2;
srand(stime);

ins=fopen(argv[1],"r");
if(ins==NULL){ printf("can't open the data file\n"); return 1; }
fscanf(ins, " %d %d %d", &DATA_NUM, &data_num, &SELECTION);
fscanf(ins, " %d %d %d", &Q, &P, &Ycolumn);  
fclose(ins);  

// printf(" %d %d %d\n", Q, P, Ycolumn);

if(DATA_NUM+data_num>maxrow){ printf("The sample size is too large\n"); return 1; }
if(data_num>SELECTION){ printf("The sample size at the current time is too large\n"); return 1; }

sel=ivector(1,DATA_NUM+data_num);
rem=ivector(1,DATA_NUM+data_num);
aux=imatrix(1,DATA_NUM+data_num,1,Q);
dist=dmatrix(1,DATA_NUM+data_num,1,DATA_NUM+data_num);
data_mat=dmatrix(1,DATA_NUM+data_num,1,P);
data_y=dmatrix(1,DATA_NUM+data_num,1,Ycolumn);


ins=fopen(argv[2],"r");
if(ins==NULL){ printf("can't open the data file\n"); return 1; }
for(i=1; i<=DATA_NUM; i++){
    for(j=1; j<=Q; j++) fscanf(ins, " %d", &aux[i][j]);// printf(" %d", aux[i][j]); }
    for(j=1; j<=P; j++) fscanf(ins, " %lf", &data_mat[i][j]); 
    for(j=1; j<=Ycolumn; j++) fscanf(ins, " %lf", &data_y[i][j]);
    // printf("\n");
   }
fclose(ins); 

ins=fopen(argv[3],"r");
if(ins==NULL){ printf("can't open the data file\n"); return 1; }
for(i=DATA_NUM+1; i<=DATA_NUM+data_num; i++){
    for(j=1; j<=Q; j++) fscanf(ins, " %d", &aux[i][j]); // printf(" %d", aux[i][j]); } 
    for(j=1; j<=P; j++) fscanf(ins, " %lf", &data_mat[i][j]); 
    for(j=1; j<=Ycolumn; j++) fscanf(ins, " %lf", &data_y[i][j]);
    // printf("\n"); 
   }
fclose(ins);



if(DATA_NUM+data_num<SELECTION){  

    ins=fopen(argv[4],"w");
    if(ins==NULL){ printf("can't open the data file\n"); return 1; }

    for(i=1; i<=DATA_NUM+data_num; i++){
        for(j=1; j<=Q; j++) fprintf(ins, " %d", aux[i][j]);
        for(j=1; j<=P; j++) fprintf(ins, " %g", data_mat[i][j]);
        for(j=1; j<=Ycolumn; j++) fprintf(ins, " %g", data_y[i][j]);
        fprintf(ins, "\n");
      }
    fclose(ins);
  }
else{

  for(i=1; i<=DATA_NUM+data_num; i++){
    dist[i][i]=0.0;
    for(j=i+1; j<=DATA_NUM+data_num; j++){
      dist[i][j]=0.0;
      for(k=1; k<=P; k++) 
          dist[i][j]+=(data_mat[i][k]-data_mat[j][k])*(data_mat[i][k]-data_mat[j][k]);
      dist[i][j]=sqrt(dist[i][j]);
      dist[j][i]=dist[i][j];
     }
   }

  S=data_num;  
  for(i=1; i<=data_num; i++) sel[i]=DATA_NUM+i;
  R=DATA_NUM;
  for(i=1; i<=R; i++) rem[i]=i;

  for(k=S+1; k<=SELECTION; k++){ 
 
      max=0.0;
      for(i=1; i<=R; i++){
          aaa=0.0;
          for(j=1; j<=S; j++)
             if(aaa<dist[rem[i]][sel[j]]) aaa=dist[rem[i]][sel[j]];
          if(aaa>max){ max=aaa; choice=i; }
        }                 
 
      S++; R--;
      sel[S]=rem[choice]; 

      for(i=choice; i<=R-1; i++) rem[i]=rem[i+1]; 
     } 

   // for(k=1; k<=SELECTION; k++) printf(" %d", sel[k]);
   // printf("\n"); 
     
   ins=fopen(argv[4],"w");
   if(ins==NULL){ printf("can't open the data file\n"); return 1; }
   for(k=1; k<=SELECTION; k++){
       i=sel[k]; // printf(" %d\n", i);
       for(j=1; j<=Q; j++) fprintf(ins, " %d", aux[i][j]);
       for(j=1; j<=P; j++) fprintf(ins, " %g", data_mat[i][j]);
       for(j=1; j<=Ycolumn; j++) fprintf(ins, " %g", data_y[i][j]); 
       fprintf(ins, "\n");
      }
   fclose(ins);
  }
  

free_ivector(sel,1,DATA_NUM+data_num);
free_ivector(rem,1,DATA_NUM+data_num);
free_imatrix(aux,1,DATA_NUM+data_num,1,Q);
free_dmatrix(dist,1,DATA_NUM+data_num,1,DATA_NUM+data_num);
free_dmatrix(data_mat,1,DATA_NUM+data_num,1,P);
free_dmatrix(data_y,1,DATA_NUM+data_num,1,Ycolumn);

return 0;
}
