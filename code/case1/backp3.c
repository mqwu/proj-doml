#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
 
#define pi 3.14159265
 
static int shortcut=0;
static double **data_mat, **data_y;
static int data_num, Layer, P, OUT_UNIT, HD[100], dim;

#include "nrutil.h"
#include "nrutil.c"
#include "lib2.c"
#include "derivative3.c"

int main(int argc, char *argv[])
{
int stime;
long ltime;
int NV, ND, NES, v, i, j, k, m, k0,k1, k2, s1,s2,activation, generation, ID;
char frname[50], fwname[50];
double *x, *derv;
FILE *ins;


ins=fopen(argv[1], "r");
if(ins==NULL){  printf("No model information files exist\n");  return 0; }
fscanf(ins, " %d %d %d", &NV, &ND, &NES);
fscanf(ins, " %d", &generation);
fscanf(ins, " %d", &data_num);
fscanf(ins, " %d", &v);  // position of the response variable in original dataset
fscanf(ins, " %d", &activation); 
fscanf(ins, " %d", &Layer); 
if(Layer>99){ printf("The neural network has too many layers\n"); return 0; } 

for(i=0; i<Layer; i++) fscanf(ins, " %d", &HD[i]); 
fclose(ins); 

if(NV<10) k0=10;
  else if(NV<100) k0=100;
    else if(NV<1000) k0=1000;
      else if(NV<10000) k0=10000;
        else{ k0=100000; printf("NV is possibly too large\n"); }
if(ND<10) k1=10;
  else if(ND<100) k1=100;
    else if(ND<1000) k1=1000;
      else if(ND<10000) k1=10000;
        else{ k1=100000; printf("ND is possibly too large\n"); }
if(NES<10) k2=10;
  else if(NES<100) k2=100;
    else if(NES<1000) k2=1000;
      else if(NES<10000) k2=10000;
        else{ k2=100000; printf("NES is possibly too large\n"); }


P=HD[0];  OUT_UNIT=HD[Layer-1]; 
for(dim=0,j=1; j<Layer; j++) dim+=(1+HD[j-1])*HD[j];
if(shortcut==1) dim+=HD[0]*HD[Layer-1];


x=dvector(1,dim);  /* store the current samples */
derv=dvector(1,P);
data_mat=dmatrix(1,50000,1,P);
data_y=dmatrix(1,50000,1,OUT_UNIT); 


 if(generation<10){
     frname[0]='X'; frname[1]='B'; frname[2]='A'; frname[3]='S'; frname[4]='E'; frname[5]='C'; 
     frname[6]='Y'; frname[7]='0'; frname[8]='0'; frname[9]='0';
     frname[10]=(char)(generation+48);
     frname[11]='.'; frname[12]='d'; frname[13]='a'; frname[14]='t';  frname[15]='\0';  }
    else if(generation<100){
     frname[0]='X'; frname[1]='B'; frname[2]='A'; frname[3]='S'; frname[4]='E'; frname[5]='C';
     frname[6]='Y'; frname[7]='0'; frname[8]='0'; frname[9]=(char)(generation/10+48);
     frname[10]=(char)(generation%10+48);
     frname[11]='.'; frname[12]='d'; frname[13]='a'; frname[14]='t';  frname[15]='\0'; }
    else if(generation<1000){
      frname[0]='X'; frname[1]='B'; frname[2]='A'; frname[3]='S'; frname[4]='E'; frname[5]='C';
      frname[6]='Y'; frname[7]='0';
      frname[8]=(char)(generation/100+48); frname[9]=(char)((generation%100)/10+48);
      frname[10]=(char)(generation%10+48);
      frname[11]='.'; frname[12]='d'; frname[13]='a'; frname[14]='t';  frname[15]='\0'; } 
     else if(generation<10000){ 
      frname[0]='X'; frname[1]='B'; frname[2]='A'; frname[3]='S'; frname[4]='E'; frname[5]='C';
      frname[6]='Y'; frname[7]=(char)(generation/1000+48); frname[8]=(char)((generation%1000)/100+48);
      frname[9]=(char)((generation%100)/10+48); frname[10]=(char)(generation%10+48);
      frname[11]='.'; frname[12]='a'; frname[13]='a'; frname[14]='t';  frname[15]='\0'; }
     else{ printf("generation is too large\n"); return 1; }

 
ins=fopen(frname,"r");
if(ins==NULL){ printf("can't open the data file\n"); return 1; }

fscanf(ins, " %lf", &data_mat[1][1]);
data_num=0;
while(!feof(ins)){
    // for(j=1; j<=OUT_UNIT; j++) fscanf(ins, " %lf", &data_y[i][j]);
    for(j=2; j<=P; j++) fscanf(ins, " %lf", &data_mat[data_num+1][j]);
    data_num++;
    fscanf(ins, " %lf", &data_mat[data_num+1][1]);
   }
fclose(ins); 

if(data_num!=NV*NES){ printf("Data format might be wrong: # of points is not equal to NV*ND\n"); }

for(m=0, k=1; k<=Layer-1; k++){
    
   if(k<10){
     fwname[0]='w'; fwname[1]='b';
     fwname[2]=(char)(k+48);
     fwname[3]='.'; fwname[4]='t'; fwname[5]='x'; fwname[6]='t';  fwname[8]='\0';  }
    else if(k<100){
      fwname[0]='w'; fwname[1]='b';
      fwname[2]=(char)(k/10+48); fwname[3]=(char)(k%10+48);
      fwname[4]='.'; fwname[5]='t'; fwname[6]='x'; fwname[7]='t';  fwname[8]='\0'; }
     else{ printf("Too many layers in DNN\n"); return 1; }
    
   // printf(" %s\n", fwname);

   ins=fopen(fwname, "r");
   if(ins==NULL){ printf("The weight-bias file does not exist\n"); return 1; }

    for(i=1; i<=(HD[k-1]+1)*HD[k]; i++){
        m++;
        fscanf(ins, " %lf", &x[m]); 
       }
    fclose(ins); 

    // printf("dim=%d m=%d\n", dim,m); 
  }


 if(v<10){
     frname[0]='X'; frname[1]='Y'; frname[2]=(char)(v+48);
     frname[3]='.'; frname[4]='d'; frname[5]='e'; frname[6]='r';  frname[7]='\0';  }
    else if(v<100){
     frname[0]='X'; frname[1]='Y'; frname[2]=(char)(v/10+48);     frname[3]=(char)(v%10+48);
     frname[4]='.'; frname[5]='d'; frname[6]='e'; frname[7]='r';  frname[8]='\0'; }
    else if(v<1000){
      frname[0]='X'; frname[1]='Y'; frname[2]=(char)(v/100+48); frname[3]=(char)((v%100)/10+48);
      frname[4]=(char)(v%10+48);
      frname[5]='.'; frname[6]='d'; frname[7]='e'; frname[8]='r';  frname[9]='\0'; }
     else if(v<10000){
      frname[0]='X'; frname[1]='Y'; frname[2]=(char)(v/1000+48); frname[3]=(char)((v%1000)/100+48);
      frname[4]=(char)((v%100)/10+48); frname[5]=(char)(v%10+48);
      frname[6]='.'; frname[7]='d'; frname[8]='e'; frname[9]='r';  frname[10]='\0'; }
     else{ printf("v is too large\n"); return 1; }


ins=fopen(frname, "w");
for(i=1; i<=data_num; i++){ 

  s1=i%NV; // order within NV
  s2=i/NV; // order within NES
  ID=s2*k0*k1+v*k1+s1; 
   
  if(activation==0) input_derivative(data_mat[i],x,derv);  // logistic
   else if(activation==1) input_derivative_tanh(data_mat[i],x,derv);  // tanh 
      else if(activation==2) input_derivative_relu(data_mat[i],x,derv);  // RELU 

  fprintf(ins, " %d", ID); 
  for(k=1; k<=P; k++) fprintf(ins, " %10.6f", derv[k]);
  fprintf(ins, "\n");
 }
fclose(ins); 
      

free_dvector(x,1,dim);  
free_dvector(derv,1,P);
free_dmatrix(data_mat,1,50000,1,P);
free_dmatrix(data_y,1,50000,1,OUT_UNIT);


return 0;
}
