#define LOGISTIC(z) 1.0/(1.0+exp(-z))
#define TANH(z) (1.0-exp(-2*z))/(1.0+exp(-2*z))
#define RELU(z) (z>0)?z:0 

int derivative_tanh();
double cost_derivative_tanh(z,d)
double *z, *d;
{
double *pout,*dobs, sum,sum1,sum2,sum3,sum4,max,min;
int i, j,k, m1, m2, m3;
FILE *ins;

dobs=dvector(1,dim);

max=max_vector(z,dim);
min=min_vector(z,dim);
if(max>50.0 || min<-50.0) return 1.0e+100;

pout=dvector(1,OUT_UNIT);

for(j=1; j<=dim; j++) d[j]=0.0; 
 
  for(sum=0.0, i=1; i<=data_num; i++){
     derivative_tanh(data_mat[i],data_y[i],z,pout,dobs);
     for(j=1; j<=dim; j++) d[j]+=dobs[j]; 
     
     for(j=1; j<=OUT_UNIT; j++)
         sum+=0.5*(data_y[i][j]-pout[j])*(data_y[i][j]-pout[j]);
    }
printf("cost=%g\n", sum); 

free_dvector(pout,1,OUT_UNIT);
free_dvector(dobs,1,dim);

return sum;
}



int derivative_tanh(ox,oy,z,pout,dobs)
double *ox,*oy,*z, *pout,*dobs; 
{
int i, j, k, L, m, M;
double **ps, ave;
double **delta;
 
M=HD[0];
for(j=1; j<Layer; j++){
    if(M<HD[j]) M=HD[j];
   }
ps=dmatrix(0,Layer,1,M); 
delta=dmatrix(0,Layer,1,M); 

for(j=1; j<=P; j++) ps[0][j]=ox[j];


/* calculate the outputs of hidden layers */
m=0;
for(L=1; L<Layer-1; L++){
    
  for(i=1; i<=HD[L]; i++){
      k=m+(HD[L-1]+1)*(i-1);
      ps[L][i]=z[k+1];
      for(j=k+2; j<=k+HD[L-1]+1; j++) ps[L][i]+=ps[L-1][j-k-1]*z[j];
      ps[L][i]=TANH(ps[L][i]);
    }
  m+=(HD[L-1]+1)*HD[L];
 }


/* calculate the predicted mean from the mean unit */
if(shortcut==0){
   for(i=1; i<=OUT_UNIT; i++){
      k=m+(1+HD[Layer-2])*(i-1);
      ave=z[k+1];
      for(j=k+2; j<=k+HD[Layer-2]+1; j++) ave+=ps[Layer-2][j-k-1]*z[j];
      pout[i]=ave;
     }
   m+=(1+HD[Layer-2])*OUT_UNIT;
  }
 else{
   for(i=1; i<=OUT_UNIT; i++){
     k=m+(1+HD[Layer-2]+HD[0])*(i-1);
     ave=z[k+1];
     for(j=k+2; j<=k+HD[Layer-2]+1; j++) ave+=ps[Layer-2][j-k-1]*z[j];
     for(j=k+HD[Layer-2]+2; j<=k+HD[0]+HD[Layer-2]+1; j++) ave+=ps[0][j-k-HD[Layer-2]-1]*z[j];
     pout[i]=ave;
    }
   m+=(1+HD[0]+HD[Layer-2])*OUT_UNIT; 
  }
 
  

/****** calculation of delta ********/ 
 
 for(i=1; i<=OUT_UNIT; i++) delta[Layer-1][i]=(pout[i]-oy[i]);
 for(L=Layer-2; L>=1; L--){

   if(L==Layer-2 && shortcut==1) m-=(1+HD[0]+HD[L])*HD[L+1];
      else m-=(1+HD[L])*HD[L+1]; 
     
   for(i=1; i<=HD[L]; i++){
       delta[L][i]=0.0; 
       for(j=1; j<=HD[L+1]; j++) delta[L][i]+=z[m+(j-1)*(HD[L]+1)+i+1]*delta[L+1][j];
       delta[L][i]*=(1+ps[L][i])*(1-ps[L][i]); 
     }
   }

/****** calculation of gradient  **************/
 m=0;
 for(L=1; L<=Layer-2; L++){
   for(i=1; i<=HD[L]; i++){
      k=m+(HD[L-1]+1)*(i-1);
      dobs[k+1]=delta[L][i];
      for(j=k+2; j<=k+HD[L-1]+1; j++) dobs[j]=delta[L][i]*ps[L-1][j-k-1];
     }
   m+=(HD[L-1]+1)*HD[L];
  }

 if(shortcut==0){
    for(i=1; i<=OUT_UNIT; i++){
        k=m+(1+HD[Layer-2])*(i-1);

        dobs[k+1]=delta[Layer-1][i];
        for(j=k+2; j<=k+HD[Layer-2]+1; j++) dobs[j]=delta[Layer-1][i]*ps[Layer-2][j-k-1];
       }
     }
  else{
   for(i=1; i<=OUT_UNIT; i++){
     k=m+(1+HD[Layer-2]+HD[0])*(i-1);

     dobs[k+1]=delta[Layer-1][i];
     for(j=k+2; j<=k+HD[Layer-2]+1; j++) dobs[j]=delta[Layer-1][i]*ps[Layer-2][j-k-1];
     for(j=k+HD[Layer-2]+2; j<=k+HD[0]+HD[Layer-2]+1; j++) dobs[j]=delta[Layer-1][i]*ps[0][j-k-HD[Layer-2]-1];
    }
  }

 free_dmatrix(ps,0,Layer,1,M);
 free_dmatrix(delta,0,Layer,1,M);

return 0;
}

 




int derivative_relu();
double cost_derivative_relu(z,d)
double *z, *d;
{
double *pout,*dobs, sum,sum1,sum2,sum3,sum4,max,min;
int i, j,k, m1, m2, m3;
FILE *ins;

dobs=dvector(1,dim);

max=max_vector(z,dim);
min=min_vector(z,dim);
if(max>50.0 || min<-50.0) return 1.0e+100;

pout=dvector(1,OUT_UNIT);

for(j=1; j<=dim; j++) d[j]=0.0; 
 
  for(sum=0.0, i=1; i<=data_num; i++){
     derivative_relu(data_mat[i],data_y[i],z,pout,dobs);
     for(j=1; j<=dim; j++) d[j]+=dobs[j]; 
     
     for(j=1; j<=OUT_UNIT; j++)
         sum+=0.5*(data_y[i][j]-pout[j])*(data_y[i][j]-pout[j]);
    } 
printf("cost=%g\n", sum); 

free_dvector(pout,1,OUT_UNIT);
free_dvector(dobs,1,dim);

return sum;
}



int derivative_relu(ox,oy,z,pout,dobs)
double *ox,*oy,*z, *pout,*dobs; 
{
int i, j, k, L, m, M;
double **ps, ave;
double **delta;
 
M=HD[0];
for(j=1; j<Layer; j++){
    if(M<HD[j]) M=HD[j];
   }
ps=dmatrix(0,Layer,1,M); 
delta=dmatrix(0,Layer,1,M); 

for(j=1; j<=P; j++) ps[0][j]=ox[j];


/* calculate the outputs of hidden layers */
m=0;
for(L=1; L<Layer-1; L++){
    
  for(i=1; i<=HD[L]; i++){
      k=m+(HD[L-1]+1)*(i-1);
      ps[L][i]=z[k+1];
      for(j=k+2; j<=k+HD[L-1]+1; j++) ps[L][i]+=ps[L-1][j-k-1]*z[j];
      ps[L][i]=RELU(ps[L][i]);
    }
  m+=(HD[L-1]+1)*HD[L];
 }


/* calculate the predicted mean from the mean unit */
if(shortcut==0){
   for(i=1; i<=OUT_UNIT; i++){
      k=m+(1+HD[Layer-2])*(i-1);
      ave=z[k+1];
      for(j=k+2; j<=k+HD[Layer-2]+1; j++) ave+=ps[Layer-2][j-k-1]*z[j];
      pout[i]=ave;
     }
   m+=(1+HD[Layer-2])*OUT_UNIT;
  }
 else{
   for(i=1; i<=OUT_UNIT; i++){
     k=m+(1+HD[Layer-2]+HD[0])*(i-1);
     ave=z[k+1];
     for(j=k+2; j<=k+HD[Layer-2]+1; j++) ave+=ps[Layer-2][j-k-1]*z[j];
     for(j=k+HD[Layer-2]+2; j<=k+HD[0]+HD[Layer-2]+1; j++) ave+=ps[0][j-k-HD[Layer-2]-1]*z[j];
     pout[i]=ave;
    }
   m+=(1+HD[0]+HD[Layer-2])*OUT_UNIT; 
  }
 
  

/****** calculation of delta ********/ 
 
 for(i=1; i<=OUT_UNIT; i++) delta[Layer-1][i]=(pout[i]-oy[i]);
 for(L=Layer-2; L>=1; L--){

   if(L==Layer-2 && shortcut==1) m-=(1+HD[0]+HD[L])*HD[L+1];
      else m-=(1+HD[L])*HD[L+1]; 
     
   for(i=1; i<=HD[L]; i++){
       delta[L][i]=0.0; 
       for(j=1; j<=HD[L+1]; j++) delta[L][i]+=z[m+(j-1)*(HD[L]+1)+i+1]*delta[L+1][j];
       delta[L][i]*=(ps[L][i]>0)?1:0; 
     }
   }

/****** calculation of gradient  **************/
 m=0;
 for(L=1; L<=Layer-2; L++){
   for(i=1; i<=HD[L]; i++){
      k=m+(HD[L-1]+1)*(i-1);
      dobs[k+1]=delta[L][i];
      for(j=k+2; j<=k+HD[L-1]+1; j++) dobs[j]=delta[L][i]*ps[L-1][j-k-1];
     }
   m+=(HD[L-1]+1)*HD[L];
  }

 if(shortcut==0){
    for(i=1; i<=OUT_UNIT; i++){
        k=m+(1+HD[Layer-2])*(i-1);

        dobs[k+1]=delta[Layer-1][i];
        for(j=k+2; j<=k+HD[Layer-2]+1; j++) dobs[j]=delta[Layer-1][i]*ps[Layer-2][j-k-1];
       }
     }
  else{
   for(i=1; i<=OUT_UNIT; i++){
     k=m+(1+HD[Layer-2]+HD[0])*(i-1);

     dobs[k+1]=delta[Layer-1][i];
     for(j=k+2; j<=k+HD[Layer-2]+1; j++) dobs[j]=delta[Layer-1][i]*ps[Layer-2][j-k-1];
     for(j=k+HD[Layer-2]+2; j<=k+HD[0]+HD[Layer-2]+1; j++) dobs[j]=delta[Layer-1][i]*ps[0][j-k-HD[Layer-2]-1];
    }
  }

 free_dmatrix(ps,0,Layer,1,M);
 free_dmatrix(delta,0,Layer,1,M);

return 0;
}
   


/*** for the default case: logistic activation function ***/

int derivative();
double cost_derivative(z,d)
double *z, *d;
{
double *pout,*dobs, sum,sum1,sum2,sum3,sum4,max,min;
int i, j,k, m1, m2, m3;
FILE *ins;

dobs=dvector(1,dim);

max=max_vector(z,dim);
min=min_vector(z,dim);
if(max>50.0 || min<-50.0) return 1.0e+100;

pout=dvector(1,OUT_UNIT);

for(j=1; j<=dim; j++) d[j]=0.0; 
 
  for(sum=0.0, i=1; i<=data_num; i++){
     derivative(data_mat[i],data_y[i],z,pout,dobs);
     for(j=1; j<=dim; j++) d[j]+=dobs[j]; 
     
     for(j=1; j<=OUT_UNIT; j++)
         sum+=0.5*(data_y[i][j]-pout[j])*(data_y[i][j]-pout[j]);
    }
printf("cost=%g\n", sum);

free_dvector(pout,1,OUT_UNIT);
free_dvector(dobs,1,dim);

return sum;
}



int derivative(ox,oy,z,pout,dobs)
double *ox,*oy,*z, *pout,*dobs; 
{
int i, j, k, L, m, M;
double **ps, ave;
double **delta;
 
M=HD[0];
for(j=1; j<Layer; j++){
    if(M<HD[j]) M=HD[j];
   }
ps=dmatrix(0,Layer,1,M); 
delta=dmatrix(0,Layer,1,M); 

for(j=1; j<=P; j++) ps[0][j]=ox[j];


/* calculate the outputs of hidden layers */
m=0;
for(L=1; L<Layer-1; L++){
    
  for(i=1; i<=HD[L]; i++){
      k=m+(HD[L-1]+1)*(i-1);
      ps[L][i]=z[k+1];
      for(j=k+2; j<=k+HD[L-1]+1; j++) ps[L][i]+=ps[L-1][j-k-1]*z[j];
      ps[L][i]=LOGISTIC(ps[L][i]);
    }
  m+=(HD[L-1]+1)*HD[L];
 }


/* calculate the predicted mean from the mean unit */
if(shortcut==0){
   for(i=1; i<=OUT_UNIT; i++){
      k=m+(1+HD[Layer-2])*(i-1);
      ave=z[k+1];
      for(j=k+2; j<=k+HD[Layer-2]+1; j++) ave+=ps[Layer-2][j-k-1]*z[j];
      // pout[i]=LOGISTIC(ave);
      pout[i]=ave;
     }
   m+=(1+HD[Layer-2])*OUT_UNIT;
  }
 else{
   for(i=1; i<=OUT_UNIT; i++){
     k=m+(1+HD[Layer-2]+HD[0])*(i-1);
     ave=z[k+1];
     for(j=k+2; j<=k+HD[Layer-2]+1; j++) ave+=ps[Layer-2][j-k-1]*z[j];
     for(j=k+HD[Layer-2]+2; j<=k+HD[0]+HD[Layer-2]+1; j++) ave+=ps[0][j-k-HD[Layer-2]-1]*z[j];
     // pout[i]=LOGISTIC(ave);
     pout[i]=ave;
    }
   m+=(1+HD[0]+HD[Layer-2])*OUT_UNIT; 
  }
 
  

/****** calculation of delta ********/ 
 
 for(i=1; i<=OUT_UNIT; i++) delta[Layer-1][i]=(pout[i]-oy[i]);
 for(L=Layer-2; L>=1; L--){

   if(L==Layer-2 && shortcut==1) m-=(1+HD[0]+HD[L])*HD[L+1];
      else m-=(1+HD[L])*HD[L+1]; 
     
   for(i=1; i<=HD[L]; i++){
       delta[L][i]=0.0; 
       for(j=1; j<=HD[L+1]; j++) delta[L][i]+=z[m+(j-1)*(HD[L]+1)+i+1]*delta[L+1][j];
       delta[L][i]*=ps[L][i]*(1-ps[L][i]); 
     }
   }

/****** calculation of gradient  **************/
 m=0;
 for(L=1; L<=Layer-2; L++){
   for(i=1; i<=HD[L]; i++){
      k=m+(HD[L-1]+1)*(i-1);
      dobs[k+1]=delta[L][i];
      for(j=k+2; j<=k+HD[L-1]+1; j++) dobs[j]=delta[L][i]*ps[L-1][j-k-1];
     }
   m+=(HD[L-1]+1)*HD[L];
  }

 if(shortcut==0){
    for(i=1; i<=OUT_UNIT; i++){
        k=m+(1+HD[Layer-2])*(i-1);

        dobs[k+1]=delta[Layer-1][i];
        for(j=k+2; j<=k+HD[Layer-2]+1; j++) dobs[j]=delta[Layer-1][i]*ps[Layer-2][j-k-1];
       }
     }
  else{
   for(i=1; i<=OUT_UNIT; i++){
     k=m+(1+HD[Layer-2]+HD[0])*(i-1);

     dobs[k+1]=delta[Layer-1][i];
     for(j=k+2; j<=k+HD[Layer-2]+1; j++) dobs[j]=delta[Layer-1][i]*ps[Layer-2][j-k-1];
     for(j=k+HD[Layer-2]+2; j<=k+HD[0]+HD[Layer-2]+1; j++) dobs[j]=delta[Layer-1][i]*ps[0][j-k-HD[Layer-2]-1];
    }
  }

 free_dmatrix(ps,0,Layer,1,M);
 free_dmatrix(delta,0,Layer,1,M);

return 0;
}
   
