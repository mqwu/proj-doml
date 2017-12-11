#define LOGISTIC(z) 1.0/(1.0+exp(-z))
#define TANH(z) (1.0-exp(-2*z))/(1.0+exp(-2*z))
#define RELU(z) (z>0)?z:0 


int input_derivative_tanh(ox,z,dobs)
double *ox,*z,*dobs; 
{
int i, j, k, L, m, M, **node, hindex;
double **ps, ave, **D, un;
double **delta;
 
M=HD[0]; hindex=0;
for(j=1; j<Layer; j++){
    if(M<HD[j]) M=HD[j]; 
    hindex+=HD[j]; 
   }

ps=dmatrix(0,Layer,1,M); 
delta=dmatrix(0,Layer,1,M); 
D=dmatrix(1,hindex,1,P); 
node=imatrix(0,Layer,1,M); 

for(k=0, L=1; L<Layer; L++)
   for(i=1; i<=HD[L]; i++){
     k++;
     node[L][i]=k; 
   } 

for(i=1; i<=hindex; i++) 
  for(j=1; j<=P; j++) D[i][j]=0.0;
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
  for(i=1; i<=OUT_UNIT; i++){
      k=m+(1+HD[Layer-2])*(i-1);
      ave=z[k+1];
      for(j=k+2; j<=k+HD[Layer-2]+1; j++) ave+=ps[Layer-2][j-k-1]*z[j];
      ps[Layer-1][i]=ave;
     }
  
/******* calculation of derivative with respect to inputs ****/
 
// first hidden layer
 L=1; m=0;
 for(i=1; i<=HD[L]; i++){
    un=(1+ps[L][i])*(1-ps[L][i]);
    for(j=m+2; j<=m+HD[L-1]+1; j++) D[node[L][i]][j-m-1]=un*z[j]; 
    m+=HD[L-1]+1; 
   }

// other hidden layer
 for(L=2; L<Layer-1; L++)
   for(i=1; i<=HD[L]; i++){
      un=(1+ps[L][i])*(1-ps[L][i]);
      for(k=1; k<=HD[0]; k++)  
         for(j=m+2; j<=m+HD[L-1]+1; j++)
             D[node[L][i]][k]+=un*z[j]*D[node[L-1][j-m-1]][k];
      m+=HD[L-1]+1;
    }

// output layer 
L=Layer-1;
for(i=1; i<=OUT_UNIT; i++){
    // un=(1+ps[L][i])*(1-ps[L][i]);
    un=1.0; 
    for(k=1; k<=HD[0]; k++)  
       for(j=m+2; j<=m+HD[L-1]+1; j++)
             D[node[L][i]][k]+=un*z[j]*D[node[L-1][j-m-1]][k];
      m+=HD[L-1]+1;
    }

       
 for(k=1; k<=P; k++) dobs[k]=D[node[Layer-1][1]][k];
 
 
 free_dmatrix(ps,0,Layer,1,M);
 free_dmatrix(delta,0,Layer,1,M);
 free_dmatrix(D,1,hindex,1,P);
 free_imatrix(node,0,Layer,1,M);


return 0;
}

 

int input_derivative_relu(ox,z,dobs)
double *ox,*z,*dobs; 
{
int i, j, k, L, m, M, **node, hindex;
double **ps, ave, **D, un;
double **delta;
 
M=HD[0]; hindex=0;
for(j=1; j<Layer; j++){
    if(M<HD[j]) M=HD[j]; 
    hindex+=HD[j]; 
   }

ps=dmatrix(0,Layer,1,M); 
delta=dmatrix(0,Layer,1,M); 
D=dmatrix(1,hindex,1,P); 
node=imatrix(0,Layer,1,M); 

for(k=0, L=1; L<Layer; L++)
   for(i=1; i<=HD[L]; i++){
     k++;
     node[L][i]=k; 
   } 

for(i=1; i<=hindex; i++) 
  for(j=1; j<=P; j++) D[i][j]=0.0;
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
  for(i=1; i<=OUT_UNIT; i++){
      k=m+(1+HD[Layer-2])*(i-1);
      ave=z[k+1];
      for(j=k+2; j<=k+HD[Layer-2]+1; j++) ave+=ps[Layer-2][j-k-1]*z[j];
      ps[Layer-1][i]=ave;
     }
  
/******* calculation of derivative with respect to inputs ****/
 
// first hidden layer
 L=1; m=0;
 for(i=1; i<=HD[L]; i++){
    if(ps[L][i]>0.0) un=1.0;
       else un=0.0; 
    for(j=m+2; j<=m+HD[L-1]+1; j++) D[node[L][i]][j-m-1]=un*z[j]; 
    m+=HD[L-1]+1; 
   }

// other hidden layers
 for(L=2; L<Layer-1; L++)
   for(i=1; i<=HD[L]; i++){
      if(ps[L][i]>0.0) un=1.0;
         else un=0.0;
      for(k=1; k<=HD[0]; k++)  
         for(j=m+2; j<=m+HD[L-1]+1; j++)
             D[node[L][i]][k]+=un*z[j]*D[node[L-1][j-m-1]][k];
      m+=HD[L-1]+1;
    }

// output layer 
L=Layer-1;
for(i=1; i<=OUT_UNIT; i++){
    // un=(1+ps[L][i])*(1-ps[L][i]);
    un=1.0;
    for(k=1; k<=HD[0]; k++)  
       for(j=m+2; j<=m+HD[L-1]+1; j++)
             D[node[L][i]][k]+=un*z[j]*D[node[L-1][j-m-1]][k];
      m+=HD[L-1]+1;
    }

 for(k=1; k<=P; k++) dobs[k]=D[node[Layer-1][1]][k];
 
 
 free_dmatrix(ps,0,Layer,1,M);
 free_dmatrix(delta,0,Layer,1,M);
 free_dmatrix(D,1,hindex,1,P);
 free_imatrix(node,0,Layer,1,M);


return 0;
}

 


/*** for the default case: logistic activation function ***/

int input_derivative(ox,z,dobs)
double *ox,*z,*dobs; 
{
int i, j, k, L, m, M, **node, hindex;
double **ps, ave, **D, un;
double **delta;
 
M=HD[0]; hindex=0;
for(j=1; j<Layer; j++){
    if(M<HD[j]) M=HD[j]; 
    hindex+=HD[j]; 
   }

ps=dmatrix(0,Layer,1,M); 
delta=dmatrix(0,Layer,1,M); 
D=dmatrix(1,hindex,1,P); 
node=imatrix(0,Layer,1,M); 

for(k=0, L=1; L<Layer; L++)
   for(i=1; i<=HD[L]; i++){
     k++;
     node[L][i]=k; 
   } 

for(i=1; i<=hindex; i++) 
  for(j=1; j<=P; j++) D[i][j]=0.0;
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
  for(i=1; i<=OUT_UNIT; i++){
      k=m+(1+HD[Layer-2])*(i-1);
      ave=z[k+1];
      for(j=k+2; j<=k+HD[Layer-2]+1; j++) ave+=ps[Layer-2][j-k-1]*z[j];
      ps[Layer-1][i]=ave;
     }
  
/******* calculation of derivative with respect to inputs ****/
 
// first hidden layer
 L=1; m=0;
 for(i=1; i<=HD[L]; i++){
    un=ps[L][i]*(1-ps[L][i]);
    for(j=m+2; j<=m+HD[L-1]+1; j++) D[node[L][i]][j-m-1]=un*z[j]; 
    m+=HD[L-1]+1; 
   }

// other hidden layers
 for(L=2; L<Layer-1; L++)
   for(i=1; i<=HD[L]; i++){
      un=ps[L][i]*(1-ps[L][i]);
      for(k=1; k<=HD[0]; k++)  
         for(j=m+2; j<=m+HD[L-1]+1; j++)
             D[node[L][i]][k]+=un*z[j]*D[node[L-1][j-m-1]][k];
      m+=HD[L-1]+1;
    }

 // output layer 
 L=Layer-1;
 for(i=1; i<=OUT_UNIT; i++){
    // un=(1+ps[L][i])*(1-ps[L][i]);
    un=1.0;
    for(k=1; k<=HD[0]; k++)  
       for(j=m+2; j<=m+HD[L-1]+1; j++)
             D[node[L][i]][k]+=un*z[j]*D[node[L-1][j-m-1]][k];
      m+=HD[L-1]+1;
    }

       
 for(k=1; k<=P; k++) dobs[k]=D[node[Layer-1][1]][k];
 
 free_dmatrix(ps,0,Layer,1,M);
 free_dmatrix(delta,0,Layer,1,M);
 free_dmatrix(D,1,hindex,1,P);
 free_imatrix(node,0,Layer,1,M);


return 0;
}

 

