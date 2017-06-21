/* The cost function of the conditional density network */

#define LOGISTIC(z) 1.0/(1.0+exp(-z))

int fden();
double cost(z)
double *z;
{
double *ps2,sum,sum1,sum2,sum3,sum4,max,min; 
int i, j,k, m1, m2, m3;
FILE *ins;

max=max_vector(z,dim);
min=min_vector(z,dim);
if(max>50.0 || min<-50.0) return 1.0e+100;

ps2=dvector(1,OUT_UNIT);

 for(sum=0.0, i=1; i<=data_num; i++){
     fden(data_mat[i],z,ps2);
     for(j=1; j<=OUT_UNIT; j++)
         sum+=(data_y[i][j]-ps2[j])*(data_y[i][j]-ps2[j]);
    }

free_dvector(ps2,1,OUT_UNIT);


return sum;
}



/* calculate the likelihood for each observation data (ox[1...p], oy)*/
int fden(ox,z,ps2)
double *ox,*z,*ps2;
{
int i, j, k;
double ps1[hd1+1],ave; 

/* calculate the output of the first hidden layer */
for(i=1; i<=hd1; i++){
    k=(P+1)*(i-1);
    ps1[i]=z[k+1];
    for(j=k+2; j<=k+P+1; j++) ps1[i]+=ox[j-k-1]*z[j];
    ps1[i]=LOGISTIC(ps1[i]);
  }

/* calculate the predicted mean from the mean unit */
if(shortcut==0){
   for(i=1; i<=OUT_UNIT; i++){
      k=(P+1)*hd1+(1+hd1)*(i-1);
      ave=z[k+1];
      for(j=k+2; j<=k+hd1+1; j++) ave+=ps1[j-k-1]*z[j];
      // ps2[i]=LOGISTIC(ave);
      ps2[i]=ave;
     }
  }
 else{
   for(i=1; i<=OUT_UNIT; i++){
     k=(P+1)*hd1+(1+hd1+P)*(i-1);
     ave=z[k+1];
     for(j=k+2; j<=k+hd1+1; j++) ave+=ps1[j-k-1]*z[j];
     for(j=k+hd1+2; j<=k+P+hd1+1; j++) ave+=ox[j-k-hd1-1]*z[j];
     // ps2[i]=LOGISTIC(ave);
     ps2[i]=ave;
    }
  }

return 0;
}
 
