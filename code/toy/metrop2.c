int metropolis_post(x, fx)
double *x, *fx;
{
int i, j, k, l,m,ok,accept, regiony, iter,Best=5, M=0;
double *y,*e,**solution,*min,fy,r, un, sum, max;
double amet, tmet, t=1.0e-7;
FILE *ins;

y=dvector(1,dim);
e=dvector(1,dim);
solution=dmatrix(0,Best,0,dim);
min=dvector(0, Best);

for(i=1; i<=Best; i++){
   min[i]=*fx;
   for(j=1; j<=dim; j++) solution[i][j]=x[j];
  }

for(i=1; i<=dim; i++) y[i]=x[i];
amet=tmet=1;

for(iter=1; iter<M; iter++){

  un=rand()*1.0/RAND_MAX;
  if(un<=0.5){
     i=dim+1;
     while(i>dim) i=(int)(rand()*1.0/RAND_MAX*dim)+1;
     y[i]=x[i]+0.005*gasdev();
    }
   else{
     uniform_direction(e,dim);
     un=0.01*gasdev();
     for(i=1; i<=dim; i++) y[i]=x[i]+e[i]*un;
    }
  fy=cost(y);

  r=(-fy+(*fx))/t;
	    
  if(r>0.0) accept=1;
     else{
        un=0.0;
        while(un<=0.0) un=rand()*1.0/RAND_MAX;
        if(un<=exp(r)) accept=1;
           else accept=0;
       }

  if(accept==1){ 
     for(k=1; k<=dim; k++) x[k]=y[k];
     *fx=fy;
     amet+=1.0;
    
     j=1; ok=1;
     while(j<=Best && (*fx)<min[j]) j++;
     if(*fx==min[j]) ok=0;
     j--;
     if(j>=1 && ok==1){
         for(l=1; l<j; l++){
             for(m=1; m<=dim; m++) solution[l][m]=solution[l+1][m];
              min[l]=min[l+1];
             }
         for(m=1; m<=dim; m++) solution[j][m]=x[m];
         min[j]=*fx;
        }
     }
  tmet+=1.0;
  
  if(iter%10000==0)  printf("iter=%d  min=%13.7f\n", iter/10000, min[Best]);
}

/*
ins=fopen("uu0.log","a");
if(ins==NULL){ printf("can't write to file\n"); return 1; }
for(i=1; i<=Best; i++){
   fprintf(ins, "energy=%15.8f\n", min[i]);
   for(j=1; j<=dim; j++) fprintf(ins, " %g", solution[i][j]);
   fprintf(ins, "\n");
  }
fclose(ins);
*/

*fx=min[Best];
for(j=1; j<=dim; j++) x[j]=solution[Best][j];


printf("acceptance rate=%g\n", 1.0*amet/tmet);

free_dvector(y,1,dim);
free_dvector(e,1,dim);
free_dmatrix(solution,0,Best,0,dim);
free_dvector(min,0,Best);

return 0;
}
