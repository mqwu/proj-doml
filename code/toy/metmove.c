int MetmoveSA(x,fvalue,hist,tem)
double **x,*fvalue,**hist,tem;
{
int i, j,k1, k2,l,m,accept,ok, rep;
double *y,*e,fy,r, un, dw, dwde,maxe,logrho,var1,var2;
FILE *ins;

y=dvector(1,dim);
e=dvector(1,dim);

for(rep=1; rep<=N; rep++){
  
    if(maxE<1.0) range=2.5;
    sze=floor((maxE+range-lowE)*scale);

    if(fvalue[rep]>maxE+range) k1=sze;
       else if(fvalue[rep]<lowE) k1=0;
          else k1=floor((fvalue[rep]-lowE)*scale);
    var1=sample_variance(x[rep],dim);

    for(i=1; i<=dim; i++) y[i]=x[rep][i];
    un=rand()*1.0/RAND_MAX;
    if(un<=0.5){
       i=dim+1;
       while(i>dim) i=(int)(rand()*1.0/RAND_MAX*dim)+1;
       /*
       un=gasdev()*sqrt(var1*stepsize);
       */
       un=gasdev()*stepsize;
       y[i]=x[rep][i]+un;
      }
     else{
       uniform_direction(e,dim);
       /*
       un=gasdev()*sqrt(var1*stepsize);
       */
       un=gasdev()*stepsize;
       for(i=1; i<=dim; i++) y[i]=x[rep][i]+e[i]*un;
      }
     fy=cost(y);
     var2=sample_variance(y,dim);


     if(fy>maxE+range) k2=sze;
        else if(fy<lowE) k2=0;
          else k2=floor((fy-lowE)*scale);

     r=hist[k1][2]-hist[k2][2]-fy/tem+fvalue[rep]/tem;
     /* +0.5*log(var1/var2)-0.5*un*un*(1.0/var2-1.0/var1)/stepsize; */

     if(k2>=sze) accept=0;
       else{
         if(r>0.0) accept=1;
           else{
              un=rand()*1.0/RAND_MAX;
              if(un<exp(r)) accept=1;
                 else accept=0;
	      }
          }
	

     if(accept==1){ 
        for(j=1; j<=dim; j++) x[rep][j]=y[j];
        fvalue[rep]=fy;  
         
        for(i=1; i<=sze; i++){
           if(i!=k2) hist[i][2]+=-delta*refden[i];
              else hist[i][2]+=delta*(1.0-refden[i]);
          }
        hist[k2][3]+=1.0;
        accept_loc+=1.0; total_loc+=1.0;
       }
      else{  /* accept=0 */

         for(i=1; i<=sze; i++){
             if(i!=k1) hist[i][2]+=-delta*refden[i];
                else hist[i][2]+=delta*(1.0-refden[i]);
              }
         hist[k1][3]+=1.0;
         total_loc+=1.0;
       }


      if(accept==1){
         j=1; ok=1;
         while(j<=Best && fy<min[j]) j++;
         if(fy==min[j]) ok=0;
         j--;
         if(j>=1 && ok==1){
            for(l=1; l<j; l++){
                for(m=1; m<=dim; m++) solution[l][m]=solution[l+1][m];
                min[l]=min[l+1];
               }
            for(m=1; m<=dim; m++) solution[j][m]=y[m];
            min[j]=fy;
          }
        maxe=max_vector(fvalue,N); 
        if(maxe<maxE) maxE=maxe;
      }
}

free_dvector(y,1,dim);
free_dvector(e,1,dim);


return 0;
}
