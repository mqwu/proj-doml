int MetmoveSA(z,fvalue,hist,tem,pos,accept)
double *z,*fvalue,**hist,tem;
int *pos,*accept;
{
int i, j,k1, k2,l,m,ok, rep;
double *y,*e, fz, fy,r, un, dw, dwde,maxe,logrho,var1,var2;
FILE *ins;

y=dvector(1,dim);
e=dvector(1,dim);

    // sze=floor((maxEE+range-lowE)*scale); 
    fz=*fvalue;

    if(fz>maxEE+range) k1=sze;
       else if(fz<lowE) k1=0;
          else k1=floor((fz-lowE)*scale);
    var1=sample_variance(z,dim);

    for(i=1; i<=dim; i++) y[i]=z[i];
    un=rand()*1.0/RAND_MAX;
    if(un<=0.5){
       i=dim+1;
       while(i>dim) i=(int)(rand()*1.0/RAND_MAX*dim)+1;
       /*
       un=gasdev()*sqrt(var1*stepsize);
       */
       un=gasdev()*stepsize;
       y[i]=z[i]+un;
      }
     else{
       uniform_direction(e,dim);
       /*
       un=gasdev()*sqrt(var1*stepsize);
       */
       un=gasdev()*stepsize;
       for(i=1; i<=dim; i++) y[i]=z[i]+e[i]*un;
      }
     fy=cost(y);
     var2=sample_variance(y,dim);


     if(fy>maxEE+range) k2=sze;
        else if(fy<lowE) k2=0;
          else k2=floor((fy-lowE)*scale);

     r=hist[k1][2]-hist[k2][2]-fy/tem+fz/tem;
     /* +0.5*log(var1/var2)-0.5*un*un*(1.0/var2-1.0/var1)/stepsize; */

     if(r>0.0) *accept=1;
        else{
            un=rand()*1.0/RAND_MAX;
            if(un<exp(r)) *accept=1;
                else *accept=0;
           }
	

     if(*accept==1){ 
         for(j=1; j<=dim; j++) z[j]=y[j];
         fz=fy; 
         *fvalue=fz;  
         *pos=k2;
        }
      else{ *pos=k1; }
        
        
free_dvector(y,1,dim);
free_dvector(e,1,dim);


return 0;
}
