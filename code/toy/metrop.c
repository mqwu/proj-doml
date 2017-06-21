int metropolis(x, fx)
double *x, *fx;
{
int i, j, k, accept;
double *y,*e,fy,r, un, sum, max;

y=dvector(1,dim);
e=dvector(1,dim);

   uniform_direction(e,dim);
   for(k=1; k<=dim; k++) y[k]=x[k]+e[k]*gasdev()*0.1;
   fy=cost(y);

   r=-fy+(*fx);
	    
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
        }

free_dvector(y,1,dim);
free_dvector(e,1,dim);

return 0;
}
