
/* calculate log(Gamma(s))  */
static double Gamma_table[10001];
double loggamma(s)
double s;
{
double value, ss;
int k;

value=.0;  ss=s;
while(ss>1.0){
   ss=ss-1.0;
   value+=log(ss); }

k=(int)(ss*10000);

value+=log((k+1-ss*10000)*Gamma_table[k]+(ss*10000-k)*Gamma_table[k+1]);

return value;
}


/* calculate log(k!) */
double logpum(k)
int k;
{
double value;
int i;

for(value=0.0, i=1; i<=k; i++) value+=log(1.0*i);

return value;
}


/* generate the random variable form Gamma(a,b) */
double Rgamma(a,b)
double a, b;
{
int ok;
double d, q, un, u1, y, z; 

if(a<=0.0 || b<=0.0) { printf("Gamma parameter error (<0.0)\n"); return; }

if(a<1.0){  /* Ahrens, P.213 */
 ok=0;
while(ok==0){ 
  un=0.0;
  while(un<=0.0 || un>=1.0) un=rand()*1.0/RAND_MAX;
  d=(2.718282+a)/2.718282;
  q=d*un;
  
  if(q<=1.0){ 
    z=exp(1.0/a*log(q));
    u1=rand()*1.0/RAND_MAX;
    if(u1<exp(-z)) ok=1;
             }
   else{
     z=-1.0*log((d-q)/a); 
     u1=rand()*1.0/RAND_MAX;
     if(u1<exp((a-1)*log(z))) ok=1;
       }
            } /* end ok */
   }
 else {  /* a>=1.0 Fishman, P.214 */
  ok=0;
  while(ok==0){
    un=0.0;
    while(un<=0.0 || un>=1.0) un=rand()*1.0/RAND_MAX;
    y=-1.0*log(un);

    u1=rand()*1.0/RAND_MAX;
    if(u1<exp((a-1)*(log(y)-(y-1)))) { z=a*y; ok=1; }
             }
      }

return z/b;
}


/* calculated the log-density of  z~gamma(a,b) */
double dgamma(z,a,b)
double z,a,b;
{
double logcon, den;

logcon=loggamma(a);
den=log(b)-b*z+(a-1)*log(b*z)-logcon;

return den;
}



double dloggauss(z,mu,var)
double z,mu,var;
{
double sum;

sum=-0.5*log(2.0*pi*var);
sum+=-0.5*(z-mu)*(z-mu)/var;

return sum;
}


double dlogstudent(z,k) 
double z;
int k;  /* the degree of freedom */
{
double logprob;

logprob=-0.5*(k+1)*log(1.0+z*z/k);

return logprob;
}


double gasdev()
{
        static int iset=0;
        static double gset;
        double fac,r,v1,v2;

        if  (iset == 0) {
                do {
                        v1=rand()*2.0/RAND_MAX-1.0;
                        v2=rand()*2.0/RAND_MAX-1.0;
                        r=v1*v1+v2*v2;
                } while (r >= 1.0);
                fac=sqrt(-2.0*log(r)/r);
                gset=v1*fac;
                iset=1;
                return v2*fac;
        } else {
                iset=0;
                return gset;
        }
}




int uniform_direction(d, n)
double *d;
int n;
{
double sum;
int k;
 
for(sum=0, k=1; k<=n; k++){
   d[k]=gasdev();
   sum+=d[k]*d[k];
 }
   
 for(k=1; k<=n; k++)
    d[k]=d[k]/sqrt(sum);
 
return 0;
 }

/*
int dmaxclass(z)
double *z;
{
int i, maxi;
double maxd;

maxd=z[1]; maxi=1;
for(i=2; i<=OUT_UNIT; i++)
  if(z[i]>maxd){ maxd=z[i]; maxi=i; }

return maxi;
}

int imaxclass(z)
int *z;
{
int i, maxi;
int maxd;

maxd=z[1]; maxi=1;
for(i=2; i<=OUT_UNIT; i++)
  if(z[i]>maxd){ maxd=z[i]; maxi=i; }

return maxi;
}
*/

int binary_trans(k,l,d)
int k, l,*d;
{
int i, j;

for(i=1; i<=l; i++) d[i]=0;
j=l;
while(k>0){
   d[j]=k%2;
   k=k/2;
   j--;
 } 

return 0;
}


double logsum(a,b)
        double a, b;
{
        double sum;
 
        if(a>b) sum=a+log(1.0+exp(b-a));
           else sum=b+log(1.0+exp(a-b));
 
        return sum;
}



double max_vector(x,n)
double *x;
int n;
{
double max;
int i;

  max=x[1];
  for(i=2; i<=n; i++) 
     if(x[i]>max) max=x[i];

  return max;
 }

double min_vector(x,n)
 double *x;
 int n;
{
double min;
int i;

  min=x[1];
  for(i=2; i<=n; i++)
     if(x[i]<min) min=x[i];

  return min;
 }



double sample_variance(x,n)
double *x;
int n;
{
int i;
double sum1, sum2, mean, var;

sum1=sum2=0.0;
for(i=1; i<=n; i++){
   sum1+=x[i];
   sum2+=x[i]*x[i];
  }
mean=sum1/n;
var=(sum2-n*mean*mean)/(n-1);

return var;
}
