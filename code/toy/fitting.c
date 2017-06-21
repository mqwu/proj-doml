// #define LOGISTIC(z) 1.0/(1.0+exp(-z))

/* calculate the fitted/predicted value for each observation data (ox[1...p], oy)*/
int fitting(z,fit)
double *z, **fit;
{
int i, j, k;
double ps1[hd1+1],ave, *ps2; 


ps2=dvector(1,OUT_UNIT);

 for(i=1; i<=data_num; i++){
     fden(data_mat[i],z,ps2);
     for(j=1; j<=OUT_UNIT; j++) fit[i][j]=ps2[j];
    }

free_dvector(ps2,1,OUT_UNIT);

return 0;
}
 
