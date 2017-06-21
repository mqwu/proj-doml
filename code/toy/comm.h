#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
 
#define pi 3.14159265
 
#define data_num  194
 
#define P 2            /* the number of input units */
#define hd1 30         /* 15 the number of first hidden layer units */
#define OUT_UNIT 1     /* the number of output units */
static int shortcut=0;

#define total_iteration 10000000
#define stepscale  10000
#define WARM           1
#define scale          5 
#define N              1         /* population size */

static double lowE=0.0, maxEE=50.0, maxE=50.0, range=5.0;
static double tau=0.6, rho=1.0, tem=1.0, delta,stepsize; 
static double **solution,*min,*refden,**data_mat, **data_y, VARIANCE;
static double accept_loc, total_loc,accept_dir, total_dir;
static int dim, sze, Best=5;
