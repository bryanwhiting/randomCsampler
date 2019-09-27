#include <stdlib.h>            // malloc
#include <stdio.h>             // printf
#include <math.h>              // fabs, sqrt, etc.
#include <time.h>              // time
#include <unistd.h>            // getpid
#include <gsl/gsl_rng.h>       // GNU Scientific Library
#include <gsl/gsl_cdf.h>       // GNU Scientific Library
#include <gsl/gsl_randist.h>   // GNU Scientific Library

double dvonmises(double* x, double mu, double nu,
               double kappa1, double kappa2, double lambda) {
  double phi = x[0];
  double psi = x[1];
  double kernal = exp( kappa1*cos(phi - mu) + kappa2*cos(psi-nu) +
  lambda*sin(phi-mu)*sin(psi-nu));
  return kernal; 
}

void rvonmises2(int* ptrn, double* ptrmu, double* ptrnu,
              double* ptrkappa1, double* ptrkappa2, double* ptrlambda,
              double* phiout, double* psiout)
{
  // Assign the variables
  int n = *ptrn;
  double mu = *ptrmu;
  double nu = *ptrnu;
  double kappa1 = *ptrkappa1;
  double kappa2 = *ptrkappa2;
  double lambda = *ptrlambda;
  double modeh;
  // If Unimodal...
   if(lambda*lambda<(kappa1*kappa2)) { 
     double x[] = {mu, nu};
     modeh = dvonmises(x,mu,nu,kappa1,kappa2,lambda);
   }
   //If Bimodal...
   if(lambda*lambda>(kappa1*kappa2)) {
     double psi =
     acos(kappa2/fabs(lambda)*sqrt((lambda*lambda+kappa1*kappa1)/(lambda*lambda+kappa2*kappa2)));
     double phi =
     acos(kappa1/fabs(lambda)*sqrt((lambda*lambda+kappa2*kappa2)/(lambda*lambda+kappa1*kappa1)));
     //Caclulate the 4 possible maximum if lambda >0
     if(lambda>0.0){
       double x1[] = {mu+psi , nu+phi};
       double x2[] = {mu-psi , nu-phi};
       double mode1 = dvonmises(x1,mu,nu,kappa1,kappa2,lambda);
       double mode2 = dvonmises(x2,mu,nu,kappa1,kappa2,lambda);
       modeh = fmax(mode1,mode2);
     }
     if(lambda<0.0){
       double x1[] = {mu-phi , nu+psi};
       double x2[] = {mu+psi , nu-psi};
       double mode1 = dvonmises(x1,mu,nu,kappa1,kappa2,lambda);
       double mode2 = dvonmises(x2,mu,nu,kappa1,kappa2,lambda);
       modeh = fmax(mode1,mode2);
     }
   }
  //implement the rejection sampler now that we can use the mode to get alpha
  double pi = acos(-1);
  double height_of_sample_density = 1/(4*pi*pi);
  double alpha = height_of_sample_density/modeh;
  //printf("%s","mode height\n");
  //printf("%f\n",modeh);
  
  gsl_rng* r = gsl_rng_alloc (gsl_rng_taus);
  int seed = time(NULL)*getpid();
  gsl_rng_set(r,seed);
  
  int i = 0;
  while( i < n){
    //sample from the sample envelope: the uniform
    double two_uniforms[0];
    //evaluate the density at those uniforms
    two_uniforms[0] = gsl_ran_flat(r, -pi, pi);
    two_uniforms[1] = gsl_ran_flat(r, -pi, pi);
    double f_x = dvonmises(two_uniforms,mu,nu,kappa1,kappa2,lambda);
    double e_x = height_of_sample_density/alpha;

    // Test if U <= f_x/e_x and store the uniforms if they work
    double U = gsl_ran_flat(r,0,1);
    if( U <= f_x/e_x){
      phiout[i] = two_uniforms[0]; 
      psiout[i] = two_uniforms[1];
      i +=1;  
    }
  }
}


int main(int argc, char* argv[]) {
  if ( argc != 7 ) {
    fprintf(stderr,"usage: %s nreps mu nu kappa1 kappa2 lambda\n",argv[0]);
   
    return 1;
  }
  if (argv[7]=0){
    fprintf(stderr,"lambda cannot be 0");
    return 1;
  }
  // Simulation parameters
  int n = atoi(argv[1]);
  double mu = atof(argv[2]);
  double nu = atof(argv[3]);
  double kappa1 = atof(argv[4]);
  double kappa2 = atof(argv[5]);
  double lambda = atof(argv[6]);
//because our rvonmises funciton is void, we must allocate memor to psi and phi
//so that we can get results and call them later
  double* phiout = (double*) malloc(n*sizeof(double));
  double* psiout = (double*) malloc(n*sizeof(double));
  rvonmises2(&n,&mu,&nu,&kappa1,&kappa2,&lambda,phiout,psiout);
  printf("%s %s\n","phi","psi");
  for(int i = 0; i < n; i++){
    printf("%f %f\n",phiout[i],psiout[i]);
  }
  free(phiout);
  free(psiout);
  return 0;//to get it to work, comment this out
}

