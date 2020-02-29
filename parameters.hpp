#ifndef parameters_hpp
#define parameters_hpp

#include <vector>

/* simulation settings */
// setting of arrays
static double r_max = 50.0;
static int num_steps = 10000;
static double h = r_max/num_steps;

// setting of SCF loop
static const int max_SCF_loop = 200;
static const double SCF_convergence = 5.0e-8;


// setting of DFT loop
static long int max_DFT_loop = 1000;
static const double DFT_convergence = 1.0e-8;

// other parameters
static double n_4pi_minimum = 1.0e-30; // the minimum value of n_4pi : if n_4pi < n_4pi_minimum, we use this minimum value when calculating r_s


/* initialization */
double u_1s_init(double r, int Z){
  return 2.0 * pow(Z, 1.5) * r * exp(-Z*r);
}
double u_2s_init(double r, int Z){
  return 1.0/(2.0*sqrt(2.0)) * pow(Z, 1.5) * r * (2.0-Z*r) * exp(-Z*r/2.0);
}
double u_2p_init(double r, int Z){
  return 1.0/(2.0*sqrt(6.0)) * pow(Z, 2.5) * r*r * exp(-Z*r/2.0);
}
double u_3s_init(double r, int Z){
  return 2.0/(81.0*sqrt(3.0)) * pow(Z, 1.5) * r*(27.0 - 18.0*Z*r + 2.0*pow(Z*r, 2.0)) * exp(-Z*r/3.0);
}
double u_3p_init(double r, int Z){
  return 4.0/(81.0*sqrt(6.0)) * pow(Z, 2.5) * r*r * (6.0 - Z*r) * exp(-Z*r/3.0);
}
/*double u_3d_init(double r, int Z){
  return 4.0/(81.0*sqrt(30.0)) * pow(Z, 3.5) * pow(r, 3.0) * exp(-Z*r/3.0);
}*/
#endif /* parameters_hpp */
