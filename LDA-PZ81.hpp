#ifndef LDA_hpp
#define LDA_hpp

#include <vector>
#include <cmath>

/* energy/potential functions of LDA */
double e_xc(double r_s){
  if(r_s < 1.0){
    return - 0.458165/r_s - 0.0480 + 0.0311*log(r_s) - 0.0116*r_s + 0.0020*r_s*log(r_s);
    // -3.0/4.0 * pow(9.0/(4.0*pow(M_PI,2)), 1.0/3.0) * 1/r_s - 0.0480 + 0.0311*log(r_s) - 0.0116*r_s + 0.0020*r_s*log(r_s)
  }
  else{
    return - 0.458165/r_s - 0.1423/(1.0 + 1.0529*pow(r_s, 0.5) + 0.3334*r_s);
    // -3.0/4.0 * pow(9.0/(4.0*pow(M_PI,2)), 1.0/3.0) * 1/r_s  - 0.1423/(1.0 + 1.0529*pow(r_s, 0.5) + 0.3334*r_s)
  }
}

// exchange-correlation potential
double V_xc(double r_s){
  if(r_s < 1.0){
    return -0.610887/r_s + 1.0/3.0 * (-0.1751 + 0.0933*log(r_s) - 0.0252*r_s + 0.0040*r_s*log(r_s));
  }
  else{
    return -0.610887/r_s + 1.0/3.0 * (-0.4269/(1.0 + 1.0529*pow(r_s, 0.5) + 0.3334*r_s) - 0.1423 * r_s * (1.0529/(2.0*pow(r_s, 0.5)) + 0.3334)/pow(1.0 + 1.0529*pow(r_s, 0.5) + 0.3334*r_s, 2));
  }
}

/* CORRELATION PART ONLY */
// These functionals are used to check program using the example of He atom in Thijssen.
/*double e_xc(double r_s){
  return -0.458/r_s;
}

double V_xc(double r_s){
  return -0.458 * 4.0/(3.0 * r_s);
}*/
#endif /* LDA_hpp */
