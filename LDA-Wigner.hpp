#ifndef LDA_hpp
#define LDA_hpp

#include <vector>
#include <cmath>

/* energy/potential functions of LDA */
double e_xc(double r_s){
  return -0.458/r_s - 0.44/(r_s+7.8);
}

// exchange-correlation potential
double V_xc(double r_s){
  return -0.458 * 4.0/(3.0 * r_s) - 0.44/(7.8+r_s) - 0.44 * r_s/(3.0 * pow(r_s+7.8, 2));
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
