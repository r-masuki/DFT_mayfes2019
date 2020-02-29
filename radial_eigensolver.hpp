#ifndef radial_eigensolver_hpp
#define radial_eigensolver_hpp

#include <vector>
#include "parameters.hpp"
#include "LDA-PZ81.hpp"

double calculate_u0(double E, int Z, std::vector<double> &n_4pi, std::vector<double> &rV_Coulomb, std::vector<double> &r, int l, std::vector<double> &u_tmp);
double solve_radial_eigen_eq(double E_left, double search_step, int Z, std::vector<double> &n_4pi, std::vector<double> &rV_Coulomb, std::vector<double> &r, int l, std::vector<double> &u_tmp);

double calculate_u0(double E, int Z, std::vector<double> &n_4pi, std::vector<double> &rV_Coulomb, std::vector<double> &r, int l, std::vector<double> &u_tmp){
  // parameter
  double r_s;
  // preparation
  u_tmp[num_steps-1] = exp(-r_max+h)* 1.0e-80;
  u_tmp[num_steps-2] = exp(-r_max+2.0*h)* 1.0e-80;
  
  // solve radial equation
  for(int i = num_steps-3; i >= 0; i--){
    r_s = pow(3.0/std::max(n_4pi[i+1], n_4pi_minimum), 1.0/3.0);
    u_tmp[i] = 2.0 * u_tmp[i+1] - u_tmp[i+2] + h*h* 2.0 *(-Z/r[i+1] + rV_Coulomb[i+1]/r[i+1] + V_xc(r_s) + l*(l+1)/(2.0*pow(r[i+1], 2)) - E)*u_tmp[i+1];
  }
  
  return u_tmp[0];
}

double solve_radial_eigen_eq(double E_left, double search_step, int Z, std::vector<double> &n_4pi, std::vector<double> &rV_Coulomb, std::vector<double> &r, int l, std::vector<double> &u_tmp)
{
  double E_right, E_mid, u0_left, u0_right, u0_mid;
  
  E_right = E_left;
  u0_left = calculate_u0(E_left, Z, n_4pi, rV_Coulomb, r, l, u_tmp);
  while(std::isnan(u0_left)){
    E_left += search_step;
    u0_left = calculate_u0(E_left, Z, n_4pi, rV_Coulomb, r, l, u_tmp);
  }
  u0_right = u0_left;
  
  while(u0_left*u0_right > 0.0){
    E_right += search_step;
    u0_right = calculate_u0(E_right, Z, n_4pi, rV_Coulomb, r, l, u_tmp);
  }
  // bisection method
  for(int DFT_loop = 0; DFT_loop < max_DFT_loop; DFT_loop++){
    // print section
    std::cerr << "DFT_loop: " << DFT_loop << std::scientific << " section: [" << E_left << ", " << E_right << "]" << std::endl;
    // calculate u(0) of E_mid
    E_mid = (E_right + E_left)/2.0;
    u0_mid = calculate_u0(E_mid, Z, n_4pi, rV_Coulomb, r, l, u_tmp);
    // update section
    if(u0_mid == 0.0){
      break;
    }
    else if(u0_mid * u0_left > 0.0){
      E_left = E_mid;
      u0_left = u0_mid;
    }
    else{
      E_right = E_mid;
      u0_right = u0_mid;
    }
    // check convergence
    if(E_right - E_left < DFT_convergence){
      std::cerr << "The radial eigenequation for 1s-orbital has converged:" << std::endl;
      break;
    }
  }
  
  return E_mid;
}

#endif /* radial_eigensolver_hpp */
