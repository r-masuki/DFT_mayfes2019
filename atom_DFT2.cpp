// headers
#include <iostream>
#include <iomanip>
#include <complex>
#include <cstdio>
#include <cmath>
#include <vector>
#include "physics.hpp"
#include "units.hpp"

// simulation setting
#include "system_settings.hpp"
#include "parameters.hpp"
#include "LDA-PZ81.hpp"
#include "radial_eigensolver.hpp"

double w = 0.5; // converger parameter

int main(int argc, char **argv){
  
  /* control variables */
  int ret;
  double normalize_u[5], normalize_u_tmp;
  double E_left, E_tmp, r_s, delta_E;
  
  /* preparations of physical variables */
  int Z;
  int N[5];
  std::vector<double> r;
  std::vector<double> u[5];
  std::vector<double> n_4pi;
  std::vector<double> rV_Coulomb;
  std::vector<double> u_tmp;
  double E[5], E_new[5], E_diff[5], E_gnd, E_xc, E_correction;
  
  // initialization
  ret = scan_Z(argv[1], &Z); // definition is written in "system_settings.hpp". assert 1 <= Z <= 18
  std::cerr << "The total number of protons/electrons: " << std::endl;
  std::cerr << "Z = " << Z << std::endl << std::endl;
  if(ret > 0){ // error when scanning Z
    std::cerr << "error: Z is out of range" << std::endl;
    return 1;
  }
  set_N(Z, N); // the definition is written in "system_settings.hpp"
  std::cerr << "The numbers of electrons in each orbitals:" << std::endl;
  std::cerr << "N_1s : " << N[0] << std::endl;
  std::cerr << "N_2s : " << N[1] << std::endl;
  std::cerr << "N_2p : " << N[2] << std::endl;
  std::cerr << "N_3s : " << N[3] << std::endl;
  std::cerr << "N_3p : " << N[4] << std::endl << std::endl;
  
  for(int i = 0; i < num_steps; i++){
    r.push_back(h*i);
    // initialization of u : the definitions are in "parameters.hpp"
    u[0].push_back(u_1s_init(h*i, Z));
    u[1].push_back(u_2s_init(h*i, Z));
    u[2].push_back(u_2p_init(h*i, Z));
    u[3].push_back(u_3s_init(h*i, Z));
    u[4].push_back(u_3p_init(h*i, Z));
    // push temporary entries to each vectors
    n_4pi.push_back(0.0);
    rV_Coulomb.push_back(0.0);
    u_tmp.push_back(0.0);
  }
  for(int i = 0; i < 5; i++){
    E[i] = 0.0;
  }
  
  /* SCF equation loop : solve Self-Consistent-Field equation */
  std::cerr << "SOLVE SCF-EQUATION:" << std::endl << std::endl;
  
  for(int SCF_loop = 0; SCF_loop < max_SCF_loop; SCF_loop++){
    //std::cerr << "SCF_loop: " << SCF_loop << std::endl; // we print this after printing E_gnd
    
    /* update current state */
    // copy E_new to E
    for(int i = 0; i < 5; i++){
      E[i] = E_new[i];
    }
    
    // normalization of u
    for(int i = 0; i < 5; i++){
      normalize_u[i] = 0.0;
      for(int j = 0; j < num_steps; j++){
        normalize_u[i] += u[i][j]*u[i][j] * h;
      }
      for(int j = 0; j < num_steps; j++){
        u[i][j] /= sqrt(normalize_u[i]);
      }
    }
    
    // calculate n_4pi
    for(int i = 0; i < num_steps; i++){
      n_4pi[i] = 0.0;
    }
    for(int i = 0; i < 5; i++){
      for(int j = 1; j < num_steps; j++){
        n_4pi[j] += N[i] * pow(u[i][j]/r[j], 2);
      }
    }n_4pi[0] = n_4pi[1]; // this is temporary : we do this because we cannot calculate u/r for r = 0. There is no problem for doing this because n_4pi[0] is used only for plotting and not used in solving radial differential equations
    
    // calculate rV_Coulomb
    rV_Coulomb[0] = 0.0; // these initial conditions are temporary.
    rV_Coulomb[1] = h;
    for(int i = 2; i < num_steps; i++){ // solve d^2/dr^2 rV_Coulomb(r) = - r*4pi*n(r)
      rV_Coulomb[i] = 2.0*rV_Coulomb[i-1] - rV_Coulomb[i-2] - h*h * r[i-1] * n_4pi[i-1];
    }
    for(int i = 1; i < num_steps; i++){ // determine uncertainty term
      rV_Coulomb[i] += (Z - rV_Coulomb[num_steps-1]) /((double)num_steps-1.0) * i;
    } // we impose rV_Coulomb(r_max) = Z : there is an correction term of exp(-r)
    
    /* print current state */
    // calculate ground state energy
    E_gnd = 0.0;
    // the sum of one particle energies
    for(int i = 0; i < 5; i++){
      E_gnd += N[i] * E[i];
    }
    // calculate exchange-correlation energy and correction energy
    E_correction = 0.0;
    for(int i = 1; i < num_steps; i++){
      r_s = pow(3.0/(std::max(n_4pi[i], n_4pi_minimum)), 1.0/3.0);
      E_correction += n_4pi[i] * r[i]*r[i] * (e_xc(r_s) - V_xc(r_s) - 0.5 * rV_Coulomb[i]/r[i]) * h;
    }
    E_gnd += E_correction;
    std::cerr << "Ground state energy:\nE_gnd: " << E_gnd << std::endl << std::endl;
    std::cerr << "SCF_loop: " << SCF_loop << std::endl;
    
    /* solve radial equation */
    // solve 1s-orbital
    E_left = -0.5*Z*Z;
    std::cerr << "SOLVE RADIAL EIGENEQUATION for 1s-orbital:" << std::endl;
    E_tmp = solve_radial_eigen_eq(E_left, 0.01, Z, n_4pi, rV_Coulomb, r, 0, u_tmp); // search_step = 0.01, l=0.
    std::cerr << "E_1s = " << E_tmp << std::endl << std::endl;
    // update 1s-orbital energy
    E_new[0] = E_tmp;
    // normalization of u_tmp
    normalize_u_tmp = 0.0;
    for(int i = 0; i < num_steps; i++){
      normalize_u_tmp += u_tmp[i]*u_tmp[i] * h;
    }
    // update 1s-orbital
    for(int i = 0; i < num_steps; i++){
      u_tmp[i] /= sqrt(normalize_u_tmp);
      u[0][i] = u[0][i]*(1-w) + u_tmp[i]*w;
    }
    
    // solve 2s-orbital
    E_left = E_new[0]+ 0.001; // start bisection from a energy which is slighly larger than 1s-orbital
    std::cerr << "SOLVE RADIAL EIGENEQUATION for 2s-orbital:" << std::endl;
    E_tmp = solve_radial_eigen_eq(E_left, 0.01, Z, n_4pi, rV_Coulomb, r, 0, u_tmp); // search_step = 0.01, l=0.
    std::cerr << "E_2s = " << E_tmp << std::endl << std::endl;
    // update 2s-orbital energy
    E_new[1] = E_tmp;
    // normalization of u_tmp
    normalize_u_tmp = 0.0;
    for(int i = 0; i < num_steps; i++){
      normalize_u_tmp += u_tmp[i]*u_tmp[i] * h;
    }
    // update 2s-orbital
    for(int i = 0; i < num_steps; i++){
      u_tmp[i] /= sqrt(normalize_u_tmp);
      u[1][i] = u[1][i]*(1-w) + u_tmp[i]*w;
    }
    
    // solve 2p-orbital
    E_left = -0.5*Z*Z/4.0;
    std::cerr << "SOLVE RADIAL EIGENEQUATION for 2p-orbital:" << std::endl;
    E_tmp = solve_radial_eigen_eq(E_left, 0.002, Z, n_4pi, rV_Coulomb, r, 1, u_tmp); // search_step = 0.01, l=1.
    std::cerr << "E_2p = " << E_tmp << std::endl << std::endl;
    // update 2p-orbital energy
    E_new[2] = E_tmp;
    // normalization of u_tmp
    normalize_u_tmp = 0.0;
    for(int i = 0; i < num_steps; i++){
      normalize_u_tmp += u_tmp[i]*u_tmp[i] * h;
    }
    // update 2p-orbital
    for(int i = 0; i < num_steps; i++){
      u_tmp[i] /= sqrt(normalize_u_tmp);
      u[2][i] = u[2][i]*(1-w) + u_tmp[i]*w;
    }
    
    // solve 3s-orbital
    E_left = E_new[1]+ 0.001; // start bisection from a energy which is slightly larger than 2s-orbital
    std::cerr << "SOLVE RADIAL EIGENEQUATION for 3s-orbital:" << std::endl;
    E_tmp = solve_radial_eigen_eq(E_left, 0.002, Z, n_4pi, rV_Coulomb, r, 0, u_tmp); // search_step = 0.005, l=0.
    std::cerr << "E_3s = " << E_tmp << std::endl << std::endl;

    // update 3s-orbital energy
    E_new[3] = E_tmp;
    // normalization of u_tmp
    normalize_u_tmp = 0.0;
    for(int i = 0; i < num_steps; i++){
      normalize_u_tmp += u_tmp[i]*u_tmp[i] * h;
    }
    // update 3s-orbital
    for(int i = 0; i < num_steps; i++){
      u_tmp[i] /= sqrt(normalize_u_tmp);
      u[3][i] = u[3][i]*(1-w) + u_tmp[i]*w;
    }
    
    // solve 3p-orbital
    E_left = E_new[2]+ 0.001; // start bisection from a energy which is slightly larger than 2p-orbital
    std::cerr << "SOLVE RADIAL EIGENEQUATION for 3p-orbital:" << std::endl;
    E_tmp = solve_radial_eigen_eq(E_left, 0.002, Z, n_4pi, rV_Coulomb, r, 1, u_tmp); // search_step = 0.002, l=1.
    std::cerr << "E_3p = " << E_tmp << std::endl << std::endl;
    // update 3p-orbital
    E_new[4] = E_tmp;
    // normalization of u_tmp
    normalize_u_tmp = 0.0;
    for(int i = 0; i < num_steps; i++){
      normalize_u_tmp += u_tmp[i]*u_tmp[i] * h;
    }
    // update 3p-orbital
    for(int i = 0; i < num_steps; i++){
      u_tmp[i] /= sqrt(normalize_u_tmp);
      u[4][i] = u[4][i]*(1-w) + u_tmp[i]*w;
    }
    
    // check convergence
    delta_E = 0.0;
    for(int i = 0; i < 5; i++){
      E_diff[i] = E_new[i] - E[i];
      delta_E = std::max(delta_E, std::abs(E_diff[i]));
    }
    
    // print current state
    std::cerr << "The radial equations have been solved:" << std::endl;
    std::cerr << "Obtained one particle energies: " << std::endl;
    std::cerr << "E_1s = " << E_new[0] << "(" << E_diff[0] << ")" << std::endl;
    std::cerr << "E_2s = " << E_new[1] << "(" << E_diff[1] << ")" << std::endl;
    std::cerr << "E_2p = " << E_new[2] << "(" << E_diff[2] << ")" << std::endl;
    std::cerr << "E_3s = " << E_new[3] << "(" << E_diff[3] << ")" << std::endl;
    std::cerr << "E_3p = " << E_new[4] << "(" << E_diff[4] << ")" << std::endl << std::endl;

    if(delta_E <= SCF_convergence){
      std::cerr << "The SCF loop has converged in " << SCF_loop << "-th loop:" << std::endl << std::endl;
      break;
    }
  }
  
  /* update state to the final answer */
  for(int i = 0; i < 5; i++){
    E[i] = E_new[i];
  }
  
  // normalization of u
  for(int i = 0; i < 5; i++){
    normalize_u[i] = 0.0;
    for(int j = 0; j < num_steps; j++){
      normalize_u[i] += u[i][j]*u[i][j] * h;
    }
    for(int j = 0; j < num_steps; j++){
      u[i][j] /= sqrt(normalize_u[i]);
    }
  }
  
  // calculate n_4pi
  for(int i = 0; i < num_steps; i++){
    n_4pi[i] = 0.0;
  }
  for(int i = 0; i < 5; i++){
    for(int j = 1; j < num_steps; j++){
      n_4pi[j] += N[i] * pow(u[i][j]/r[j], 2);
    }
  }n_4pi[0] = n_4pi[1]; // this is temporary : we do this because we cannot calculate u/r for r = 0. There is no problem for doing this because n_4pi[0] is used only for plotting and not used in solving radial differential equations
  
  // calculate rV_Coulomb
  rV_Coulomb[0] = 0.0; // these initial conditions are temporary.
  rV_Coulomb[1] = h;
  for(int i = 2; i < num_steps; i++){ // solve d^2/dr^2 rV_Coulomb(r) = - r*4pi*n(r)
    rV_Coulomb[i] = 2.0*rV_Coulomb[i-1] - rV_Coulomb[i-2] - h*h * r[i-1] * n_4pi[i-1];
  }
  for(int i = 1; i < num_steps; i++){ // determine uncertainty term
    rV_Coulomb[i] += (Z - rV_Coulomb[num_steps-1]) /((double)num_steps-1.0) * i;
  } // we impose rV_Coulomb(r_max) = Z : there is an correction term of exp(-r)
  
  /* calculate ground state energy */
  E_gnd = 0.0;
  // the sum of one particle energies
  for(int i = 0; i < 5; i++){
    E_gnd += N[i] * E[i];
  }
  // calculate exchange-correlation energy and correction energy
  E_xc = E_correction = 0.0;
  for(int i = 1; i < num_steps; i++){
    r_s = pow(3.0/(std::max(n_4pi[i], n_4pi_minimum)), 1.0/3.0);
    E_xc += n_4pi[i] * r[i]*r[i] * e_xc(r_s) * h;
    E_correction += n_4pi[i] * r[i]*r[i] * (e_xc(r_s) - V_xc(r_s) - 0.5 * rV_Coulomb[i]/r[i]) * h;
  }
  E_gnd += E_correction;
  
  // print result
  print_element_symbol(Z); // print the atom
  std::cerr << "The total number of protons/electrons: " << std::endl;
  std::cerr << "Z = " << Z << std::endl << std::endl;
  
  set_N(Z, N); // the definition is written in "system_settings.hpp"
  std::cerr << "The numbers of electrons in each orbitals:" << std::endl;
  std::cerr << "N_1s : " << N[0] << std::endl;
  std::cerr << "N_2s : " << N[1] << std::endl;
  std::cerr << "N_2p : " << N[2] << std::endl;
  std::cerr << "N_3s : " << N[3] << std::endl;
  std::cerr << "N_3p : " << N[4] << std::endl << std::endl;
  
  // Kohn-Sham energy eigenvalues
  std::cerr << "Kohn-Sham energy eigenvalues: " << std::endl;
  std::cerr << "E_1s = " << E_new[0] << std::endl;
  std::cerr << "E_2s = " << E_new[1] << std::endl;
  std::cerr << "E_2p = " << E_new[2] << std::endl;
  std::cerr << "E_3s = " << E_new[3] << std::endl;
  std::cerr << "E_3p = " << E_new[4] << std::endl << std::endl;
  
  std::cerr << "E_xc  = " << E_xc << std::endl;
  std::cerr << "E_gnd = " << E_gnd << std::endl << std::endl;
  
  // print obtained arrays
  for(int i = 0; i < num_steps; i++){
    std::cout << std::scientific << r[i] << " " << u[0][i] << " " << u[1][i] << " " << u[2][i] << " " << u[3][i] << " " << u[4][i] << " " << n_4pi[i] << " " << rV_Coulomb[i] << std::endl;
  }
  
  return 0;
}
