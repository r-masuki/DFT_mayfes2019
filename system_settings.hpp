#ifndef system_settings_hpp
#define system_settings_hpp

#include <vector>

/* scan Z */
int scan_Z(char *arg, int *Z){
  // scanf Z
  sscanf(arg, "%d", Z);
  if((*Z) >= 1 && (*Z) <= 18){ // if the Z is in the range
    return 0;
  }
  else return 1; // error : Z is out of range
}

/* set N */
// calculate the numbers of electrons in each orbitals from total electron number Z.
int set_N(int Z, int *N){
  // set N_1s
  if(Z == 1) N[0] = 1;
  else N[0] = 2;
  
  // set N_2s
  if(Z <= 2) N[1] = 0;
  else if(Z == 3) N[1] = 1;
  else N[1] = 2;
  
  // set N_2p
  if(Z <= 4) N[2] = 0;
  else if(Z >= 5 && Z <= 10) N[2] = Z-4;
  else N[2] = 6;
  
  // set N_3s
  if(Z <= 10) N[3] = 0;
  else if(Z == 11) N[3] = 1;
  else N[3] = 2;
  
  // set N_3p
  if(Z <= 12) N[4] = 0;
  else if(Z >= 13 && Z <= 18) N[4] = Z - 12;
  else N[4] = 0;
  
  return 0;
}

int print_element_symbol(int Z){
  std::cerr << "The element symbol: ";
  if(Z == 1) std::cerr << "H";
  else if(Z == 2) std::cerr << "He";
  else if(Z == 3) std::cerr << "Li";
  else if(Z == 4) std::cerr << "Be";
  else if(Z == 5) std::cerr << "B";
  else if(Z == 6) std::cerr << "C";
  else if(Z == 7) std::cerr << "N";
  else if(Z == 8) std::cerr << "O";
  else if(Z == 9) std::cerr << "F";
  else if(Z == 10) std::cerr << "Ne";
  else if(Z == 11) std::cerr << "Na";
  else if(Z == 12) std::cerr << "Mg";
  else if(Z == 13) std::cerr << "Al";
  else if(Z == 14) std::cerr << "Si";
  else if(Z == 15) std::cerr << "P";
  else if(Z == 16) std::cerr << "S";
  else if(Z == 17) std::cerr << "Cl";
  else if(Z == 18) std::cerr << "Ar";
  
  std::cerr << std::endl << std::endl;
  
  return 0;
}

#endif /* system_settings_hpp */
