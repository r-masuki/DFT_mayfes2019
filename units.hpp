#ifndef units_hpp
#define units_hpp

namespace units{
  // length
  static const double um = 1.0e-6; //[m]
  static double const nm = 1.0e-9; //[m]
  static const double Angstroem = 1.0e-10; //[m]
  static const double pm = 1.0e-12; //[m]
  
  // mass
  static const double u_atomic = 1.660539040e-27; //[kg]
  
  
  //energy [J]
  static const double eV = 1.6021766208e-19; //[J]
  static const double meV = 1.6021766208e-22; //[J]
  
  // Hartree Atomic Units
  namespace Hartree{
    // length
    static const double a_Bohr =  0.52917721067e-10; //[m]
    // mass
    static const double m_e =9.10938356e-31; //[kg]
    // charge
    static const double e_charge = 1.6021766208e-19; //[C]
    // energy
    static const double energy_unit = 4.359744650e-18; //[J]
  }
  
}

#endif /* units_hpp */
