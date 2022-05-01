#ifndef Kinematics_h
#define Kinematics_h

//#include "Constants.h"


class Kinematics {
 public:
  static double Q2(double E1, double E2, double cth);
  static double Pt(double Px, double Py);
  static double cth(double Px, double Py, double Pz);
  static double y(double E1, double E2);
  static double nu(double E1, double E2);
  static double W(double Q2, double mT, double nu);
};
#endif
