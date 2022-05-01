#include "Kinematics.h"
#include <math.h>
using namespace std;

double Kinematics::Q2(double E1, double E2, double cth){
  return 2.0 * E1 * E2 * (1.0 - cth);
}

double Kinematics::Pt(double Px, double Py){
  return pow(Px*Px + Py*Py,0.5);
}

double Kinematics::cth(double Px, double Py, double Pz){
  double Pt = Kinematics::Pt(Px,Py);
  return cos(Pz / (pow(Pz*Pz+Pt*Pt,0.5)));
}

double Kinematics::y(double E1, double E2){
  return (E1-E2)/E1;
}

double Kinematics::nu(double E1, double E2){
  return E1-E2;
}

double Kinematics::W(double Q2, double mT, double nu){
  return -Q2 + pow(mT,2) + 2 * mT * nu;
} 
