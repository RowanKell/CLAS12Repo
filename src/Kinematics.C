#include "Kinematics.h"
#include <math.h>
using namespace std;

double Kinematics::Q2(double E1, double E2, double cth){
  return 2.0 * E1 * E2 * (1.0 - cth);
}

double Kinematics::x(double Q2, double s, double y){
  return Q2/s/y;
}

double Kinematics::Px(double P, double th, double phi){
  return P * sin(th) * cos(phi);
}

double Kinematics::Py(double P, double th, double phi){
  return P * sin(th) * sin(phi);
}

double Kinematics::Pz(double P, double th, double phi){
  return P * cos(th);
}

double Kinematics::Pt(double Px, double Py){
  return pow(Px*Px + Py*Py,0.5);
}

double Kinematics::P(double Px, double Py, double Pz){
  return pow(Px*Px + Py*Py + Pz*Pz,0.5);
}

double Kinematics::E(double M, double P){
  return pow(M*M+P*P,0.5);
}

double Kinematics::cth(double Px, double Py, double Pz){
  double Pt = Kinematics::Pt(Px,Py);
  return Pz / (pow(Pz*Pz+Pt*Pt,0.5));
}

double Kinematics::y(double E1, double E2){
  return (E1-E2)/E1;
}

double Kinematics::nu(double E1, double E2){
  return E1-E2;
}

double Kinematics::W(double Q2, double mT, double nu){
  return sqrt(-Q2 + pow(mT,2) + 2 * mT * nu);
} 

double Kinematics::th(double Pt, double Pz){
  return abs(atan(Pt/Pz));
}

double Kinematics::eta(double th){
  return -1.0 * log ( tan ( th / 2.0 ) );
}

double Kinematics::phi(double Px, double Py){
  return atan2(Py,Px);
}

