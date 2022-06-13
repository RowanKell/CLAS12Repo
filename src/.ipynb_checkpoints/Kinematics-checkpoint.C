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

double Kinematics::y(double E1, double E2){ //E1 and E2 here are the initial energy and final energy of the lepton (electron) see page 85 of probing nucleon spin doc
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

double Kinematics::phi_h(TLorentzVector Q, TLorentzVector L, TLorentzVector p1, TLorentzVector p2){
  TLorentzVector ph = p1 + p2;
  TLorentzVector r = 0.5*(p1-p2);

  TVector3 q(Q.Px(), Q.Py(), Q.Pz());
  TVector3 l(L.Px(), L.Py(), L.Pz());
  TVector3 Ph(ph.Px(), ph.Py(), ph.Pz());

  TVector3 qcrossl = q.Cross(l);
  TVector3 qcrossPh = q.Cross(Ph);

  double factor1 = (qcrossl*Ph)/abs(qcrossl*Ph);
double factor2 = (qcrossl*qcrossPh)/qcrossl.Mag()/qcrossPh.Mag();

  return factor1*acos(factor2);
}

double Kinematics::phi_R(TLorentzVector Q, TLorentzVector L, TLorentzVector p1, TLorentzVector p2){
  TLorentzVector ph = p1 + p2;
  TLorentzVector r = 0.5*(p1-p2);

  TVector3 q(Q.Px(), Q.Py(), Q.Pz());
  TVector3 l(L.Px(), L.Py(), L.Pz());
  TVector3 R(r.Px(), r.Py(), r.Pz());

  TVector3 Rperp = R-(q*R)/(q*q)*q;
  TVector3 qcrossl = q.Cross(l);
  TVector3 qcrossRperp = q.Cross(Rperp);

  double factor1 = (qcrossl*Rperp)/abs(qcrossl*Rperp);
  double factor2 = (qcrossl*qcrossRperp)/qcrossl.Mag()/qcrossRperp.Mag();

  return factor1*acos(factor2);
}

double Kinematics::com_th(TLorentzVector P1, TLorentzVector P2){
  TLorentzVector Ptotal = P1+P2;
  TVector3 comBOOST = Ptotal.BoostVector();
  Ptotal.Boost(-comBOOST);
  P1.Boost(-comBOOST);
  return P1.Angle(comBOOST);
}
