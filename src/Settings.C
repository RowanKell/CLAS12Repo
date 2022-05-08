#include "Settings.h"
#include <fstream>
#include <sstream>
#include <cstdlib>
Settings::Settings()
{
  
}

bool Settings::doMC() const {return _doMC;}
bool Settings::doReco() const {return _doReco;}
bool Settings::connectMC2Reco() const {return _connectMC2Reco;}
bool Settings::ignoreOtherRecoParticles() const {return _ignoreOtherRecoParticles;}
Settings::eventRecoMethod Settings::getEventRecoMethod() const {return _eventRecoMethod;}
double Settings::electronBeamEnergy() const {return _electronBeamEnergy;}
double Settings::Q2min() const {return _Q2min;}
double Settings::Q2max() const {return _Q2max;}
double Settings::Wmin() const {return _Wmin;}
double Settings::Wmax() const {return _Wmax;}
double Settings::ymin() const {return _ymin;}
double Settings::ymax() const {return _ymax;}

void Settings::setdoMC(bool b){
  _doMC = b;
  return;
}

void Settings::setdoReco(bool b){
  _doReco = b;
  return;
}

void Settings::setconnectMC2Reco(bool b){
  _connectMC2Reco = b;
  return;
}

void Settings::setignoreOtherRecoParticles(bool b){
  _ignoreOtherRecoParticles = b;
  return;
}

void Settings::setEventRecoMethod(eventRecoMethod eRM){
  _eventRecoMethod = eRM;
  return;
}

void Settings::setElectronBeamEnergy(double E){
  _electronBeamEnergy = E;
  return;
}

void Settings::setQ2range(double Q2min, double Q2max) {
  _Q2min = Q2min;
  _Q2max = Q2max;
  return;
}

void Settings::setWrange(double Wmin, double Wmax) {
  _Wmin = Wmin;
  _Wmax = Wmax;
  return;
}

void Settings::setyrange(double ymin, double ymax) {
  _ymin = ymin;
  _ymax = ymax;
  return;
}


void Settings::addFinalState(int pid, int n, bool exact=false) {
  // Make sure pid isn't a duplicate
  if(std::count(_fPID.begin(), _fPID.end(), pid)){
    std::cout << "ERROR in Settings::addFinalState -- " << pid << " is a repeat! Skipping..." << std::endl;
  }
  else if(n<0)
    {
      std::cout << "ERROR in Settings::addFinalState -- " << "n cannot be negative! Skipping..." << std::endl;
    }
  else
    {
      _fPID.push_back(pid);
      _fNpart.push_back(n);
      _fExact.push_back(exact);
    }
  return;
}

void Settings::addHipoFile(std::string filename){
  _hipoFileStrings.push_back(filename);
}

void Settings::addPIDforEmin(int pid, double Emin){
  _EminPID.push_back(pid);
  _Emin.push_back(Emin);
}

void Settings::addPIDforPmin(int pid, double Pmin){
  _PminPID.push_back(pid);
  _Pmin.push_back(Pmin);
}

void Settings::addPIDforVzrange(int pid, double Vzmin, double Vzmax) {
  _VzPID.push_back(pid);
  _Vzmin.push_back(Vzmin);
  _Vzmax.push_back(Vzmax);
  return;
}

void Settings::addPIDforChi2max(int pid, double chi2max){
  _chi2PID.push_back(pid);
  _chi2max.push_back(chi2max);
  return;
}

std::vector<int> Settings::getFinalStatePIDs(){ return _fPID; }

int Settings::getN_fromPID(int pid){
  for(unsigned int idx = 0 ; idx < _fPID.size() ; idx++)
    {
      if(_fPID.at(idx)==pid)
	return _fNpart.at(idx);
    }
  
  return 0;
}

bool Settings::isExact_fromPID(int pid){
  for(unsigned int idx = 0 ; idx < _fPID.size() ; idx++)
    {
      if(_fPID.at(idx)==pid)
	return _fExact.at(idx);
    }
  return false;
}

double Settings::getEmin_fromPID(int pid){
 for(unsigned int idx = 0 ; idx < _EminPID.size() ; idx++)
   {
     if(_EminPID.at(idx)==pid)
       return _Emin.at(idx);
   }
 return 0;
} 

double Settings::getPmin_fromPID(int pid){
 for(unsigned int idx = 0 ; idx < _PminPID.size() ; idx++)
   {
     if(_PminPID.at(idx)==pid)
       return _Pmin.at(idx);
   }
 return 0;
} 


double Settings::getVzmin_fromPID(int pid){
 for(unsigned int idx = 0 ; idx < _VzPID.size() ; idx++)
   {
     if(_VzPID.at(idx)==pid)
       return _Vzmin.at(idx);
   }
 return -99999;
} 

double Settings::getVzmax_fromPID(int pid){
 for(unsigned int idx = 0 ; idx < _VzPID.size() ; idx++)
   {
     if(_VzPID.at(idx)==pid)
       return _Vzmax.at(idx);
   }
 return 99999;
} 

double Settings::getChi2max_fromPID(int pid){
 for(unsigned int idx = 0 ; idx < _chi2PID.size() ; idx++)
   {
     if(_chi2PID.at(idx)==pid)
       return _chi2max.at(idx);
   }
 return 99999;
} 



std::vector<std::string> Settings::hipoFileStrings() {return _hipoFileStrings;}

bool Settings::needsChi2PidCut(int pid){
  for(unsigned int idx = 0 ; idx < _chi2PID.size() ; idx++)
    {
      if(_chi2PID.at(idx)==pid)
	return true;
    }
  return false;
}

