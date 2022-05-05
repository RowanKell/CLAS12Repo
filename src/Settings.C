#include "Settings.h"
#include <fstream>
#include <sstream>
#include <cstdlib>
Settings::Settings()
{
  
}

bool Settings::doMC() const {return _doMC;}
bool Settings::doReco() const {return _doReco;}
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

std::vector<int> Settings::getFinalStatePIDs(){ return _fPID; }

int Settings::getN_fromPID(int pid){
  for(unsigned int idx = 0 ; idx < _fPID.size() ; idx++)
    {
      if(_fPID.at(idx)==pid)
	return _fNpart.at(idx);
    }
  
  return -1;
}

bool Settings::isExact_fromPID(int pid){
  for(unsigned int idx = 0 ; idx < _fPID.size() ; idx++)
    {
      if(_fPID.at(idx)==pid)
	return _fExact.at(idx);
    }
  
  return false;
}

std::vector<std::string> Settings::hipoFileStrings() {return _hipoFileStrings;}
