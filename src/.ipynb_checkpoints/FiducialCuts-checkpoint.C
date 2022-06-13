#include "FiducialCuts.h"

using namespace std;

FiducialCuts::FiducialCuts(){}

FiducialCuts::FiducialCuts(const std::unique_ptr<clas12::clas12reader>& _c12){
  _idx_RECCal = _c12->addBank("REC::Calorimeter");
  _ipindex = _c12->getBankOrder(_idx_RECCal,"pindex");
  _ilv = _c12->getBankOrder(_idx_RECCal,"lv");
  _ilw = _c12->getBankOrder(_idx_RECCal,"lw");
}

bool FiducialCuts::FidCutParticle(const std::unique_ptr<clas12::clas12reader>& _c12, int pid, int pindex, float theta){
  
  // Loop over entries in the Calorimeter bank
  for(auto i = 0 ; i < _c12->getBank(_idx_RECCal)->getRows() ; i++){
    // Continue loop if the pindex in the calo bank does not match 
    if(_c12->getBank(_idx_RECCal)->getInt(_ipindex,i)!=pindex)
      continue;
    
    // Perform standard fiducial cuts
    if(theta<5*PI/180.0 || theta>35*PI/180.0)
      return false;

    float lv = _c12->getBank(_idx_RECCal)->getFloat(_ilv,i);
    float lw = _c12->getBank(_idx_RECCal)->getFloat(_ilw,i);

    // Perform electron (pid==11) fiducial cuts
    if(pid==11){
      if(lv<9 || lw<9)
	return false;
    }
   
    // Perform photon (pid==22) fiducial cuts
    if(pid==22){
      if(lv<14 || lw <14)
	return false;
    }
  }
  
  return true;
}
