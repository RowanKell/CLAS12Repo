#ifndef FiducialCuts_h
#define FiducialCuts_h

/* Src Includes */
#include "Constants.h"
/* CLAS12 Includes */
#include "HipoChain.h"
#include "clas12reader.h"

using namespace std;

class FiducialCuts{
 public:
  FiducialCuts();
  FiducialCuts(const std::unique_ptr<clas12::clas12reader>&);
  
  bool FidCutParticle(const std::unique_ptr<clas12::clas12reader>&, int , int, float);

 protected:
  int _idx_RECCal;
  int _ipindex;
  int _ilv;
  int _ilw;

};
#endif
