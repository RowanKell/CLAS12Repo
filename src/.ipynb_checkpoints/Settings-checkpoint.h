#ifndef Settings_h
#define Settings_h
#include <iostream>         
#include <string>         
#include <vector>  
#include <algorithm>
using namespace std;

class Settings{
 

 public:
  Settings();
  
  enum eventRecoMethod
  {
    useRecKinematicsBank = 1,
    useLargestPinFD = 2
  };

  bool doMC() const;
  bool doReco() const;
  bool doFiducialCuts() const;
  bool doPostProcess() const;
  bool connectMC2Reco() const;
  bool ignoreOtherRecoParticles() const;
  eventRecoMethod getEventRecoMethod() const;
  const char * postProcessMethod() const;
  double electronBeamEnergy() const;
  double Q2min() const;
  double Q2max() const;
  double Wmin() const;
  double Wmax() const;
  double ymin() const;
  double ymax() const;

  void setdoMC(bool);
  void setdoReco(bool);
  void setdoFiducialCuts(bool);
  void setdoPostProcess(bool);
  void setconnectMC2Reco(bool);
  void setignoreOtherRecoParticles(bool);
  void setEventRecoMethod(eventRecoMethod);
  void setPostProcessMethod(const char *);
  void setElectronBeamEnergy(double);
  void setQ2range(double, double);
  void setWrange(double, double);
  void setyrange(double, double);
  void setabschi2pidmax(double);
  void addFinalState(int, int, bool);
  void addPIDforEmin(int,double);
  void addPIDforPmin(int,double);
  void addPIDforVzrange(int,double, double);
  void addPIDforBetarange(int,double, double);
  void addPIDforChi2max(int,double);
  void addHipoFile(std::string);

  std::vector<int> getFinalStatePIDs();
  int getN_fromPID(int);
  bool isExact_fromPID(int);
  double getEmin_fromPID(int);
  double getPmin_fromPID(int);
  double getVzmin_fromPID(int);
  double getVzmax_fromPID(int);
  double getBetamin_fromPID(int);
  double getBetamax_fromPID(int);
  double getChi2max_fromPID(int);
  std::vector<std::string> hipoFileStrings();

  bool needsChi2PidCut(int);
 private:

  bool _doMC = false;
  bool _doReco = false;
  bool _doFiducialCuts = false;
  bool _doPostProcess = false;
  bool _connectMC2Reco = false;
  bool _ignoreOtherRecoParticles = false;
  eventRecoMethod _eventRecoMethod;
  const char *_postProcessMethod;
  double _electronBeamEnergy = 10.6;
  double _Q2min = -999;
  double _Q2max = 999;
  double _Wmin  = -999;
  double _Wmax  = 999;
  double _ymin  = 0;
  double _ymax  = 1;

  // Vectors for final state
  // _fPID --> {11, 211, 22} = {e-, pi+, gamma}
  // _fNpart --> {1, 1, 2}   = {1 e-, 1 pi+, 2 gammas}
  // _fExact --> {True, True, False}
  //    "True" = Exactly 'n' per event
  //    "False" = Strictly '>=n' per event
  // ------------------------
  std::vector<int> _fPID;
  std::vector<int> _fNpart;
  std::vector<bool> _fExact;
  
  // std::vector of Hipo filename strings
  std::vector<std::string> _hipoFileStrings;
  
  // std::vector of Emin
  std::vector<int> _EminPID;
  std::vector<float> _Emin;

  // std::vector of Pmin
  std::vector<int> _PminPID;
  std::vector<float> _Pmin;

  // std::vector of Vzmin and Vzmax
  std::vector<int> _VzPID;
  std::vector<float> _Vzmin;
  std::vector<float> _Vzmax;

  // std::vector of Betamin and Betamax
  std::vector<int> _BetaPID;
  std::vector<float> _Betamin;
  std::vector<float> _Betamax;

  // std::vector of chi2max
  std::vector<int> _chi2PID;
  std::vector<float> _chi2max;

};
#endif
