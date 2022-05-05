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
  
  bool doMC() const;
  bool doReco() const;
  double Q2min() const;
  double Q2max() const;
  double Wmin() const;
  double Wmax() const;
  double ymin() const;
  double ymax() const;

  void setdoMC(bool);
  void setdoReco(bool);
  void setQ2range(double, double);
  void setWrange(double, double);
  void setyrange(double, double);
  void addFinalState(int, int, bool);
  void addHipoFile(std::string);

  std::vector<int> getFinalStatePIDs();
  int getN_fromPID(int);
  bool isExact_fromPID(int);
  
  std::vector<std::string> hipoFileStrings();

 private:

  bool _doMC = false;
  bool _doReco = false;

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
  
};
#endif
