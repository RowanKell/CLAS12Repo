#ifndef Settings_h
#define Settings_h
#include <iostream>         
#include <string>         
#include <vector>  
#include <map>
using namespace std;

class Settings{
 

 public:
  Settings();
  bool readCard(const char*);
  bool loadSettings();
  
  string Filename() const;
  int IsMC() const;
  double Q2min() const;
  double Wmin() const;
  double ymax() const;
  int picharge() const;
  int piPID() const;

 private:
  vector<string> _lines;
  map<string,string> _settings;
  string _cardname;

  string _Filename = "";
  bool _IsMC;
  double _Q2min = -999;
  double _Wmin  = -999;
  double _ymax  = 1;
  int    _picharge = 0;
};
#endif
