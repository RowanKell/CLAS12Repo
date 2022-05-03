#include "Settings.h"
#include <fstream>
#include <sstream>
#include <cstdlib>
Settings::Settings()
{
  
}

bool Settings::readCard(const char *file)
{
  if(!file)
    return false;
  _cardname = string(file);
  
  
  // Open the file
  // ------------------
  ifstream ifs(file);
  if(!ifs){
    cout << "Cannot open filename --> " << file << endl;
    return false;
  }
  
  // Read the file
  // ------------------
  while (ifs.good() && !ifs.eof()){
    string line;
    getline(ifs,line);
    
    // If the line is empty and at eof, break
    if(ifs.eof() && line.empty()) break;
    
    // Skip empty lines
    if(line.find_first_not_of(" \n\t\v\b\r\f\a") == string::npos) continue;
    
    // Skip lines which do not begin with a letter or digit
    int firstChar = line.find_first_not_of(" \n\t\v\b\r\f\a");
    if(!isalnum(line[firstChar])) continue;

    // Save the useful line
    _lines.push_back(line);
  }
  
  // Close the file
  // ------------------
  ifs.close();

  // Parse through the lines
  // ------------------------
  for(string l: _lines)
    {
      string key,val;
      istringstream ss(l);
      ss>>key;
      ss>>val;
      _settings[key] = val;
    }
  return true;
}  

bool Settings::loadSettings()
{
  for( map<string,string>::iterator it = _settings.begin(); it!=_settings.end(); it++)
    {
      string key = (*it).first;
      const char *val = ((*it).second).c_str();

      if(key.compare("Filename") == 0){
	_Filename = val;
      }else if(key.compare("IsMC") == 0){
	_IsMC = atoi(val);
      }else if(key.compare("Q2min") == 0){
	_Q2min = atof(val);
      }else if(key.compare("Wmin") == 0){
	_Wmin  = atof(val);
      }else if(key.compare("ymax") == 0){
	_ymax  = atof(val);
      }else if(key.compare("picharge") == 0){
	string temp = string(val);
	if(temp.compare("plus") == 0 || temp.compare("+") == 0 || temp.compare("1") == 0)
	  _picharge = 1;
	else if(temp.compare("minus") == 0 || temp.compare("-") == 0 || temp.compare("-1") == 0)
	  _picharge = -1;
	else
	  {
	    cout << "Pion charge error, see Settings.C for options" << endl;
	    return false;
	  }
      }
    }
  return true;
}
string Settings::Filename() const{return _Filename;}
int Settings::IsMC() const {return _IsMC;}
double Settings::Q2min() const {return _Q2min;}
double Settings::Wmin() const {return _Wmin;}
double Settings::ymax() const {return _ymax;}
int Settings::picharge() const {return _picharge;}
int Settings::piPID() const {return _picharge*211;}
