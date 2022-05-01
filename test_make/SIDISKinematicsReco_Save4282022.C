// ----------------------------
//
// SIDISKinematicsReco.C()
// 
// Written by Gregory Matousek
//
// Purpose
// ----------------
// Read in MC or data hipo file and extract events
// with scattered electron, pi+/pi-, and at least 2 photons
// 
// The extracted data is stored within a ROOT Tree for later cuts
//
// Also, basic SIDIS cuts are applied to the data
// 
// Changelog
// ----------------
// 4/27/2022: Program was created for the first time
// 4/28/2022: Program was renamed
//            Adding support for Makefile
// 
// ----------------------------

#include "./src/Settings.C"
#include "./src/Constants.h"

//int SIDISKinematicsReco(int argc, char *argv[])
int SIDISKinematicsReco(char *argv)
{
  int argc=1;
  // Read cut settings
  // -------------------------
  Settings settings;
  if(argc<1)
    {
      std::cout << " Missing arguments. Please provide runcard name " << std::endl;
      return -1;
    }
  else
    {
      if(!settings.readCard(argv)){
	return -1;
      }
      settings.loadSettings();
    }


  // Grab Hipofile name
  // -------------------------  
  string hipofile = settings.Filename();
  
  // Create HipoChain
  // -------------------------  
  HipoChain chain;
  chain.Add(hipofile);

  // Create C12reader
  // -------------------------  
  auto config_c12=chain.GetC12Reader();

  // Get specific final states
  // i) At least 1 electron
  // ii) One charged pion
  // iii) At least 2 photons
  // -------------------------  
  config_c12->addAtLeastPid(11,1);

  if(settings.piPID()==0){
    std::cout << " Missing pion charge in CutCard. Must be included (e.g. picharge        1) " << std::endl;
    return -1;
  }
  else{
    config_c12->addExactPid(settings.piPID(),1);
  }
  config_c12->addAtLeastPid(22,2);
  
  auto& c12=chain.C12ref();
  
  // Parse through HIPO file
  // -------------------------  
  while(chain.Next()==true){
    if(c12->getDetParticles().empty())
      continue;
    
    //    auto parts=c12->getDetParticles();
    std::cout << c12->getNPid(11) << std::endl;
  }
  

  return 0;
}
