#pragma once
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)
// If anybody can figure out how to have root automically search in the include
// directory for this file, it would be GREATLY appreciated! 
#include "/work/clas12/users/gmat/CLAS12Analysis/include/clas12ana/SIDISKinematicsReco.h"
#include "/work/clas12/users/gmat/CLAS12Analysis/include/clas12ana/Settings.h"
R__LOAD_LIBRARY(libclas12ana.so)
#endif

using namespace std;

int pipi0_process(
			const char * outputFile = "/work/clas12/users/gmat/CLAS12Analysis/data/raw/test.root"
		  )
{
  //---------------
  // Load libraries
  //---------------
  gSystem->Load("libclas12ana.so"); 
 
  //-----------------------------------
  // Create SIDISKinematicsReco Object 
  //-----------------------------------
  SIDISKinematicsReco *ana = new SIDISKinematicsReco(outputFile);
  //-----------------------------------
  // Create Settings for Study
  //-----------------------------------
  Settings settings;
  settings.setQ2range(1,100);
  settings.setWrange(2,100);
  settings.setyrange(0,1);
  settings.addFinalState(11,1,true);  // Exactly 1 scattered electron
  settings.addFinalState(211,1,true); // Exactly 1 pi+
  settings.addFinalState(22,2,false); // 2 or more gammas
  settings.setdoMC(true);             // Analyze MC::Lund
  settings.setdoReco(true);           // Analyze REC::Particle
  settings.setconnectMC2Reco(true);   // Connect pindex of REC::Particle to pindex of MC::Lund
  
  std::string dirname = "/w/hallb-scshelf2102/clas12/users/gmat/CLAS12Analysis/data/raw/sample/";
  settings.addHipoFile(dirname+std::string("Out_DIS_pass1_915_920.hipo_skim23.hipo"));
  //-----------------------------------
  // Import Settings into Processing Framework
  //-----------------------------------
  ana->ImportSettings(settings);
  //-----------------------------------
  // Initialize the Analysis Module
  //-----------------------------------
  ana->Init();
  //-----------------------------------
  // Process all events
  //-----------------------------------
  ana->process_events();
  //-----------------------------------
  // End processing
  //-----------------------------------
  ana->End();
  return 0;
}
