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
		  const char * hipoFile = "/cache/clas12/rg-a/production/montecarlo/clasdis/fall2018/torus-1/v1/bkg45nA_10604MeV/45nA_job_3301_3.hipo",
		  //const char * hipoFile = "/w/hallb-scshelf2102/clas12/users/gmat/CLAS12Analysis/data/raw/sample/Out_DIS_pass1_915_920.hipo_skim23.hipo",
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
  ana->set_verbosity(1);
  ana->set_printFrequency(100000);
  //-----------------------------------
  // Create Settings for Study
  //-----------------------------------
  Settings settings;
  
  settings.setElectronBeamEnergy(10.6);
  
  settings.setQ2range(1,100);
  settings.setWrange(2,100);
  settings.setyrange(0,0.8);
  settings.addFinalState(11,1,true);  // Exactly 1 electron
  settings.addFinalState(211,1,true); // Exactly 1 pi+
  settings.addFinalState(22,2,false); // 2 or more gammas
  settings.setignoreOtherRecoParticles(true); // Do not save to tree particle info of uninterested PIDs

  settings.addPIDforEmin(22,0.6);     // Gammas must have minimum energy of 0.6 GeV
  settings.addPIDforPmin(211,1.25);   // Pi+ must have minimum momentum of 1.25 GeV
  settings.addPIDforVzrange(11,-8,3); // e- must have vertex 'z' between [-13,12] cm
  settings.addPIDforBetarange(22,0.9,1.1); // Beta range for photon
  settings.addPIDforChi2max(211,3);        // Pi+ must have abs(chi2pid) < 3 

  settings.setdoMC(true);             // Analyze MC::Lund
  settings.setdoReco(true);           // Analyze REC::Particle
  settings.setdoPostProcess(true);    // Apply further cuts
  settings.setPostProcessMethod("pipluspi0_MC"); // Perform pipluspi0 default processing
  settings.setconnectMC2Reco(true);   // Connect pindex of REC::Particle to pindex of MC::Lund
  settings.setEventRecoMethod(Settings::eventRecoMethod::useLargestPinFD);
  settings.addHipoFile(hipoFile);
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
  // Perform post-processing
  //-----------------------------------
  ana->PostProcessReco();
  //-----------------------------------
  // End processing
  //-----------------------------------
  ana->End();
  return 0;
}
