#pragma once
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)
// If anybody can figure out how to have root automically search in the include
// directory for this file, it would be GREATLY appreciated! 
#include "/work/clas12/users/gmat/CLAS12Analysis/include/clas12ana/PostProcess.h"
R__LOAD_LIBRARY(libclas12ana.so)
#endif

using namespace std;

int pipi0_postprocess_only(
		  const char * inputFile = "/work/clas12/users/gmat/CLAS12Analysis/data/raw/test.root"
		  )
{
  //---------------
  // Load libraries
  //---------------
  gSystem->Load("libclas12ana.so"); 
 
  
  //-----------------------
  // Load TFile and TTree's
  //-----------------------
  TFile *fIn = new TFile(inputFile,"UPDATE");
  TTree *tReco = (TTree*)fIn->Get("tree_reco");
  if(fIn->Get("tree_postprocess"))
    gDirectory->Delete("tree_postprocess");  
  TTree *tPostProcess = new TTree("tree_postprocess","Tree created after all hipo events are filtered");    
  //-----------------------------------
  // Create PostProcess Object 
  //-----------------------------------
  PostProcess *pp = new PostProcess();
  pp->Init(tReco,10.6);
  pp->setPostProcessMethod("pipluspi0");
  
  //-----------------------------------
  // Post process events and save results in tPostProcess 
  //-----------------------------------
  pp->Process(tPostProcess);
  
  //-----------------------------------
  // Write tPostProcess to .root file
  //-----------------------------------
  tPostProcess->Write();

  //-----------------------------------
  // Close .root file
  //-----------------------------------
  fIn->Close();
  return 0;
}
