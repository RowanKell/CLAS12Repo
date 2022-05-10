#ifndef PostProcess_h
#define PostProcess_h

#include <TTree.h>
#include <TLorentzVector.h>
#include <iostream>
#include <fstream>
#include "Constants.h"
#include <math.h>
using namespace std;
class PostProcess{

 public:
  PostProcess();

  enum PROCESS_ID{
    
    pipluspi0 = 0,
    pipluspi0_MC = 1,
    piminuspi0 = 10,
    piminuspi0_MC = 11

  };
  int Init(TTree *, double);
  
  int Process(TTree *);

  int setPostProcessMethod(const char * method){
  
    if(strcmp(method, "pipluspi0")==0){
      _process_id = PROCESS_ID::pipluspi0;
      _piPID=211;
    }

    else if(strcmp(method, "piminuspi0")==0){
      _process_id = PROCESS_ID::piminuspi0;
      _piPID=-211;
    }

    else if(strcmp(method, "pipluspi0_MC")==0){
      _process_id = PROCESS_ID::pipluspi0_MC;
      _piPID=211;
    }

    else if(strcmp(method, "piminuspi0_MC")==0){
      _process_id = PROCESS_ID::piminuspi0_MC;
      _piPID=-211;
    }

    else{
      std::cout << "ERROR in PostProcess::setPostProcessMethod -- Unrecognized method " << method << " . Breaking . . . " << std::endl;
      return -1;
    }
    return 0;
  }

  int pipi0(TTree *);
  int pipi0_MC(TTree *);
  

  // private:
  Float_t         Q2;
  Float_t         W;
  Float_t         nParticles;
  Float_t         nPhotons;
  Float_t         nu;
  Float_t         x;
  Float_t         y;
  vector<float>   *pid=0;
  vector<float>   *px=0;
  vector<float>   *py=0;
  vector<float>   *pz=0;
  vector<float>   *pt=0;
  vector<float>   *p=0;
  vector<float>   *E=0;
  vector<float>   *theta=0;
  vector<float>   *eta=0;
  vector<float>   *phi=0;
  vector<float>   *vz=0;
  vector<float>   *pindex=0;
  vector<float>   *beta=0;
  vector<float>   *chi2=0;
  vector<float>   *parentID=0;
  vector<float>   *parentPID=0; // Declaration of leaf types                                                                                                                    
  // List of branches      
  TBranch        *b_Q2;        
  TBranch        *b_W;
  TBranch        *b_nParticles; 
  TBranch        *b_nPhotons; 
  TBranch        *b_nu; 
  TBranch        *b_x;
  TBranch        *b_y;
  TBranch        *b_pid;
  TBranch        *b_px;
  TBranch        *b_py;
  TBranch        *b_pz;
  TBranch        *b_pt;
  TBranch        *b_p; 
  TBranch        *b_E;                                                                                                                                                    
  TBranch        *b_theta;
  TBranch        *b_eta;
  TBranch        *b_phi; 
  TBranch        *b_vz;
  TBranch        *b_pindex;
  TBranch        *b_beta;   
  TBranch        *b_chi2;
  TBranch        *b_parentID;
  TBranch        *b_parentPID;

  TTree *_tree_Reco;
  Long64_t _nentries;
  
  int _piPID;

  PROCESS_ID _process_id;

  double _electron_beam_energy;
};
#endif
