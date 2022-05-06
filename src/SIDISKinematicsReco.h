#ifndef SIDISKinematicsReco_h
#define SIDISKinematicsReco_h

/* STL Includes */
#include <math.h>
#include <map>
#include <utility>
/* Src Includes */
#include "Settings.h"
#include "Constants.h"
#include "Kinematics.h"
#include "SIDISParticle.h"
#include "SIDISParticlev1.h"
/* ROOT includes */
#include <TFile.h>
#include <TTree.h>
/* CLAS12 includes */
#include "HipoChain.h"
#include "clas12reader.h"
#include "region_particle.h"

class SIDISParticle;

typedef std::map<int, SIDISParticle*> type_map_part;

class SIDISKinematicsReco{
 
 public:
  
  SIDISKinematicsReco(std::string);

  int Init();

  int InitHipo();

  void ImportSettings(Settings theSettings) {
    _settings = theSettings;
  }

  int process_events();

  int End();
  
  void set_beam(double E){
    _electron_beam_energy = E;
  }
 
 private:
 
  int _ievent;
  Settings _settings;
  Kinematics _kin;

  int _idx_RECKin;
  int _ix;
  int _iQ2;
  int _iy;
  int _inu;
  int _iW;


  std::string _outfilename="";
  TFile *_tfile;
  TTree *_tree_MC;
  TTree *_tree_Reco;

  double _electron_beam_energy;

  std::map<std::string, float> _map_event;
  std::map<SIDISParticle::PROPERTY,std::vector<float>> _map_particle;

  /* Get true particle info from HIPO Bank */
  int CollectParticlesFromTruth( const std::unique_ptr<clas12::clas12reader>&,
				 type_map_part& );


  /* Get reco particle info from HIPO Bank */
  int CollectParticlesFromReco( const std::unique_ptr<clas12::clas12reader>&,
				type_map_part& );

  /* Connect the pindex of the Reco to the MC map */
  int ConnectTruth2Reco( type_map_part&, type_map_part&);

  /* Write particle information to the tree */
  int WriteParticlesToTree ( type_map_part& );


  /* Add truth event information */
  int AddTruthEventInfo(const std::unique_ptr<clas12::clas12reader>&);

  /* Add reco event information */
  int AddRecoEventInfo(const std::unique_ptr<clas12::clas12reader>&);

  
  /* Reset branch maps for each event */
  void ResetBranchMap();

  /* Clas12root objects */
  clas12root::HipoChain _chain;
  clas12::clas12reader *_config_c12{nullptr};

};
#endif
