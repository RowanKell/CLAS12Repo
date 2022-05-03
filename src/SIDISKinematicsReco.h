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

#include <TFile.h>
#include <TTree.h>

//class TFile;
//class TTree;
class SIDISParticle;

typedef std::map<int, SIDISParticle*> type_map_part;

class SIDISKinematicsReco{
 
 public:
  
  SIDISKinematicsReco(std::string);

  int Init();

  void ImportSettings(Settings theSettings) {
    _settings = theSettings;
  }

  int process_event();

  int End();
  
  void set_beam(double E){
    _electron_beam_energy = E;
  }
 
 private:
 
  int _ievent;
  Settings _settings;
  Kinematics _kin;

  std::string _outfilename="";
  TFile *_tfile;
  TTree *_tree_MC;
  TTree *_tree_Reco;

  double _electron_beam_energy;

  std::map<std::string, double> _map_event;
  std::map<SIDISParticle::PROPERTY,std::vector<double>> _map_particle;

  /* Get true particle info from HIPO Bank */
  int CollectParticlesFromTruth( type_map_part& );


  /* Get reco particle info from HIPO Bank */
  int CollectParticlesFromReco( type_map_part& );


  /* Write particle information to the tree */
  int WriteParticlesToTree ( type_map_part& );


  /* Add truth event information */
  int AddTruthEventInfo();

  
  /* Reset branch maps for each event */
  void ResetBranchMap();

};
#endif
