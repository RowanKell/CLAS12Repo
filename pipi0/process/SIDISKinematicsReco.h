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

typedef std::map<double, SIDISParticle*> type_map_part;

class SIDISKinematicsReco{
 
 public:
  
  SIDISKinematicsReco(std::string cutcardname);

  int Init();

  int process_event();

  int End();
  
  void set_beam(double E){
    _electron_beam_energy = E;
  }
  
  void set_do_MC(bool select){
    _do_MC = select;
  }

  void set_do_Reco(bool select){
    _do_Reco = select;
  }

 private:
 
  int _ievent;
  std::string _cutcardname;
  Settings _settings;
  Kinematics _kin;
  
  TFile *_tfile;
  TTree *_tree_MC;
  TTree *_tree_Reco;

  double _electron_beam_energy;

  bool _do_MC;
  bool _do_Reco;

  std::map<std::string, double> _map_event;
  std::map< SIDISParticle::PROPERTY , std::vector<double>> _map_particle;

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
