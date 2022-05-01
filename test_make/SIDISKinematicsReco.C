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

#include "SIDISKinematicsReco.h"

SIDISKinematicsReco::SIDISKinematicsReco(std::string cutcardname):
  _ievent(0),
  _cutcardname(cutcardname),
  _tfile(nullptr),
  _tree_MC(nullptr),
  _tree_Reco(nullptr),
  _electron_beam_energy(12),
  _do_MC(false),
  _do_Reco(true)
{

}

int SIDISKinematicsReco::Init()
{

  // Open settings from CutCard
  // -------------------------
  
  if(!_settings.readCard(_cutcardname.c_str())){
    return -1;
  }
  else{
    _settings.loadSettings();
  }

  // Create TFile
  // -------------------------

  _tfile = new TFile(_settings.Filename().c_str(),"RECREATE");

  // Create event variable map 
  // -------------------------  
  double dummy = 0;
  std::vector<float> vdummy;

  _map_event.insert( make_pair( "nParticles" , dummy ) );
  _map_event.insert( make_pair( "nPhotons" , dummy ) );
  _map_event.insert( make_pair( "s" , dummy ) );
  _map_event.insert( make_pair( "x" , dummy ) );
  _map_event.insert( make_pair( "y" , dummy ) );
  _map_event.insert( make_pair( "Q2" , dummy ) );
  _map_event.insert( make_pair( "W" , dummy ) );
  _map_event.insert( make_pair( "nu" , dummy ) );
  
  // Create particle map 
  // -------------------------  
  /* _map_particle.insert( make_pair( SIDISParticle::id , vdummy) );
  _map_particle.insert( make_pair( SIDISParticle::Pt , vdummy) );
  _map_particle.insert( make_pair( SIDISParticle::Pz , vdummy) );
  _map_particle.insert( make_pair( SIDISParticle::E , vdummy) );
  */
  // Create Monte Carlo TTree
  // -------------------------

  _tree_MC = new TTree("tree_MC","A Tree with *true* event and particle information");
  _tree_MC->Branch("event_truth", &_ievent, "event/I");
  
  // Add Event Branches
  // -------------------------
  for(map< string , double >::iterator it = _map_event.begin(); it!= _map_event.end(); ++it)
    {
      _tree_MC->Branch( (it->first).c_str() , &(it->second) );
    }
  
  // Add Particle Branches
  // -------------------------
  /*for (map< SIDISParticle::PROPERTY, std::vector<double> >::iterator it = _map_particle.begin(); it!=_map_particle.end(); ++it)
    {
      _tree_MC->Branch( SIDISParticle::get_property_info( (it->first) ).first.c_str(), &(it->second));
    }
  */
  // Create Reconstructed TTree
  // as Monte Carlo copy
  // -------------------------
  _tree_Reco = _tree_MC->CloneTree();
  _tree_Reco->SetName("event_reco");
  _tree_Reco->SetTitle("A Tree with *reconstructred* event and particle information");

  return 0;
}


