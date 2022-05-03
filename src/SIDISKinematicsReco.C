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

using namespace std;

SIDISKinematicsReco::SIDISKinematicsReco(std::string outfilename):
  _ievent(0),
  _tfile(nullptr),
  _tree_MC(nullptr),
  _tree_Reco(nullptr),
  _electron_beam_energy(12)
{
  _outfilename = outfilename;
  cout << "Initialized SIDISKinematicsReco" << endl;
}

int SIDISKinematicsReco::Init()
{

  // Create TFile
  // -------------------------

  _tfile = new TFile(_outfilename.c_str(),"RECREATE");

  // Create event variable map 
  // -------------------------  
  double dummy = 0;
  std::vector<double> vdummy;

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
  _map_particle.insert( make_pair( SIDISParticle::PROPERTY::part_pid , vdummy) );
  _map_particle.insert( make_pair( SIDISParticle::PROPERTY::part_pt , vdummy) );
  _map_particle.insert( make_pair( SIDISParticle::PROPERTY::part_pz , vdummy) );
  _map_particle.insert( make_pair( SIDISParticle::PROPERTY::part_E , vdummy) );
   
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
  for (map< SIDISParticle::PROPERTY, std::vector<double> >::iterator it = _map_particle.begin(); it!=_map_particle.end(); ++it)
    {
      _tree_MC->Branch( SIDISParticle::get_property_info( (it->first) ).first.c_str(), &(it->second));
    }
 
  // Create Reconstructed TTree
  // as Monte Carlo copy
  // -------------------------
  _tree_Reco = _tree_MC->CloneTree();
  _tree_Reco->SetName("event_reco");
  _tree_Reco->SetTitle("A Tree with *reconstructred* event and particle information");

  return 0;
}

int SIDISKinematicsReco::process_event()
{
  // Parse through true Monte Carlo data
  if(_settings.doMC())
    {
      // Reset the branch map
      ResetBranchMap();

      // Using the particle's pindex as a key (should be unique for each reco particle, but could be wrong)
      type_map_part particleMap;
      
      /* Add particle information */
      //      CollectParticlesFromTruth( particleMap );

      /* Write particle information to Tree */
      //WriteParticlesToTree( particleMap );
      
      /* Add event information */
      //      AddTruthEventInfo();
      
      /* Fill MC tree */
      _tree_MC->Fill();
    }
  
  /* Increase event # */
  _ievent++;

  return 0;
}

int SIDISKinematicsReco::WriteParticlesToTree( type_map_part& particleMap )
{
  /* Get number of particles */
  ( _map_event.find("nParticles") )->second = particleMap.size();
  
  /* Loop over all Particles and add them to tree */
  for (type_map_part::iterator it = particleMap.begin(); it!= particleMap.end(); ++it)
    {
      
      for (map< SIDISParticle::PROPERTY , std::vector<double> >::iterator it_prop = _map_particle.begin(); it_prop!= _map_particle.end(); ++it_prop)
	{
	  // Not confident about this section
	  //(it_prop->second).push_back( (it->second) ); 
	}
    }
  
  return 0;
}
void SIDISKinematicsReco::ResetBranchMap()
{
  // Event branch
  for(map<std::string,double>::iterator it = _map_event.begin(); it!=_map_event.end();++it)
    {
      (it->second) = NAN;
    } 
 
  // Particle branch
  for(map<SIDISParticle::PROPERTY , std::vector<double>>::iterator it = _map_particle.begin(); it!=_map_particle.end(); ++it)
    {
      (it->second).clear();
    }
  
  return;
}

int SIDISKinematicsReco::End()
{
  _tfile->cd();
  
  if(_settings.doMC())
    _tree_MC->Write();
  
  if(_settings.doReco())
    _tree_Reco->Write();

  _tfile->Close();
  
  return 0;
}












/*
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
*/

