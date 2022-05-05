// ----------------------------
//
// SIDISKinematicsReco.C()
// 
// Written by Gregory Matousek
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


  // Load in HipoFiles
  // -------------------------
  if(_settings.hipoFileStrings().size()==0){
    std::cout << "ERROR in SIDISKinematicsReco::Init() -- No files added to hipoChain. Double-check that 'addHipoFile' in Settings.h was called by processing script" << std::endl;
    return -1;
  }

  for(unsigned int idx = 0 ; idx < _settings.hipoFileStrings().size() ; idx ++){
    _chain.Add(_settings.hipoFileStrings().at(idx).c_str());
  }

  // Configure CLAS12Reader
  // -------------------------
  _config_c12=_chain.GetC12Reader();
  
  // Initialize Hipo file settings
  // -------------------------------
  InitHipo();
  return 0;
}

int SIDISKinematicsReco::InitHipo()
{
  // -----------------------------------------------------
  // Using the _config_c12 object, cut PIDs for each event
  // -----------------------------------------------------
  std::vector<int> finalStatePIDs = _settings.getFinalStatePIDs();
  for(unsigned int idx = 0 ; idx < finalStatePIDs.size(); idx++){
    int pid = finalStatePIDs.at(idx);
    int npid = _settings.getN_fromPID(pid);
    bool exactpid = _settings.isExact_fromPID(pid);
    if(exactpid==true)
      _config_c12->addExactPid(pid,npid);
    else
      _config_c12->addAtLeastPid(pid,npid);
  }


  return 0;
}
int SIDISKinematicsReco::process_events()
{
  // Establish CLAS12 event parser
  // -------------------------
  auto &_c12= _chain.C12ref();
  
  // Move to the next event in the Hipo chain
  while(_chain.Next()==true){
    if(_c12->getDetParticles().empty())
      continue;
    
    // Only analyze event if we have the desired number of detected particles
    //  if(test->getDetParticles().empty())
    //  return 1;

    // Parse through true Monte Carlo data
    if(_settings.doMC())
      {
	// Reset the branch map
	ResetBranchMap();

	// Using the particle's pindex as a key
	type_map_part particleMap;
      
	/* Add particle information */
	CollectParticlesFromTruth( _c12, particleMap );

	/* Write particle information to Tree */
	//WriteParticlesToTree( _c12, particleMap );
      
	/* Add event information */
	//      AddTruthEventInfo( _c12 );
      
	/* Fill MC tree */
	_tree_MC->Fill();
      }
  
    /* Increase event # */
    _ievent++;
  }
  
  std::cout << " All events completed \n Done." << std::endl;

  return 0;
}

int SIDISKinematicsReco::CollectParticlesFromTruth(const std::unique_ptr<clas12::clas12reader>& _c12,
						   type_map_part& particleMap )
{
  // Get std::vector<> of particles in event
  auto particles=_c12->getDetParticles();

  // Loop over all particles
  for(unsigned int idx = 0 ; idx < particles.size() ; idx++){
    auto particle = particles.at(idx);
    // Skip over particles with PIDs that do not match Settings
    particle->getPid();
    //    int part_pid = parts.at(idx)
  }
  return 0;
}


int SIDISKinematicsReco::WriteParticlesToTree(const std::unique_ptr<clas12::clas12reader>& _c12,
					      type_map_part& particleMap )
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

