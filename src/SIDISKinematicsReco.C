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
  _electron_beam_energy(0)
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
  float dummy = 0;
  std::vector<float> vdummy;

  _map_event.insert( make_pair( "nParticles" , dummy ) );
  _map_event.insert( make_pair( "nPhotons" , dummy ) );
  _map_event.insert( make_pair( "x" , dummy ) );
  _map_event.insert( make_pair( "y" , dummy ) );
  _map_event.insert( make_pair( "Q2" , dummy ) );
  _map_event.insert( make_pair( "W" , dummy ) );
  _map_event.insert( make_pair( "nu" , dummy ) );
  
  // Create particle map 
  // -------------------------  
  _map_particle.insert( make_pair( SIDISParticle::PROPERTY::part_pid , vdummy) );
  _map_particle.insert( make_pair( SIDISParticle::PROPERTY::part_px , vdummy) );
  _map_particle.insert( make_pair( SIDISParticle::PROPERTY::part_py , vdummy) );
  _map_particle.insert( make_pair( SIDISParticle::PROPERTY::part_pz , vdummy) );
  _map_particle.insert( make_pair( SIDISParticle::PROPERTY::part_pt , vdummy) );
  _map_particle.insert( make_pair( SIDISParticle::PROPERTY::part_p , vdummy) );
  _map_particle.insert( make_pair( SIDISParticle::PROPERTY::part_E , vdummy) );
  _map_particle.insert( make_pair( SIDISParticle::PROPERTY::part_evtgen_E , vdummy) );
  _map_particle.insert( make_pair( SIDISParticle::PROPERTY::part_theta , vdummy) );
  _map_particle.insert( make_pair( SIDISParticle::PROPERTY::part_eta , vdummy) );
  _map_particle.insert( make_pair( SIDISParticle::PROPERTY::part_phi , vdummy) );
  _map_particle.insert( make_pair( SIDISParticle::PROPERTY::part_vz , vdummy) );
  _map_particle.insert( make_pair( SIDISParticle::PROPERTY::part_pindex , vdummy) );
  _map_particle.insert( make_pair( SIDISParticle::PROPERTY::part_beta , vdummy) );
  _map_particle.insert( make_pair( SIDISParticle::PROPERTY::part_chi2 , vdummy) );
  _map_particle.insert( make_pair( SIDISParticle::PROPERTY::part_parentID , vdummy) );
  _map_particle.insert( make_pair( SIDISParticle::PROPERTY::part_parentPID , vdummy) );
   
  // Create Monte Carlo TTree
  // -------------------------

  _tree_MC = new TTree("tree_MC","A Tree with *true* event and particle information");
  _tree_MC->Branch("event_truth", &_ievent, "event/I");
  
  // Add Event Branches
  // -------------------------
  for(map< string , float >::iterator it = _map_event.begin(); it!= _map_event.end(); ++it)
    {
      _tree_MC->Branch( (it->first).c_str() , &(it->second) );
    }
  
  // Add Particle Branches
  // -------------------------
  for (map< SIDISParticle::PROPERTY, std::vector<float> >::iterator it = _map_particle.begin(); it!=_map_particle.end(); ++it)
    {
      _tree_MC->Branch( SIDISParticle::get_property_info( (it->first) ).first.c_str(), &(it->second));
    }
 
  // Create Reconstructed TTree
  // as Monte Carlo copy
  // -------------------------
  _tree_Reco = _tree_MC->CloneTree();
  _tree_Reco->SetName("tree_reco");
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

  // Set beam energy
  // -------------------------
  _electron_beam_energy = _settings.electronBeamEnergy();

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

  // -----------------------------------------------------
  // Add banks to the _config_c12 object
  // -----------------------------------------------------
  if(_settings.getEventRecoMethod() == Settings::eventRecoMethod::useRecKinematicsBank){
    _idx_RECKin = _config_c12->addBank("REC::Kinematics");
    _ix = _config_c12->getBankOrder(_idx_RECKin,"x");
    _iQ2 = _config_c12->getBankOrder(_idx_RECKin,"Q2"); 
    _iy = _config_c12->getBankOrder(_idx_RECKin,"y");
    _inu = _config_c12->getBankOrder(_idx_RECKin,"nu");
    _iW = _config_c12->getBankOrder(_idx_RECKin,"W");
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
    if(_verbosity > 0 && (_ievent)%_printEvery==0 && _ievent!=0){
      std::cout << _ievent << " events completed | " << _tree_Reco->GetEntriesFast() << " passed cuts --> " << _tree_Reco->GetEntriesFast()*100.0/_ievent << "%" << std::endl;
    }
    
    /* Increase event # */
    _ievent++;

    if(_c12->getDetParticles().empty())
      continue;
    
    // Create map to store reconstructed SIDIS particles
    // Use the pindex as a unique key
    type_map_part recoparticleMap;

    // Create map to store true SIDIS particles
    type_map_part particleMap;

    
    // Parse through reconstructed particle data
    if(_settings.doReco())
      {
		
	/* Add reco particle information */
	/* Skip event if certain cuts are not satisfied */
	if(CollectParticlesFromReco( _c12, recoparticleMap )!=0)
       	  continue;
      }
    
    // Parse through true Monte Carlo particle data
    if(_settings.doMC())
      {
	
	/* Add particle information */
	CollectParticlesFromTruth( _c12, particleMap );
	
	/* Connect MC and Reco particle info in TTrees */
	if(_settings.connectMC2Reco())
	  {
	    ConnectTruth2Reco( particleMap, recoparticleMap );
	  }
      }

    // ---------------------------
    // Writing to TTrees 
    // ---------------------------

    if(_settings.doReco())
      {
	/* Reset branch map */
	ResetBranchMap();

	/* Add event information */
	/* Skip event if it fails cuts */
	if(AddRecoEventInfo( _c12 )!=0)
	  continue;
	/* Write particle information to Tree */
       	WriteParticlesToTree( recoparticleMap );

	/* fill reco tree */
       	_tree_Reco->Fill();
      }
    if(_settings.doMC())
      {
	/* Reset branch map */
	ResetBranchMap();

	/* Add event information */
	AddTruthEventInfo( _c12 );
	/* Write particle information to Tree */	
	WriteParticlesToTree( particleMap );
	
	/* Fill MC tree */
	_tree_MC->Fill();
      }
    
  }
  
  std::cout << " All events completed \n Done." << std::endl;

  return 0;
}

int SIDISKinematicsReco::CollectParticlesFromTruth(const std::unique_ptr<clas12::clas12reader>& _c12,
						   type_map_part& particleMap )
{
  // Get pointer to Monte Carlo particles in event
  auto mcparticles=_c12->mcparts();

  // Loop over all particles
  for(int idx = 0 ; idx < mcparticles->getRows() ; idx++){
    
    // Get particle in MC::Lund
    mcparticles->setEntry(idx);
    if(mcparticles->getType()!=1) // Reject non-final state
      continue; 
    
    // Create new SIDISParticle
    SIDISParticlev1 *sp = new SIDISParticlev1();
    sp->set_candidate_id( particleMap.size() );

    int pid = mcparticles->getPid();
    float px = mcparticles->getPx();
    float py = mcparticles->getPy();
    float pz = mcparticles->getPz();
    float m = mcparticles->getMass();
    
    float pt = _kin.Pt(px,py);
    float p  = _kin.P(px,py,pz);
    float E  = _kin.E(m,p);

    float theta = _kin.th(pt,pz);
    float eta = _kin.eta(theta);
    float phi   = _kin.phi(px,py);

    float vz = mcparticles->getVz();

    sp->set_property( SIDISParticle::part_pid, pid);
    sp->set_property( SIDISParticle::part_px,  px);
    sp->set_property( SIDISParticle::part_py,  py);
    sp->set_property( SIDISParticle::part_pz,  pz);
    sp->set_property( SIDISParticle::part_pt,  pt);
    sp->set_property( SIDISParticle::part_p,  p);
    sp->set_property( SIDISParticle::part_E,   E);
    sp->set_property( SIDISParticle::part_evtgen_E,   E);
    sp->set_property( SIDISParticle::part_theta,   theta);
    sp->set_property( SIDISParticle::part_eta,   eta);
    sp->set_property( SIDISParticle::part_phi,   phi);

    sp->set_property( SIDISParticle::part_vz,   vz);

    sp->set_property( SIDISParticle::part_pindex,   -999);
    sp->set_property( SIDISParticle::part_beta,   -999);
    sp->set_property( SIDISParticle::part_chi2,   -999);
    sp->set_property( SIDISParticle::part_ID,   mcparticles->getIndex());
    sp->set_property( SIDISParticle::part_parentID,  mcparticles->getParent());
    sp->set_property( SIDISParticle::part_parentPID, mcparticles->getPid(mcparticles->getParent()-1));
    
    // Add SIDISParticle to the collection
    particleMap.insert ( make_pair( sp->get_candidate_id() , sp) );
    
  }
  return 0;
}


int SIDISKinematicsReco::CollectParticlesFromReco(const std::unique_ptr<clas12::clas12reader>& _c12,
						   type_map_part& recoparticleMap )
{
  // Get std::vector<> of particles in event
  auto particles=_c12->getDetParticles();

  // Loop over all particles
  for(unsigned int idx = 0 ; idx < particles.size() ; idx++){

    // Extract each particle from event one-at-a-time
    auto particle = particles.at(idx);
    int pid = particle->getPid();

    // CUT pid -------------------------------------------------------------
    // Skip over particles that are not interesting in the final state
    if(_settings.ignoreOtherRecoParticles() && _settings.getN_fromPID(pid)==0)
      continue;

    float chi2 = particle->getChi2Pid();

    // CUT chi2 -------------------------------------------------------------
    // Skip over particles that both need a chi2pid cut, and do not satisfy it
    if(abs(chi2) > _settings.getChi2max_fromPID(pid))
      {
	// Skip event if this lone particle NEEDED to pass the cut
	if(_settings.getN_fromPID(pid)==1 && _settings.isExact_fromPID(pid))
	  return -1;
	else
	  continue;
      }

    float beta = particle->getBeta();

    // CUT beta -------------------------------------------------------------
    // Skip over particles that both need a beta cut, and do not satisfy it
    if(abs(beta) > _settings.getBetamax_fromPID(pid) || abs(beta) < _settings.getBetamin_fromPID(pid))
      {
	// Skip event if this lone particle NEEDED to pass the cut
	if(_settings.getN_fromPID(pid)==1 && _settings.isExact_fromPID(pid))
	  return -1;
	else
	  continue;
      }

    float p = particle->getP();
    // CUT p -------------------------------------------------------------
    // Skip over particles that do not satisfy minimum momentum cut

    if(_settings.getPmin_fromPID(pid) > p)
      {
	// Skip event if this lone particle NEEDED to pass the cut
	if(_settings.getN_fromPID(pid)==1 && _settings.isExact_fromPID(pid))
	  return -1;
	else
	  continue;

      }
    float theta = particle->getTheta();
    float eta = _kin.eta(theta);
    float phi = particle->getPhi();
    float px = _kin.Px(p,theta,phi);
    float py = _kin.Py(p,theta,phi);
    float pz = _kin.Pz(p,theta,phi);
    float pt = _kin.Pt(px,py);
    float m = 0.0;
    if(pid!=22)
      m = particle->getPdgMass();

    float E  = _kin.E(m,p);
    // CUT E -------------------------------------------------------------
    // Skip over particles that do not satisfy minimum energy cut
    if(_settings.getEmin_fromPID(pid) > E)
      continue;

    int pindex = particle->getIndex();
    float vz = particle->par()->getVz();

    // CUT vz -------------------------------------------------------------
    // Skip over particles that do not satisfy Vz range
    if(vz < _settings.getVzmin_fromPID(pid) || vz > _settings.getVzmax_fromPID(pid))
      {
	// Skip event if this lone particle NEEDED to pass the cut
	if(_settings.getN_fromPID(pid)==1 && _settings.isExact_fromPID(pid))
	  return -1;
	else
	  continue;
      }

    SIDISParticlev1 *sp = new SIDISParticlev1();
    sp->set_candidate_id( pindex );
    
    sp->set_property( SIDISParticle::part_pid, pid);
    sp->set_property( SIDISParticle::part_px,  px);
    sp->set_property( SIDISParticle::part_py,  py);
    sp->set_property( SIDISParticle::part_pz,  pz);
    sp->set_property( SIDISParticle::part_pt,  pt);
    sp->set_property( SIDISParticle::part_p,  p);
    sp->set_property( SIDISParticle::part_E,   E);
    sp->set_property( SIDISParticle::part_evtgen_E,  -999);
    sp->set_property( SIDISParticle::part_theta,   theta);
    sp->set_property( SIDISParticle::part_eta,   eta);
    sp->set_property( SIDISParticle::part_phi,   phi);

    sp->set_property( SIDISParticle::part_vz,   vz);

    sp->set_property( SIDISParticle::part_pindex,   pindex);
    sp->set_property( SIDISParticle::part_beta,   beta);
    sp->set_property( SIDISParticle::part_chi2,   chi2);
    sp->set_property( SIDISParticle::part_ID, pindex);
    sp->set_property( SIDISParticle::part_parentID, -999);
    sp->set_property( SIDISParticle::part_parentPID, -999);
    // Add SIDISParticle to the collection
    recoparticleMap.insert ( make_pair( sp->get_candidate_id() , sp) );

  }
  return 0;
}

int SIDISKinematicsReco::ConnectTruth2Reco( type_map_part& particleMap,
					    type_map_part& recoparticleMap )
{
  /* Loop over all reco particles */
  for (type_map_part::iterator it_reco = recoparticleMap.begin(); it_reco!= recoparticleMap.end() ; ++it_reco){

    double reco_theta = (it_reco->second)->get_property_float(SIDISParticle::part_theta); 
    double reco_phi = (it_reco->second)->get_property_float(SIDISParticle::part_phi); 
    //    double reco_E = (it_reco->second)->get_property_float(SIDISParticle::part_E);
      /* Loop over all MC particles */
      for(type_map_part::iterator it_mc = particleMap.begin(); it_mc!= particleMap.end(); ++it_mc){
	double mc_theta = (it_mc->second)->get_property_float(SIDISParticle::part_theta); 
	double mc_phi = (it_mc->second)->get_property_float(SIDISParticle::part_phi); 
	//	double mc_E = (it_reco->second)->get_property_float(SIDISParticle::part_E);
	/* Match the *theta* and *phi* of two particles. For details, see https://www.jlab.org/Hall-B/general/thesis/THayward_thesis.pdf */
	double dth = abs(reco_theta-mc_theta);
	double dphi = abs(reco_phi-mc_phi);
	//	double dE = abs(reco_E - mc_E);

	if( (dth < 6*degtorad) && 
	    (dphi < 2*degtorad || abs(dphi - 2*PI) < 2*degtorad)){
	  (it_mc->second)->set_property( SIDISParticle::part_pindex, (it_reco->second)->get_property_int(SIDISParticle::part_pindex));
	  (it_reco->second)->set_property( SIDISParticle::part_parentID , (it_mc->second)->get_property_int(SIDISParticle::part_parentID));
	  (it_reco->second)->set_property( SIDISParticle::part_parentPID , (it_mc->second)->get_property_int(SIDISParticle::part_parentPID));
	  (it_reco->second)->set_property( SIDISParticle::part_evtgen_E , (it_mc->second)->get_property_float(SIDISParticle::part_evtgen_E));
	}
      }
  }
  
  return 0;
  
}

int SIDISKinematicsReco::AddTruthEventInfo(const std::unique_ptr<clas12::clas12reader>& _c12)
{
  // Get pointer to Monte Carlo particles in event
  auto mcparticles=_c12->mcparts();

  // Loop over all particles
  for(int idx = 0 ; idx < mcparticles->getRows() ; idx++){
   
    // Get particle in MC::Lund
    mcparticles->setEntry(idx);
    
    // Test if scattered lepton
    if(mcparticles->getParent()==1 && mcparticles->getPid()==11){
      float px = mcparticles->getPx();
      float py = mcparticles->getPy();
      float pz = mcparticles->getPz();
      float m = mcparticles->getMass();
      float p = _kin.P(px,py,pz);
      float E = _kin.E(m,p);
      float cth = _kin.cth(px,py,pz);
      float Q2 = _kin.Q2(_electron_beam_energy,E,cth);
      float y  = _kin.y(_electron_beam_energy, E);
      float nu  = _kin.nu(_electron_beam_energy, E);
      float W  = _kin.W(Q2,protonMass,nu);
      float s  = protonMass*protonMass+electronMass*electronMass+2*protonMass*_electron_beam_energy;
      float x = _kin.x(Q2,s,y);

      (_map_event.find("x"))->second = x;
      (_map_event.find("Q2"))->second = Q2;
      (_map_event.find("y"))->second = y;
      (_map_event.find("nu"))->second = nu;
      (_map_event.find("W"))->second = W;

      break;
    }
  }
  
  return 0;
}


int SIDISKinematicsReco::AddRecoEventInfo(const std::unique_ptr<clas12::clas12reader>& _c12)
{

  double reco_x = -999;
  double reco_Q2 = -999;
  double reco_y = -999;
  double reco_nu = -999;
  double reco_W = -999;

  /* METHOD 1: If available, use REC::Kinematics bank for event reco */
  if(_settings.getEventRecoMethod() == Settings::eventRecoMethod::useRecKinematicsBank){
    reco_x=_c12->getBank(_idx_RECKin)->getDouble(_ix,0);
    reco_Q2=_c12->getBank(_idx_RECKin)->getDouble(_iQ2,0);
    reco_y=_c12->getBank(_idx_RECKin)->getDouble(_iy,0);
    reco_nu=_c12->getBank(_idx_RECKin)->getDouble(_inu,0);
    reco_W=_c12->getBank(_idx_RECKin)->getDouble(_iW,0);
  }
  /* METHOD 2: Use the scattered electron to reconstruct event variables */
  else if(_settings.getEventRecoMethod() == Settings::eventRecoMethod::useLargestPinFD){
    // Get std::vector<> of particles in event
    auto particles=_c12->getDetParticles();
    
    // Loop over all particles
    for(unsigned int idx = 0 ; idx < particles.size() ; idx++){
      // Extract each particle from event one-at-a-time
      auto particle = particles.at(idx);
      if(particle->getStatus()<0){
	float p = particle->getP();
	float theta = particle->getTheta();
	float phi = particle->getPhi();
	float px = _kin.Px(p,theta,phi);
	float py = _kin.Py(p,theta,phi);
	float pz = _kin.Pz(p,theta,phi);
    	float m = particle->getCalcMass();
	float E = _kin.E(m,p);
	float cth = _kin.cth(px,py,pz);
	reco_Q2 = _kin.Q2(_electron_beam_energy,E,cth);
	reco_y  = _kin.y(_electron_beam_energy, E);
	reco_nu  = _kin.nu(_electron_beam_energy, E);
	reco_W  = _kin.W(reco_Q2,protonMass,reco_nu);
	float s  = protonMass*protonMass+electronMass*electronMass+2*protonMass*_electron_beam_energy;
	reco_x = _kin.x(reco_Q2,s,reco_y);
	
	// Found scattered electron, just break out of here
	break;
      }
    }
  }
  
  if(reco_W < _settings.Wmin() || reco_W > _settings.Wmax())
    return -1;
  else if(reco_y < _settings.ymin() || reco_y > _settings.ymax())
    return -1;
  else if(reco_Q2 < _settings.Q2min() || reco_Q2 > _settings.Q2max())
    return -1;
  else{
    (_map_event.find("x"))->second = reco_x; 
    (_map_event.find("Q2"))->second = reco_Q2;
    (_map_event.find("y"))->second = reco_y;
    (_map_event.find("nu"))->second = reco_nu;
    (_map_event.find("W"))->second = reco_W;
  }
  
  return 0;
}


int SIDISKinematicsReco::WriteParticlesToTree(type_map_part& particleMap )
{
  /* Get number of particles */
  ( _map_event.find("nParticles") )->second = particleMap.size();
  
  /* Loop over all Particles and add them to tree */
  for (type_map_part::iterator it = particleMap.begin(); it!= particleMap.end(); ++it)
    {
      
      for (map< SIDISParticle::PROPERTY , std::vector<float> >::iterator it_prop = _map_particle.begin(); it_prop!= _map_particle.end(); ++it_prop)
	{
	  switch ( SIDISParticle::get_property_info( (it_prop->first) ).second ) {
	    
          case SIDISParticle::type_float :
            (it_prop->second).push_back( (it->second)->get_property_float( (it_prop->first) ) );
            break;

          case SIDISParticle::type_int :
            (it_prop->second).push_back( (it->second)->get_property_int( (it_prop->first) ) );
            break;
	  }
	}
    }
  
  return 0;
}
void SIDISKinematicsReco::ResetBranchMap()
{
  // Event branch
  for(map<std::string,float>::iterator it = _map_event.begin(); it!=_map_event.end();++it)
    {
      (it->second) = NAN;
    }   

  // Particle branch
  for(map<SIDISParticle::PROPERTY , std::vector<float>>::iterator it = _map_particle.begin(); it!=_map_particle.end(); ++it)
    {
      (it->second).clear();
    }
  return;
}


int SIDISKinematicsReco::PostProcessReco()
{

  std::cout << "----------------------------------------" << std::endl;
  std::cout << "Beginning of PostProcessing. Total events = " << _tree_Reco->GetEntriesFast() << std::endl;
  std::cout << "----------------------------------------" << std::endl;
 
  _tree_PostProcess = new TTree("tree_postprocess","Tree created after all hipo events are filtered");

  _postprocess.Init(_tree_Reco,_electron_beam_energy);
  
  if(_postprocess.setPostProcessMethod(_settings.postProcessMethod())!=0)
    return -1;
  
  // This line is quite the mouthful...
  _postprocess.Process(_tree_PostProcess);

  
  return 0;
}




int SIDISKinematicsReco::End()
{

  std::cout << "----------------------------------------" << std::endl;
  std::cout << "Number of events saved to TTRee = " << _tree_Reco->GetEntriesFast() << std::endl;
  std::cout << "----------------------------------------" << std::endl;
  _tfile->cd();
  
  if(_settings.doMC())
    _tree_MC->Write();
  
  if(_settings.doReco())
    _tree_Reco->Write();

  if(_settings.doPostProcess())
    _tree_PostProcess->Write();

  _tfile->Close();
  
  return 0;
}
