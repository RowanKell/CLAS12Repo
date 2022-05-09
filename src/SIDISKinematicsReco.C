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
    if(abs(beta) > _settings.getBetamax_fromPID(pid))
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
    
      /* Loop over all MC particles */
      for(type_map_part::iterator it_mc = particleMap.begin(); it_mc!= particleMap.end(); ++it_mc){
	double mc_theta = (it_mc->second)->get_property_float(SIDISParticle::part_theta); 
	double mc_phi = (it_mc->second)->get_property_float(SIDISParticle::part_phi); 

	/* Match the *theta* and *phi* of two particles. For details, see https://www.jlab.org/Hall-B/general/thesis/THayward_thesis.pdf */
	double dth = abs(reco_theta-mc_theta);
	double dphi = abs(reco_phi-mc_phi);
	
	if( (dth < 6*degtorad) && 
	    (dphi < 2*degtorad || abs(dphi - 2*PI) < 2*degtorad)){
	  (it_mc->second)->set_property( SIDISParticle::part_pindex, (it_reco->second)->get_property_int(SIDISParticle::part_pindex));
	  (it_reco->second)->set_property( SIDISParticle::part_parentID , (it_mc->second)->get_property_int(SIDISParticle::part_parentID));
	  (it_reco->second)->set_property( SIDISParticle::part_parentPID , (it_mc->second)->get_property_int(SIDISParticle::part_parentPID));
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
  
  std::vector<float> Mgg;
  _tree_PostProcess->Branch("Mdiphoton",&Mgg);

  // Declaration of leaf types                                                                                                                                                Int_t           event_truth;
  Float_t         Q2;
  Float_t         W;
  Float_t         nParticles;
  Float_t         nPhotons;
  Float_t         nu;
  Float_t         x;
  Float_t         y;
  vector<float>   *pid;
  vector<float>   *px;
  vector<float>   *py;
  vector<float>   *pz;
  vector<float>   *pt;
  vector<float>   *p;
  vector<float>   *E;
  vector<float>   *theta;
  vector<float>   *eta;
  vector<float>   *phi;
  vector<float>   *vz;
  vector<float>   *pindex;
  vector<float>   *beta;
  vector<float>   *chi2;
  vector<float>   *parentID;
  vector<float>   *parentPID; // Declaration of leaf types                                                                                                                    
  // List of branches      
  // TBranch        *b_event;                                                                                                                  
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

  // Set object pointer
  pid = 0;
  px = 0;
  py = 0;
  pz = 0;
  pt = 0;
  p = 0;
  E = 0;
  theta = 0;
  eta = 0;
  phi = 0;
  vz = 0;
  pindex = 0;
  beta = 0;
  chi2 = 0;
  parentID = 0;
  parentPID = 0;

  _tree_Reco->SetBranchAddress("Q2", &Q2, &b_Q2);
  _tree_Reco->SetBranchAddress("W", &W, &b_W);
  _tree_Reco->SetBranchAddress("nParticles", &nParticles, &b_nParticles);
  _tree_Reco->SetBranchAddress("nPhotons", &nPhotons, &b_nPhotons);
  _tree_Reco->SetBranchAddress("nu", &nu, &b_nu);
  _tree_Reco->SetBranchAddress("x", &x, &b_x);
  _tree_Reco->SetBranchAddress("y", &y, &b_y);
  _tree_Reco->SetBranchAddress("pid", &pid, &b_pid);
  _tree_Reco->SetBranchAddress("px", &px, &b_px);
  _tree_Reco->SetBranchAddress("py", &py, &b_py);
  _tree_Reco->SetBranchAddress("pz", &pz, &b_pz);
  _tree_Reco->SetBranchAddress("pt", &pt, &b_pt);
  _tree_Reco->SetBranchAddress("p", &p, &b_p);
  _tree_Reco->SetBranchAddress("E", &E, &b_E);
  _tree_Reco->SetBranchAddress("theta", &theta, &b_theta);
  _tree_Reco->SetBranchAddress("eta", &eta, &b_eta);
  _tree_Reco->SetBranchAddress("phi", &phi, &b_phi);
  _tree_Reco->SetBranchAddress("vz", &vz, &b_vz);
  _tree_Reco->SetBranchAddress("pindex", &pindex, &b_pindex);
  _tree_Reco->SetBranchAddress("beta", &beta, &b_beta);
  _tree_Reco->SetBranchAddress("chi2", &chi2, &b_chi2);
  _tree_Reco->SetBranchAddress("parentID", &parentID, &b_parentID);
  _tree_Reco->SetBranchAddress("parentPID", &parentPID, &b_parentPID);

  Long64_t nentries = _tree_Reco->GetEntriesFast();
  
  TLorentzVector init_electron;
  init_electron.SetPxPyPzE(0,0,sqrt(_electron_beam_energy*_electron_beam_energy - electronMass * electronMass),_electron_beam_energy);

  TLorentzVector electron;
  TLorentzVector q;
  TLorentzVector piplus;
  TLorentzVector pi0;
  TLorentzVector gamma1;
  TLorentzVector gamma2;

  double vz_electron = 0.0;
  double vz_piplus = 0.0;
  double vz_pi0 = 0.0;

  double xF_piplus = 0.0;
  double xF_pi0 = 0.0;
  double zpair = 0.0;

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    //    Long64_t ientry = LoadTree(jentry);
    //    if (ientry < 0) break;
    nb = _tree_Reco->GetEntry(jentry);   nbytes += nb;
    
    Mgg.clear();

    for(unsigned int i = 0 ; i < pid->size() ; i++){
      // Identify the scattered electron
      if(pid->at(i)==11){
	electron.SetPxPyPzE(px->at(i),py->at(i),pz->at(i),E->at(i));
	vz_electron=vz->at(i);
      }
      if(pid->at(i)==211){
	piplus.SetPxPyPzE(px->at(i),py->at(i),pz->at(i),E->at(i));
	vz_piplus=vz->at(i);      
      }
    }
    // Set virtual photon
    q = init_electron - electron;
    
    // Get xF of Pi+
    xF_piplus=2*(piplus*q)/(q.P()*W);

    // Next, identify pairs of photons
    for(unsigned int i = 0 ; i < pid->size() ; i++){
      if(pid->at(i)==22){
	for(unsigned int j = i+1 ; j < pid->size() ; j++){
	  if(pid->at(j)==22){
	    gamma1.SetPxPyPzE(px->at(i),py->at(i),pz->at(i),E->at(i));
	    gamma2.SetPxPyPzE(px->at(j),py->at(j),pz->at(j),E->at(j));
	    pi0=gamma1+gamma2;
	    vz_pi0=(vz->at(i)+vz->at(j))/2.0;
	    xF_pi0=2*(pi0*q)/(q.P()*W);
	    zpair=(piplus.E()+pi0.E())/(nu);
	    if(gamma1.Angle(electron.Vect())>8*PI/180.0 &&
	       gamma2.Angle(electron.Vect())>8*PI/180.0 &&
	       xF_piplus>0 && xF_pi0>0 &&
	       zpair<0.95 &&
	       abs(vz_electron-vz_piplus)<20 && abs(vz_electron-vz_pi0)<20){
	      
	      
	      // All cuts are addressed, now appended interesting quantities
	      Mgg.push_back((gamma1+gamma2).M());
	    }
	  }
	}
      }
    }
    _tree_PostProcess->Fill();
  }
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
