#include "PostProcess.h"
PostProcess::PostProcess(){}
int PostProcess::Init(TTree * tree_Reco, double eBeam)
{
  _electron_beam_energy = eBeam;
  _tree_Reco = tree_Reco;
  _tree_Reco->SetBranchAddress("Q2", &Q2, &b_Q2);
  _tree_Reco->SetBranchAddress("W", &W, &b_W);
  _tree_Reco->SetBranchAddress("nParticles", &nParticles, &b_nParticles);
  _tree_Reco->SetBranchAddress("nPhotons", &nPhotons, &b_nPhotons);
  _tree_Reco->SetBranchAddress("nu", &nu, &b_nu);
  _tree_Reco->SetBranchAddress("x", &x, &b_x);
  _tree_Reco->SetBranchAddress("y", &y, &b_y);
  _tree_Reco->SetBranchAddress("helicity", &helicity, &b_helicity);
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
  _nentries = _tree_Reco->GetEntriesFast();
  return 0;
}

int PostProcess::Process(TTree *_tree_postprocess){
  switch(_process_id){
  case PROCESS_ID::pipluspi0:
    pipi0(_tree_postprocess);
    break;
  case PROCESS_ID::piminuspi0:
    pipi0(_tree_postprocess);
    break;
  case PROCESS_ID::pipluspi0_MC:
    pipi0_MC(_tree_postprocess);
    break;
  case PROCESS_ID::piminuspi0_MC:
    pipi0_MC(_tree_postprocess);
    break;
  }
  return 0;
}

int PostProcess::pipi0(TTree * _tree_postprocess){
  std::vector<float> Mgg;
  std::vector<float> E1;
  std::vector<float> E2;
  std::vector<float> Mdihadron;
  std::vector<float> beta1;
  std::vector<float> beta2;

  _tree_postprocess->Branch("Mdiphoton",&Mgg);
  _tree_postprocess->Branch("E1",&E1);
  _tree_postprocess->Branch("E2",&E2);
  _tree_postprocess->Branch("Mdihadron",&Mdihadron);
  _tree_postprocess->Branch("beta1",&beta1);
  _tree_postprocess->Branch("beta2",&beta2);
  _tree_postprocess->Branch("helicity",&helicity);

  TLorentzVector init_electron;
  init_electron.SetPxPyPzE(0,0,sqrt(_electron_beam_energy*_electron_beam_energy - electronMass * electronMass),_electron_beam_energy);
  TLorentzVector init_target;
  init_target.SetPxPyPzE(0,0,0,protonMass);

  TLorentzVector electron;
  TLorentzVector q;
  TLorentzVector pi;
  TLorentzVector pi0;
  TLorentzVector gamma1;
  TLorentzVector gamma2;
  TLorentzVector dihadron;

  double vz_electron = 0.0;
  double vz_pi = 0.0;
  double vz_pi0 = 0.0;

  double xF_pi = 0.0;
  double xF_pi0 = 0.0;
  double zpair = 0.0;

  double zpi = 0.0;
  double zpi0 = 0.0;

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<_nentries;jentry++) {

    nb = _tree_Reco->GetEntry(jentry);   nbytes += nb;
    
    Mgg.clear();
    E1.clear();
    E2.clear();
    Mdihadron.clear();
    beta1.clear();
    beta2.clear();

    for(unsigned int i = 0 ; i < pid->size() ; i++){
      // Identify the scattered electron
      if(pid->at(i)==11){
	electron.SetPxPyPzE(px->at(i),py->at(i),pz->at(i),E->at(i));
	vz_electron=vz->at(i);
      }
      if(pid->at(i)==_piPID){
	pi.SetPxPyPzE(px->at(i),py->at(i),pz->at(i),E->at(i));
	vz_pi=vz->at(i);      
      }
    }

   
    // Set virtual photon
    q = init_electron - electron;
    
    // Get z of Pi+/-
    zpi = (init_target*pi)/(init_target*q);

    // Get xF of Pi+/-
    //    xF_pi=2*(pi*q)/(abs(q.M())*W);
    xF_pi = zpi/W;

    // Next, identify pairs of photons
    for(unsigned int i = 0 ; i < pid->size() ; i++){
      if(pid->at(i)==22){
	for(unsigned int j = i+1 ; j < pid->size() ; j++){
	  if(pid->at(j)==22){
	    gamma1.SetPxPyPzE(px->at(i),py->at(i),pz->at(i),E->at(i));
	    gamma2.SetPxPyPzE(px->at(j),py->at(j),pz->at(j),E->at(j));
	    pi0=gamma1+gamma2;
	    vz_pi0=(vz->at(i)+vz->at(j))/2.0;
	    zpi0 = (init_target*pi0)/(init_target*q);
//	    xF_pi0=2*(pi0*q)/(abs(q.M())*W);
	    xF_pi0 = zpi0/W;
	    
	    zpair = (init_target*(pi+pi0))/(init_target * q);
	    if(gamma1.Angle(electron.Vect())>8*PI/180.0 &&
	       gamma2.Angle(electron.Vect())>8*PI/180.0 &&
	              xF_pi>0 && xF_pi0>0 &&
	              zpair<0.95 &&
	       abs(vz_electron-vz_pi)<20 && abs(vz_electron-vz_pi0)<20){
	            
	      dihadron = pi0+pi;
	      // All cuts are addressed, now appended interesting quantities
	      Mgg.push_back((gamma1+gamma2).M());
	      E1.push_back(gamma1.E());
	      E2.push_back(gamma2.E());
	      Mdihadron.push_back(dihadron.M());
	      beta1.push_back(beta->at(i));
	      beta2.push_back(beta->at(j));
	    }
	  }
	}
      }
    }
    if(E1.size()>0)
      _tree_postprocess->Fill();
  }
  return 0;
}








int PostProcess::pipi0_MC(TTree * _tree_postprocess){
  std::vector<float> Mgg;
  std::vector<float> E1;
  std::vector<float> E2;
  std::vector<float> Mdihadron;
  std::vector<float> beta1;
  std::vector<float> beta2;
  std::vector<int> flag;

  _tree_postprocess->Branch("Mdiphoton",&Mgg);
  _tree_postprocess->Branch("E1",&E1);
  _tree_postprocess->Branch("E2",&E2);
  _tree_postprocess->Branch("Mdihadron",&Mdihadron);
  _tree_postprocess->Branch("beta1",&beta1);
  _tree_postprocess->Branch("beta2",&beta2);
  _tree_postprocess->Branch("helicity",&helicity);
  _tree_postprocess->Branch("flag",&flag);

  TLorentzVector init_electron;
  init_electron.SetPxPyPzE(0,0,sqrt(_electron_beam_energy*_electron_beam_energy - electronMass * electronMass),_electron_beam_energy);
  TLorentzVector init_target;
  init_target.SetPxPyPzE(0,0,0,protonMass);

  TLorentzVector electron;
  TLorentzVector q;
  TLorentzVector pi;
  TLorentzVector pi0;
  TLorentzVector gamma1;
  TLorentzVector gamma2;
  TLorentzVector dihadron;

  double vz_electron = 0.0;
  double vz_pi = 0.0;
  double vz_pi0 = 0.0;

  double xF_pi = 0.0;
  double xF_pi0 = 0.0;
  double zpair = 0.0;

  double zpi = 0.0;
  double zpi0 = 0.0;

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<_nentries;jentry++) {
    nb = _tree_Reco->GetEntry(jentry);   nbytes += nb;

    Mgg.clear();
    E1.clear();
    E2.clear();
    Mdihadron.clear();
    beta1.clear();
    beta2.clear();
    flag.clear();

    for(unsigned int i = 0 ; i < pid->size() ; i++){
      // Identify the scattered electron
      if(pid->at(i)==11){
	electron.SetPxPyPzE(px->at(i),py->at(i),pz->at(i),E->at(i));
	vz_electron=vz->at(i);
      }
      if(pid->at(i)==_piPID){
	pi.SetPxPyPzE(px->at(i),py->at(i),pz->at(i),E->at(i));
	vz_pi=vz->at(i);      
      }
    }

   
    // Set virtual photon
    q = init_electron - electron;
    
    // Get z of Pi+/-
    zpi = (init_target*pi)/(init_target*q);

    // Get xF of Pi+/-
    //    xF_pi=2*(pi*q)/(abs(q.M())*W);
    xF_pi = zpi/W;

    // Next, identify pairs of photons
    for(unsigned int i = 0 ; i < pid->size() ; i++){
      if(pid->at(i)==22){
	for(unsigned int j = i+1 ; j < pid->size() ; j++){
	  if(pid->at(j)==22){
	    gamma1.SetPxPyPzE(px->at(i),py->at(i),pz->at(i),E->at(i));
	    gamma2.SetPxPyPzE(px->at(j),py->at(j),pz->at(j),E->at(j));
	    pi0=gamma1+gamma2;
	    vz_pi0=(vz->at(i)+vz->at(j))/2.0;
	    zpi0 = (init_target*pi0)/(init_target*q);
//	    xF_pi0=2*(pi0*q)/(abs(q.M())*W);
	    xF_pi0 = zpi0/W;
	    
	    zpair = (init_target*(pi+pi0))/(init_target * q);
	    if(gamma1.Angle(electron.Vect())>8*PI/180.0 &&
	       gamma2.Angle(electron.Vect())>8*PI/180.0 &&
	              xF_pi>0 && xF_pi0>0 &&
	              zpair<0.95 &&
	       abs(vz_electron-vz_pi)<20 && abs(vz_electron-vz_pi0)<20){
	            
	      dihadron = pi0+pi;
	      // All cuts are addressed, now appended interesting quantities
	      Mgg.push_back((gamma1+gamma2).M());
	      E1.push_back(gamma1.E());
	      E2.push_back(gamma2.E());
	      Mdihadron.push_back(dihadron.M());
	      beta1.push_back(beta->at(i));
	      beta2.push_back(beta->at(j));	            
	      if(parentID->at(i)==parentID->at(j) && parentPID->at(i)==111)
		flag.push_back(1);
	      else
		flag.push_back(-1);
	    }
	  }
	}
      }
    }
    if(E1.size()>0)
      _tree_postprocess->Fill();
  }
  return 0;
}


