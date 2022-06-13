#include "PostProcess.h"
PostProcess::PostProcess(){}
int PostProcess::Init(TTree * tree_Reco, double eBeam)
{
  _electron_beam_energy = eBeam;
  _tree_Reco = tree_Reco;
  _tree_Reco->SetBranchAddress("Q2", &Q2, &b_Q2);
  _tree_Reco->SetBranchAddress("W", &W, &b_W);
  _tree_Reco->SetBranchAddress("nParticles", &nParticles, &b_nParticles);
  _tree_Reco->SetBranchAddress("nu", &nu, &b_nu);
  _tree_Reco->SetBranchAddress("x", &x, &b_x);
  _tree_Reco->SetBranchAddress("y", &y, &b_y);
  _tree_Reco->SetBranchAddress("helicity", &helicity, &b_helicity);
  //  _tree_Reco->SetBranchAddress("charge", &charge, &b_charge);
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
  case PROCESS_ID::pipluspiminus:
    pipiminus(_tree_postprocess);
    break;
  }
  return 0;
}

int PostProcess::pipi0(TTree * _tree_postprocess){
  double fid = 0.0;
  double Mgg;
  double E1;
  double E2;
  double Mdihadron;
  double beta1;
  double beta2;
  double phi_R;
  double phi_h;
  double th;

  _tree_postprocess->Branch("fid",&fid);
  _tree_postprocess->Branch("Mdiphoton",&Mgg);
  _tree_postprocess->Branch("E1",&E1);
  _tree_postprocess->Branch("E2",&E2);
  _tree_postprocess->Branch("Mdihadron",&Mdihadron);
  _tree_postprocess->Branch("beta1",&beta1);
  _tree_postprocess->Branch("beta2",&beta2);
  _tree_postprocess->Branch("helicity",&helicity);
  _tree_postprocess->Branch("phi_R",&phi_R);
  _tree_postprocess->Branch("phi_h",&phi_h);
  _tree_postprocess->Branch("th",&th);
  _tree_postprocess->Branch("Q2",&Q2);
  _tree_postprocess->Branch("W",&W);
  _tree_postprocess->Branch("nu",&nu);
  _tree_postprocess->Branch("x",&x);
  _tree_postprocess->Branch("y",&y);
  //  _tree_postprocess->Branch("charge",&charge);
  
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
    
    Mgg=0.0;
    E1=0.0;
    E2=0.0;
    Mdihadron=0.0;
    beta1=0.0;
    beta2=0.0;
    phi_h=0.0;
    phi_R=0.0;
    th=0.0;

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
	      Mgg=((gamma1+gamma2).M());
	      E1=(gamma1.E());
	      E2=(gamma2.E());
	      Mdihadron=(dihadron.M());
	      beta1=(beta->at(i));
	      beta2=(beta->at(j));
	      phi_R=(_kin.phi_R(q,init_electron,pi,pi0));
	      phi_h=(_kin.phi_h(q,init_electron,pi,pi0));
     	      th = (_kin.com_th(pi,pi0));
	      _tree_postprocess->Fill();
	      fid++;
	    }
	  }
	}
      }
    }
  }
  return 0;
}








int PostProcess::pipi0_MC(TTree * _tree_postprocess){
  double fid = 0.0;
  double Mgg;
  double E1;
  double E2;
  double Mdihadron;
  double beta1;
  double beta2;
  int flag;
  double phi_h;
  double phi_R;
  double th; 

  _tree_postprocess->Branch("fid",&fid);
  _tree_postprocess->Branch("Mdiphoton",&Mgg);
  _tree_postprocess->Branch("E1",&E1);
  _tree_postprocess->Branch("E2",&E2);
  _tree_postprocess->Branch("Mdihadron",&Mdihadron);
  _tree_postprocess->Branch("beta1",&beta1);
  _tree_postprocess->Branch("beta2",&beta2);
  _tree_postprocess->Branch("helicity",&helicity);
  _tree_postprocess->Branch("flag",&flag);
  _tree_postprocess->Branch("phi_R",&phi_R);
  _tree_postprocess->Branch("phi_h",&phi_h);
  _tree_postprocess->Branch("Q2",&Q2);
  _tree_postprocess->Branch("W",&W);
  _tree_postprocess->Branch("nu",&nu);
  _tree_postprocess->Branch("x",&x);
  _tree_postprocess->Branch("y",&y);
  _tree_postprocess->Branch("th",&th);

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

    Mgg=0.0;
    E1=0.0;
    E2=0.0;
    Mdihadron=0.0;
    beta1=0.0;
    beta2=0.0;
    flag=0;
    phi_h=0.0;
    phi_R=0.0;
    th = 0.0;
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
	      Mgg=((gamma1+gamma2).M());
	      E1=(gamma1.E());
	      E2=(gamma2.E());
	      Mdihadron=(dihadron.M());
	      beta1=(beta->at(i));
	      beta2=(beta->at(j));	            
	      if(parentID->at(i)==parentID->at(j) && parentPID->at(i)==111)
		flag=(1);
	      else
		flag=(-1);
	      phi_R=(_kin.phi_R(q,init_electron,pi,pi0));
	      phi_h=(_kin.phi_h(q,init_electron,pi,pi0));	      
     	      th = (_kin.com_th(pi,pi0));
	      _tree_postprocess->Fill();
	      fid++;
	    }
	  }
	}
      }
    }
  }
  return 0;
}


int PostProcess::pipiminus(TTree * _tree_postprocess){
  double fid = 0.0;
//  double Mgg; 
//  double E1;
//  double E2;
  double Mdihadron;
//  double beta1;
//  double beta2;
  double phi_R;
  double phi_h;
  double th;
  double x_F;
  double z_h;
  double P_h;
  double M_x;
  double xF_piplus;
  double xF_piminus;

  _tree_postprocess->Branch("fid",&fid);
//  _tree_postprocess->Branch("Mdiphoton",&Mgg);
//  _tree_postprocess->Branch("E1",&E1);
//  _tree_postprocess->Branch("E2",&E2);
  _tree_postprocess->Branch("Mdihadron",&Mdihadron);
//  _tree_postprocess->Branch("beta1",&beta1);
//  _tree_postprocess->Branch("beta2",&beta2);
  _tree_postprocess->Branch("helicity",&helicity);
  _tree_postprocess->Branch("phi_R",&phi_R);
  _tree_postprocess->Branch("phi_h",&phi_h);
  _tree_postprocess->Branch("th",&th);
  _tree_postprocess->Branch("Q2",&Q2);
  _tree_postprocess->Branch("W",&W);
  _tree_postprocess->Branch("nu",&nu);
  _tree_postprocess->Branch("x",&x);
  _tree_postprocess->Branch("y",&y);
  //  _tree_postprocess->Branch("charge",&charge);
  _tree_postprocess->Branch("x_F",&x_F);
  _tree_postprocess->Branch("z_h",&z_h);
  _tree_postprocess->Branch("P_h",&P_h);
  _tree_postprocess->Branch("M_x",&M_x);
  _tree_postprocess->Branch("xF_piplus",&xF_piplus);
  _tree_postprocess->Branch("xF_piminus",&xF_piminus);
    
  TLorentzVector init_electron;
  init_electron.SetPxPyPzE(0,0,sqrt(_electron_beam_energy*_electron_beam_energy - electronMass * electronMass),_electron_beam_energy);
  TLorentzVector init_target;
  init_target.SetPxPyPzE(0,0,0,protonMass);

  TLorentzVector electron;
  TLorentzVector q;
  TLorentzVector piplus;
  TLorentzVector piminus; // introduced instead of pi0 and gammas
//  TLorentzVector gamma1; Replaced by piminus
//  TLorentzVector gamma2; ""
  TLorentzVector dihadron; //same as lv_p if looking at THaywards github
  TLorentzVector lv_Mx; // vector for missing mass
  TLorentzVector lv_p_gN; // vector for dihadron boosted to gN frame
  TLorentzVector lv_p1_gN; // vector for piplus boosted to gN frame
  TLorentzVector lv_p2_gN;
  TLorentzVector lv_q_gN;
  TLorentzVector lv_e_gN; // vector for electron boosted to gN frame
  TLorentzVector gN; // gamma nucleon vector

  TVector3 qv;
  TVector3 gNBoost;
  TVector3 gNBoostNeg;
  TVector3 lv_e_gN_unit; //Unit vector of electron boosted LV
  TVector3 lv_q_gN_unit;
  TVector3 vectPh;
  TVector3 vectPhT;
  TVector3 Pt_Q;
  TVector3 vecR;
  TVector3 vecRpre;
  TVector3 vecH;
  TVector3 VTH1;
  TVector3 vectPhT1;
  TVector3 vectPh1;
  TVector3 vT;
  TVector3 vTH;
    
  double vz_electron = 0.0;
  double vz_piplus = 0.0; // piplus instead of pi
  double vz_piminus = 0.0; // piminus instead of pi0

//  double xF_piplus = 0.0;
//  double xF_piminus = 0.0;

//  double zpair = 0.0;

  double zpiplus = 0.0;
  double zpiminus = 0.0;

  double sinPhiH = 0.0;
  double hScale = 0.0;
  double cosphih = 0.0;
    
//  double phi_R_pre = 0.0;
//  double phi_h_pre = 0.0;
    
//  double theta_gh = 0.0;//angle between dihadron momentum and virtual photon momentum (NOW UNUSED BC PERP COMMAND IS USED)
    
//  double P_piplus = 0.0;
//  double P_piminus = 0.0;

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<_nentries;jentry++) {

    nb = _tree_Reco->GetEntry(jentry);   nbytes += nb;
    
//    Mgg=0.0;
//    E1=0.0;
//    E2=0.0;
    Mdihadron=0.0;
//    beta1=0.0;
//    beta2=0.0;
    phi_h=0.0;
    phi_R=0.0;
    th=0.0;
    x_F=0.0;
    z_h=0.0;
    P_h=0.0;
    M_x=0.0;

    for(unsigned int i = 0 ; i < pid->size() ; i++){
      // Identify the scattered electron
      if(pid->at(i)==11){
	electron.SetPxPyPzE(px->at(i),py->at(i),pz->at(i),E->at(i));
	vz_electron=vz->at(i);
      }
     // Identify pi+
      if(pid->at(i)==211){ // _piplusPID is defined in PostProcess.h
	piplus.SetPxPyPzE(px->at(i),py->at(i),pz->at(i),E->at(i));
	vz_piplus=vz->at(i);      
      }
     // Identify pi-
      if(pid->at(i)==-211){ // pid -211 is pi-
	piminus.SetPxPyPzE(px->at(i),py->at(i),pz->at(i),E->at(i));
	vz_piminus=vz->at(i);      
      }
    }
    
      
    // Set virtual photon
    q = init_electron - electron;
      
    
      
    // Create Boost vector for xF calculation
    gN = q;
    gN += init_target;
    gNBoost = gN.BoostVector();
    gNBoostNeg = -gNBoost;
    
    // Set boost vectors for pi+pi-
    lv_p1_gN = piplus;
    lv_p2_gN = piminus;

    //Boosting vectors for pi+pi-
    lv_p1_gN.Boost(gNBoostNeg);
    lv_p2_gN.Boost(gNBoostNeg);


      
    //Boosting electron to gN frame
    lv_e_gN = electron;
    lv_e_gN.Boost(gNBoostNeg);
    lv_e_gN_unit.SetMagThetaPhi(1,lv_e_gN.Theta(),lv_e_gN.Phi());
      

      


    //
    //ORIGINAL CODE:
    //
    // Get z of Pi+/-
//  zpi = (init_target*pi)/(init_target*q);
    //
    // Get z of Pi+
    zpiplus = (init_target*piplus)/(init_target*q); //This seems right according to formulas in Affinity paper

    // Get z of Pi-
    zpiminus = (init_target*piminus)/(init_target*q);
      
    // Get the zPair
    z_h = zpiplus + zpiminus;

    //Setting vecH
    vecH.SetMagThetaPhi(lv_p2_gN.Vect().Mag()/zpiminus,lv_p2_gN.Vect().Theta(),lv_p2_gN.Vect().Phi());
      
    //Setting vecR
    vecRpre = lv_p1_gN.Vect();
    vecRpre.SetMagThetaPhi(lv_p1_gN.Vect().Mag()/zpiplus,lv_p1_gN.Vect().Theta(),lv_p1_gN.Vect().Phi());
    vecR = vecRpre - vecH;
      

      
      

    //
    //ORIGINAL CODE:
    //
    // Get xF of Pi+/-
    //    xF_pi=2*(pi*q)/(abs(q.M())*W);
//  xF_pi = zpi/W;
    //
      

    // Get xF pair
//    x_F = z_h / W; //This follows the original xF = z/W calculation but does not produce the desired histogram

//    zpair = (init_target*(piplus+piminus))/(init_target * q);


    if(
//        xF_piplus>0 && xF_piminus>0 && //removed to test cutting dihadron xF
	   abs(vz_electron-vz_piplus)<20 && 
       abs(vz_electron-vz_piminus)<20){
        dihadron = piplus+piminus;
        Mdihadron = (dihadron.M());
//        theta_gh = dihadron.Angle(electron.Vect());     //Need to use dot product to find angle
//        P_h = (dihadron.P()) * sin(theta_gh);
        qv.SetX(q.Px());
        qv.SetY(q.Py());
        qv.SetZ(q.Pz()); // TVector for q momentum
        P_h = dihadron.Perp(qv);
            //Boost p and q center of mass vectors
        lv_p_gN = dihadron; // lines below replaced this in an effort to mirror The thesis code
//        lv_p_gN.SetPx(dihadron.Px()); //bad code to test a theory about thesis code's use of SetPxPyPzM (M instead of E)
//        lv_p_gN.SetPy(dihadron.Py());
//        lv_p_gN.SetPz(dihadron.Pz());
//        lv_p_gN.SetE(Mdihadron);
        lv_q_gN = q;
        lv_p_gN.Boost(gNBoostNeg);
        lv_q_gN.Boost(gNBoostNeg);
        lv_q_gN_unit.SetMagThetaPhi(1,lv_q_gN.Theta(),lv_q_gN.Phi());
        //Setting Pt_Q
        Pt_Q.SetMagThetaPhi(vecR.Dot(lv_q_gN_unit),lv_q_gN_unit.Theta(),lv_q_gN_unit.Phi());
//        x_F = (2 * ((dihadron.Vect()).Dot(q.Vect()))) / ((q.Vect()).Mag() * W); // wrong xF calculation replaced/corrected by below formula
        x_F = 2*(lv_p_gN.Vect().Dot(lv_q_gN.Vect())) / (lv_q_gN.Vect().Mag()*W); // This is actually correct maybe, seems most correct and matches best
//        x_F = xF_piplus + xF_piminus;  //Old xF calculation associated with pi+ xF and pi- xF
        // Get xF of Pi+
        xF_piplus = 2*(lv_p1_gN.Vect().Dot(lv_q_gN.Vect())) / (lv_q_gN.Vect().Mag()*W);
        // Get xF of Pi-
        xF_piminus = 2*(lv_p2_gN.Vect().Dot(lv_q_gN.Vect())) / (lv_q_gN.Vect().Mag()*W);
        lv_Mx = q + init_target - piplus - piminus;
        M_x = (lv_Mx.M()); // missing mass calculation pulled from T_Haywards Github
        
        //Setting Ph vectors
        vectPh = lv_p_gN.Vect();
        vectPhT = vectPh - Pt_Q;
        //Creating Sines
        sinPhiH = lv_e_gN.Vect().Cross(vectPhT).Dot(lv_q_gN_unit);
        //Scaling
        hScale = lv_q_gN_unit.Cross(lv_e_gN.Vect()).Mag() * lv_q_gN_unit.Cross(vectPh).Mag();
        sinPhiH = sinPhiH/hScale;
        
        //Vectors needed for Phih correction
        vectPh1 = lv_p1_gN.Vect();
        vectPhT1 = vectPh1 - Pt_Q;
        VTH1 = lv_q_gN_unit.Cross(vectPhT1);
        vT = lv_q_gN_unit.Cross(lv_e_gN_unit);
        vT.Unit();
        vTH = lv_q_gN_unit.Cross(vectPhT);
        vTH.Unit();
        //Phi_h calc from THayward
        cosphih = vT.Dot(vTH);
        
        phi_R=(_kin.phi_R(q,init_electron,piplus,piminus));
//	    phi_h=(_kin.phi_h(q,init_electron,piplus,piminus));
        phi_h = acos(cosphih);
        //if statements
        if(sinPhiH < 0.0){ phi_h = 2 * TMath::Pi() - phi_h;}
//        if(phi_R_pre < 0.0){ phi_R = 2 * TMath::Pi() + phi_R_pre;}
            th = (_kin.com_th(piplus,piminus));
	    _tree_postprocess->Fill();
	    fid++;
    }
  }
  return 0;
}
