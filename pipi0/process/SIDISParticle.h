#ifndef SIDISParticle_h
#define SIDISParticle_h

class SIDISParticle:
{
 public: 
  SIDISParticle(){}
  virtual ~SIDISParticle() {}
  
  // How to add new PROPERTY tag:
  // 1. Add new tag below with a unique value
  // 2. Go to SIDISParticle.C and write a short name to SIDISParticle::get_property_info
  enum PROPERTY
  {//
    
    // -- 1-10   Generic Information -- //
    part_nParticle   = 1;
    part_nPhoton   = 1;
    // -- 10-100 Particle Kinematics and Identification -- //
    part_pid = 10;
    part_pt  = 11;
    part_pz  = 12;
    part_E   = 13;
    
  };

  enum PROPERTY_TYPE
  {//
    type_int = 1;
    type_double = 2;
  };

  static std::pair<const std::string,PROPERTY_TYPE> get_property_info(PROPERTY prop_id);
};
#endif
