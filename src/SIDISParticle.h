#ifndef SIDISParticle_h
#define SIDISParticle_h

#include <utility>
#include <iostream>
#include <climits>
class SIDISParticle
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
    part_nParticle   = 1,
    part_nPhoton   = 2,
    // -- 10-100 Particle Kinematics and Identification -- //
    part_pid = 10,
    part_pt  = 11,
    part_pz  = 12,
    part_E   = 13
    
  };

  enum PROPERTY_TYPE
  {//
    type_int = 1,
    type_float = 2
  };

  static std::pair<const std::string,PROPERTY_TYPE> get_property_info(PROPERTY prop_id);
  virtual float get_property_float(const PROPERTY prop_id) const {return -999.999;}
  virtual int   get_property_int(const PROPERTY prop_id) const {return INT_MIN;}
  virtual unsigned int   get_property_uint(const PROPERTY prop_id) const {return UINT_MAX;}
};
#endif
