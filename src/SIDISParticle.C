#include "SIDISParticle.h"


std::pair<const std::string, SIDISParticle::PROPERTY_TYPE> SIDISParticle::get_property_info(const PROPERTY prop_id)
{
  switch(prop_id)
    {
      
      // ---------------

    case part_pid:
      return std::make_pair("pid", SIDISParticle::type_int);
    case part_pt:
      return std::make_pair("pt", SIDISParticle::type_float);
    case part_pz:
      return std::make_pair("pz", SIDISParticle::type_float);
    case part_E:
      return std::make_pair("E", SIDISParticle::type_float);

      // ----------------
      
    default:
      std::cout << "SIDISParticle::get_property_info - Fatal Error - unknown prop_id " << prop_id << std::endl;
      return std::make_pair("-1", SIDISParticle::type_int);
    }
}


