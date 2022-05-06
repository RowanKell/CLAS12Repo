#include "SIDISParticle.h"


std::pair<const std::string, SIDISParticle::PROPERTY_TYPE> SIDISParticle::get_property_info(const PROPERTY prop_id)
{
  switch(prop_id)
    {
      
      // ---------------

    case part_pid:
      return std::make_pair("pid", SIDISParticle::type_int);
    case part_px:
      return std::make_pair("px", SIDISParticle::type_float);
    case part_py:
      return std::make_pair("py", SIDISParticle::type_float);
    case part_pz:
      return std::make_pair("pz", SIDISParticle::type_float);
    case part_pt:
      return std::make_pair("pt", SIDISParticle::type_float);
    case part_p:
      return std::make_pair("p", SIDISParticle::type_float);
    case part_E:
      return std::make_pair("E", SIDISParticle::type_float);
    case part_theta:
      return std::make_pair("theta", SIDISParticle::type_float);
    case part_eta:
      return std::make_pair("eta", SIDISParticle::type_float);
    case part_phi:
      return std::make_pair("phi", SIDISParticle::type_float);
      
    case part_pindex:
      return std::make_pair("pindex", SIDISParticle::type_int);
    case part_beta:
      return std::make_pair("beta", SIDISParticle::type_float);
    case part_chi2:
      return std::make_pair("chi2", SIDISParticle::type_float);

      // ----------------
      

    default:
      std::cout << "SIDISParticle::get_property_info - Fatal Error - unknown prop_id " << prop_id << std::endl;
      return std::make_pair("-1", SIDISParticle::type_int);
    }
}


