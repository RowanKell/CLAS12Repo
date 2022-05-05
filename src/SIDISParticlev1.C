#include "SIDISParticlev1.h"
#include <climits>

using namespace std;

SIDISParticlev1::SIDISParticlev1():
  _candidateid(INT_MAX)
{
}

void
SIDISParticlev1::Reset()
{
  prop_map.clear();
}

float
SIDISParticlev1::get_property_float(const PROPERTY prop_id) const
{
  prop_map_t::const_iterator i = prop_map.find(prop_id);
  if (i!=prop_map.end()) return u_property(i->second).ddata;

  return   -999.999 ;
}

int
SIDISParticlev1::get_property_int(const PROPERTY prop_id) const
{
  prop_map_t::const_iterator i = prop_map.find(prop_id);
  if (i!=prop_map.end()) return u_property(i->second).idata;

  return   INT_MIN ;
}

uint
SIDISParticlev1::get_property_uint(const PROPERTY prop_id) const
{
  prop_map_t::const_iterator i = prop_map.find(prop_id);
  if (i!=prop_map.end()) return u_property(i->second).uidata;

  return   UINT_MAX ;
}

void
SIDISParticlev1::set_property(const PROPERTY prop_id, const float value)
{
  prop_map[prop_id] = u_property(value).get_storage();
}

void
SIDISParticlev1::set_property(const PROPERTY prop_id, const int value)
{
  prop_map[prop_id] = u_property(value).get_storage();
}

void
SIDISParticlev1::set_property(const PROPERTY prop_id, const unsigned int value)
{
  prop_map[prop_id] = u_property(value).get_storage();
}

unsigned int
SIDISParticlev1::get_property_nocheck(const PROPERTY prop_id) const
{
  prop_map_t::const_iterator iter = prop_map.find(prop_id);
  if (iter != prop_map.end())
    {
      return iter->second;
    }
  return UINT_MAX;
}
