#include "SundanceCoordinateSystem.hpp"

using namespace Sundance;

CoordinateSystem 
CoordinateSystemBuilder::makeCoordinateSystem(const std::string& name)
{
  RCP<CoordinateSystemBase> rtn;

  if (name=="Cartesian")
  {
    rtn = rcp(new CartesianCoordinateSystem());
  }
  else if (name=="Meridional Cylindrical")
  {
    rtn = rcp(new MeridionalCylindricalCoordinateSystem());
  }
  else if (name=="Radial Spherical")
  {
    rtn = rcp(new RadialSphericalCoordinateSystem());
  }
  else
  {
    TEST_FOR_EXCEPTION(true, RuntimeError,
      "coordinate system type=[" << name << "] not recognized");
  }
  
  return rtn;
}


