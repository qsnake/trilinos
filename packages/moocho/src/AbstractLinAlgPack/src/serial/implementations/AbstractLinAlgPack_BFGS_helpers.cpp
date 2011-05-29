// @HEADER
// ***********************************************************************
// 
// Moocho: Multi-functional Object-Oriented arCHitecture for Optimization
//                  Copyright (2003) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#include <math.h>
#include <limits>
#include <ostream>

#include "AbstractLinAlgPack_BFGS_helpers.hpp"
#include "AbstractLinAlgPack_Vector.hpp"
#include "AbstractLinAlgPack_VectorStdOps.hpp"

bool AbstractLinAlgPack::BFGS_sTy_suff_p_d(
  const Vector    &s
  ,const Vector   &y
  ,const value_type     *sTy_in
  ,std::ostream         *out
  ,const char           func_name[]
  )
{
  const value_type
    sTy          = sTy_in ? *sTy_in : AbstractLinAlgPack::dot(s,y),
    nrm_s        = s.norm_2(),
    nrm_y        = y.norm_2(),
    sqrt_macheps = ::sqrt(std::numeric_limits<value_type>::epsilon()),
    min_sTy      = sqrt_macheps * nrm_s * nrm_y;
  // Skip update if: s'*y < sqrt(macheps)*||s||2*||y||2 (Dennis and Schnabel, A9.4.2)
  const bool
    sufficiently_p_d = sTy > min_sTy;
  if( !sufficiently_p_d && out ) {
    if(func_name)
      *out << func_name << " : ";
    *out
      << "Error, s'*y = " << sTy << " < sqrt(mach_eps) * ||s||2 * ||y||2 = "
      << sqrt_macheps << " * " << nrm_s << " * " << nrm_y << " = " << min_sTy
      << "\nTherefore the BFGS update is illdefined!\n";
  }
  return sufficiently_p_d;
}
