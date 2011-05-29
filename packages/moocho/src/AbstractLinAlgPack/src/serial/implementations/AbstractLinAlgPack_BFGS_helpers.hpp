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

#ifndef BFGS_HELPERS_H
#define BFGS_HELPERS_H

#include "AbstractLinAlgPack_Types.hpp"

namespace AbstractLinAlgPack {

/** \brief @name Functions to be used in BFGS updating.
 */
//@{

/** \brief Check that s'*y is sufficiently positive and print the result if it is not.
 *
 * @param  s           [in] DVector (size n): Secant update vector B*s=y.
 * @param  y           [in] DVector (size n): Secant update vector for B*s=y.
 * @param  sTy         [in] If sTy != NULL then *sTy must contain the value computed
 *                     from dot(s,y).  If sTy == NULL, then this value will be computed
 *                     internally.  This argument is included so that the same computation
 *                     does not have to be performed more than once by the client and
 *                     this function.
 * @param  out         [in/out] If out==NULL then no output will be printed if the
 *                     condition fails.  If out!=NULL and the function returns false
 *                     then a small amount of output will be sent to *out.
 * @param  func_name   [in] The name of the function this is being called from.
 *                     If the condition is not met then this name in included
 *                     in the printout to *out.  If func_name == NULL then this
 *                     string will obviously not be printed.
 *
 * @return If s'*y >= sqrt(mach_epsilon) * ||s||2 * ||y||2 then this function will return true.
 * Otherwise it will return false.
 */
bool BFGS_sTy_suff_p_d(
  const Vector    &s
  ,const Vector   &y
  ,const value_type     *sTy        = NULL
  ,std::ostream         *out        = NULL
  ,const char           func_name[] = NULL
  );

//@}

} // end namespace AbstractLinAlgPack

#endif // BFGS_HELPERS_H
