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

#ifndef MERIT_FUNC_NLE_SQR_RESID_H
#define MERIT_FUNC_NLE_SQR_RESID_H

#include "ConstrainedOptPack_MeritFuncNLE.hpp"

namespace ConstrainedOptPack {

/** \brief A merit function for the square of the constriant values.
  *
  * phi(x) = 1/2 * c(x)'*c(x)
  *
  * Dphi(x_k,d_k) = - c(x)'*c(x)
  *
  * Note that the definition of Dphi(x_k,d_k) assumes
  * that Gc_k'*d_k + c_k = 0.  In otherwords, d_k must
  * satisfiy the linearized equality constraints at
  * at x_k.
  *
  * Implicit copy constructor and assignment operators
  * are allowed.
  */
class MeritFuncNLESqrResid : public MeritFuncNLE {
public:

  /// Initializes deriv() = 0
  MeritFuncNLESqrResid();

  /** \brief . */
  value_type calc_deriv( const Vector& c_k );

  // ////////////////////////////////
  // Overridden from MeritFuncNLE

  /** \brief . */
  value_type value(const Vector& c) const;

  /** \brief . */
  value_type deriv() const;

  /** \brief . */
  void print_merit_func(std::ostream& out
    , const std::string& leading_str) const;

private:
  value_type deriv_;

};	// end class MeritFuncNLESqrResid

}	// end namespace ConstrainedOptPack

#endif	// MERIT_FUNC_NLE_SQR_RESID_H
