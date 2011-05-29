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

#ifndef MERIT_FUNC_NLP_DIREC_DERIV_H
#define MERIT_FUNC_NLP_DIREC_DERIV_H

#include "ConstrainedOptPack_Types.hpp"

namespace ConstrainedOptPack {

/** \brief This class provides a mix-in interface for allowing subclass merit
 * functions to compute the directional 1D derivative at a base point.
 *
 * The quantities Gf(xo) (gradient of f(xo))
 * c(xo), h(xo) and d are used by several
 * types of merit functions to calculate the derivative of:<br>
 * d(phi(x_k + alpha_k*d_k))/d(alpha_k) at alpha_k = 0.
 *
 * It is generally assumed that d satisfies Gc_k'*d_k + c_k = 0 otherwise the
 * merit function would need Gc_k to compute this directional derivative
 * properly.
 */
class MeritFuncNLPDirecDeriv {
public:

  /** \brief . */
  virtual ~MeritFuncNLPDirecDeriv() {}

  /** @name To be overridden by subclasses */
  //@{

  /** \brief Calculate d(phi(x_k + alpha_k*d_k))/d(alpha_k) at alpha_k = 0.
    *
    * The value is stored internally by the subclass are returned by its
    * deriv() member usually.  The value is also returned from this
    * function.
    *
    * If the sizes of the vectors input do not aggree then
    * #std::length_error# exception will be thrown.
    */
  virtual value_type calc_deriv(
    const Vector    &Gf_k
    ,const Vector   *c_k
    ,const Vector   *h_k
    ,const Vector   *hl
    ,const Vector   *hu
    ,const Vector   &d_k
    ) = 0;

  //@}

};	// end class MeritFuncNLPDirecDeriv

}	// end namespace ConstrainedOptPack

#endif	// MERIT_FUNC_NLP_DIREC_DERIV_H
