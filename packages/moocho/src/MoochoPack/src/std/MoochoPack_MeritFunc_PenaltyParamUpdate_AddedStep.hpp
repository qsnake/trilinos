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

#ifndef MERIT_FUNC_PENALTY_PARAM_UPDATE_ADDED_STEP_H
#define MERIT_FUNC_PENALTY_PARAM_UPDATE_ADDED_STEP_H

#include "MoochoPack_Types.hpp"
#include "IterationPack_AlgorithmStep.hpp"

namespace MoochoPack {

/** \brief Base class for steps that update penalty parameters based on the
 * Lagrange multipliers lambda_k (or some approximation to them).
 *
 * This class contains methods for setting and querying values that
 * determine how the penalty parameters are updated.
 */
class MeritFunc_PenaltyParamUpdate_AddedStep
  : public IterationPack::AlgorithmStep // doxygen needs full path
{
public:

  /** @name Pure virtual methods to be overridden by subclasses */
  //@{

  /** \brief Set the smallest value a penalty parameter is allowed to be.
    */
  virtual void small_mu( value_type small_mu ) = 0;
  /** \brief . */
  virtual value_type small_mu() const = 0;
  /** \brief Set the ratio of min(mu(i))/max(mu(i)) >= min_mu_ratio.
    *
    * If there is only one penalty parameter this is ignored.
    * The default implementation is to just return 1.
    */
  virtual void min_mu_ratio( value_type min_mu_ratio )
  {}
  /** \brief . */
  virtual value_type min_mu_ratio() const
  {	return 1.0; };
  /** \brief Set set the factor for mu = (1 + mult_factor) * abs(lambda_k(i)).
    *
    * Here it is expented that mult_factor will be very small
    * (i.e. 1e-4)
    */
  virtual void mult_factor( value_type mult_factor ) = 0;
  /** \brief . */
  virtual value_type mult_factor() const = 0;
  /** \brief Set the total KKT error ( max(||rGL||inf/max(1,||Gf||inf),||c||inf) )
    * above which the penalty parameter will be allowed to be reduced.
    *
    * If the KKT error is less than near_solution the penalty parameter will
    * only be increased.  This is usually needed to prove convergence.
    */
  virtual void kkt_near_sol( value_type kkt_near_sol ) = 0;
  /** \brief . */
  virtual value_type kkt_near_sol() const = 0;

  //@}

};	// end class MeritFunc_PenaltyParamUpdate_AddedStep

}	// end namespace MoochoPack 

#endif	// MERIT_FUNC_PENALTY_PARAM_UPDATE_ADDED_STEP_H
