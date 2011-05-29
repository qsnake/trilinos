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

#ifndef MERIT_FUNC_PENALTY_PARAM_UPDATE_MULT_FREE_ADDED_STEP_H
#define MERIT_FUNC_PENALTY_PARAM_UPDATE_MULT_FREE_ADDED_STEP_H

#include "MoochoPack_MeritFunc_PenaltyParamUpdateGuts_AddedStep.hpp"

namespace MoochoPack {

/** \brief Specializes the update of the penalty parameter for a merit function as:
  * min_mu = |(Gf_k+nu_k)'* Ypy_k| / ||c_k||1.
  */
class MeritFunc_PenaltyParamUpdateMultFree_AddedStep
  : public MeritFunc_PenaltyParamUpdateGuts_AddedStep
{
public:

  /** \brief . */
  MeritFunc_PenaltyParamUpdateMultFree_AddedStep(
    value_type   small_mu = 1e-6
    ,value_type  mult_factor = 1e-4
    ,value_type  kkt_near_sol = 1.0
    );

protected:

  /** @name Overridden from MeritFunc_PenaltyParamUpdateGuts_AddedStep */
  //@{

  /** \brief . */
  bool min_mu( NLPAlgoState& s, value_type* min_mu ) const;
  /** \brief . */
  void print_min_mu_step( std::ostream& out
    , const std::string& leading_str ) const;

  //@}
  
};	// end class MeritFunc_PenaltyParamUpdateMultFree_AddedStep

}	// end namespace MoochoPack 

#endif	// MERIT_FUNC_PENALTY_PARAM_UPDATE_MULT_FREE_ADDED_STEP_H
