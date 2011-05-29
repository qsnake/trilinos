#if 0

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

#include <ostream>
#include <typeinfo>

#include "MoochoPack_MeritFunc_PenaltyParamUpdateWithMult_AddedStep.hpp"
#include "MoochoPack_NLPAlgoState.hpp"
#include "ConstrainedOptPack/src/VectorWithNorms.h"

namespace MoochoPack {

MeritFunc_PenaltyParamUpdateWithMult_AddedStep::MeritFunc_PenaltyParamUpdateWithMult_AddedStep(
      const merit_func_ptr_t& merit_func, value_type small_mu
    , value_type mult_factor, value_type kkt_near_sol )
  : MeritFunc_PenaltyParamUpdateGuts_AddedStep(merit_func,small_mu,mult_factor,kkt_near_sol)
{}

// Overridden from MeritFunc_PenaltyParamUpdateGuts_AddedStep

bool MeritFunc_PenaltyParamUpdateWithMult_AddedStep::min_mu(
  NLPAlgoState& s, value_type* min_mu ) const
{
  if ( s.lambda().updated_k(0) ) {
    *min_mu = s.lambda().get_k(0).norm_inf();
    return true;
  }
  return false;
}

void MeritFunc_PenaltyParamUpdateWithMult_AddedStep::print_min_mu_step(
  std::ostream& out, const std::string& L ) const
{
  out
    << L << "if lambda_k is updated then\n"
    << L << "    min_mu = norm( lambda_k, inf )\n"
    << L << "    update_mu = true\n"
    << L << "else\n"
    << L << "    update_mu = false\n"
    << L << "endif\n"
    ;
}

}	// end namespace MoochoPack

#endif // 0
