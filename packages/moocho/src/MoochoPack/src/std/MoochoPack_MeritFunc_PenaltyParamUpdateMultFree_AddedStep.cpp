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

#include <ostream>
#include <typeinfo>

#include "MoochoPack_MeritFunc_PenaltyParamUpdateMultFree_AddedStep.hpp"
#include "MoochoPack_NLPAlgoState.hpp"
#include "AbstractLinAlgPack_Vector.hpp"
#include "AbstractLinAlgPack_VectorStdOps.hpp"

namespace MoochoPack {

MeritFunc_PenaltyParamUpdateMultFree_AddedStep::MeritFunc_PenaltyParamUpdateMultFree_AddedStep(
  value_type    small_mu
  ,value_type   mult_factor
  ,value_type   kkt_near_sol
  )
  :MeritFunc_PenaltyParamUpdateGuts_AddedStep(small_mu,mult_factor,kkt_near_sol)
{}

// Overridden from MeritFunc_PenaltyParamUpdateGuts_AddedStep

bool MeritFunc_PenaltyParamUpdateMultFree_AddedStep::min_mu(
  NLPAlgoState& s, value_type* min_mu
  ) const
{
  using AbstractLinAlgPack::dot;

  IterQuantityAccess<VectorMutable>
    &Gf_iq    = s.Gf(),
    &nu_iq    = s.nu(),
    &Ypy_iq   = s.Ypy(),
    &c_iq     = s.c();
  if ( Gf_iq.updated_k(0) && nu_iq.updated_k(0) && Ypy_iq.updated_k(0) && c_iq.updated_k(0) ) {
    // min_mu = abs((Gf_k+nu_k)'*Ypy_k) / norm(c_k,1)
    const value_type
      dot_Gf_Ypy = dot( Gf_iq.get_k(0), Ypy_iq.get_k(0) ),
      dot_nu_Ypy = dot( nu_iq.get_k(0), Ypy_iq.get_k(0) ),
      nrm_c      = c_iq.get_k(0).norm_1(),
      small_num  = std::numeric_limits<value_type>::min();
    *min_mu = ::fabs( dot_Gf_Ypy + dot_nu_Ypy ) / ( nrm_c + small_num );
    return true;
  }
  return false;
}

void MeritFunc_PenaltyParamUpdateMultFree_AddedStep::print_min_mu_step(
  std::ostream& out, const std::string& L ) const
{
  out
    << L << "if Gf_k, nu_k, Ypy_k and c_k are updated then\n"
    << L << "   min_mu = abs((Gf_k+nu_k)'*Ypy_k) / ( norm(c_k,1) + small_num )\n"
    << L << "   update_mu = true\n"
    << L << "else\n"
    << L << "   update_mu = false\n"
    << L << "endif\n"
    ;
}

}	// end namespace MoochoPack
