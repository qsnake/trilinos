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

#ifndef CALC_REDUCED_GRAD_LAGRANGIAN_STD_ADDED_STEP_H
#define CALC_REDUCED_GRAD_LAGRANGIAN_STD_ADDED_STEP_H

#include "MoochoPack_Types.hpp"
#include "IterationPack_AlgorithmStep.hpp"

namespace MoochoPack {

/** \brief Calculates the reduced gradient of the Lagrangian
 * <tt>rGL = rGf + Z' * nu + GcUP' * lambda(equ_undecomp) + GhUP' * lambdaI(inequ_undecomp)</tt> 
 */
class CalcReducedGradLagrangianStd_AddedStep
  : public IterationPack::AlgorithmStep // doxygen needs full path
{
public:

  // ////////////////////
  // Overridden

  /** \brief . */
  bool do_step(Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
    , poss_type assoc_step_poss);

  /** \brief . */
  void print_step( const Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
    , poss_type assoc_step_poss, std::ostream& out, const std::string& leading_str ) const;

};	// end class CalcReducedGradLagrangianStd_AddedStep

}	// end namespace MoochoPack 

#endif	// CALC_REDUCED_GRAD_LAGRANGIAN_STD_ADDED_STEP_H
