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

#ifndef CHECK_DECOMPOSITION_FROM_PY_STEP_H
#define CHECK_DECOMPOSITION_FROM_PY_STEP_H

#include "MoochoPack_NewDecompositionSelection_Strategy.hpp"
#include "IterationPack_AlgorithmStep.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

namespace MoochoPack {

/** \brief Check if the decomposition is going singular and if it is
 * select a new decomposition.
 *
 * This steps checks if the decomposition is going singular by looking
 * to see if the computation for the range space step looks bad.
 *
 * In particular we want to know how \c cond(R) is changing.  We know that:
 \verbatim
 py = -inv(R) * c
 --> ||py|| <= ||inv(R)|| * ||c||
 --> ||py|| / ||c|| <= ||inv(R)||
 --> ( ||py|| / ||c|| ) * ||R|| <= ||inv(R)|| * ||R|| = cond(R)
 \endverbatim
 *
 * If we assume <tt>||R|| > 1</tt> we know that <tt>cond(R) > beta = ||py||/||c||</tt>
 * so if \c beta is very large then cond(R) is even larger.  Since the best decomposition
 * we can find may be fairly illconditioned we may not want to use and absolute measure
 * of \c beta to determine if the decomposition is illconditioned.  Instead we will look
 * at the change in \c beta between iterations and if \c beta gets very much
 * larger than its minimum value (i.e. <tt>beta / beta_min > max_cond_change</tt> 
 * ( default = 10000 ) ) then we will select a new decomposition.  Also if
 * <tt>beta > max_cond</tt> ( default = <tt>0.01 / mach_eps</tt> ) then we know we will be
 * computing inacurate solutions so we must select a new decomposition.
 */
class CheckDecompositionFromPy_Step
  : public IterationPack::AlgorithmStep // doxygen needs full name
{
public:

  /** \brief . */
  STANDARD_COMPOSITION_MEMBERS( NewDecompositionSelection_Strategy, new_decomp_strategy );
  /** \brief . */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, max_decomposition_cond_change_frac );
  /** \brief . */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, max_cond );

  /** \brief . */
  CheckDecompositionFromPy_Step(
    const new_decomp_strategy_ptr_t   &new_decomp_strategy
    ,value_type                       max_decomposition_cond_change_frac = 100.0
    );

  /// Call the reset initialization of all defaults.
  void reset();

  /** @name Overridden from AlgorithmStep */
  //@{
  /** \brief . */
  bool do_step(Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
    , poss_type assoc_step_poss);
  /** \brief . */
  void print_step( const Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
    , poss_type assoc_step_poss, std::ostream& out, const std::string& leading_str ) const;
  //@}

private:

  value_type  beta_min_;

  // Not defined and not to be called
  CheckDecompositionFromPy_Step();

}; // end class CheckDecompositionFromPy_Step

} // end namespace MoochoPack 

#endif // CHECK_DECOMPOSITION_FROM_PY_STEP_H
