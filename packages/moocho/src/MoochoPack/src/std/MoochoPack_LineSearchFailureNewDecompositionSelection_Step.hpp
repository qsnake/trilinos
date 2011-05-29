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

#ifndef LINE_SEARCH_FAILURE_NEW_DECOMPOSITION_SELECTION_H
#define LINE_SEARCH_FAILURE_NEW_DECOMPOSITION_SELECTION_H

#include "MoochoPack_NewDecompositionSelection_Strategy.hpp"
#include "IterationPack_AlgorithmStep.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"

namespace MoochoPack {

/** \brief Directs the selection of a new decomposition if the line search fails.
  *
  * If the delegated line search Step object throws a \c LineSearchFailure
  * exception, then this object directs the selection of a new
  * decomposition.  If the very next iteration also results in a linesearch
  * failure then we must quit.
  */
class LineSearchFailureNewDecompositionSelection_Step
  : public IterationPack::AlgorithmStep // doxygen needs full path
{
public:

  /// <<std comp>> members for LineSearch object.
  STANDARD_COMPOSITION_MEMBERS( IterationPack::AlgorithmStep, line_search_step );

  /// <<std comp>> members for Decomposition Select Strategy object.
  STANDARD_COMPOSITION_MEMBERS( NewDecompositionSelection_Strategy, new_decomp_strategy );

  /** \brief . */
  LineSearchFailureNewDecompositionSelection_Step(
    const line_search_step_ptr_t        &line_search_step
    ,const new_decomp_strategy_ptr_t    &new_decomp_strategy
    );

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

  int last_ls_failure_k_;

  // not defined and not to be called
  LineSearchFailureNewDecompositionSelection_Step();
  LineSearchFailureNewDecompositionSelection_Step(
    const LineSearchFailureNewDecompositionSelection_Step&);
  LineSearchFailureNewDecompositionSelection_Step& operator=(
    const LineSearchFailureNewDecompositionSelection_Step&);

};	// end class LineSearchFailureNewDecompositionSelection_Step

}	// end namespace MoochoPack 

#endif	// LINE_SEARCH_FAILURE_NEW_DECOMPOSITION_SELECTION_H
