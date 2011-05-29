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

#ifndef NEW_DECOMPOSITION_SELECTION_STD_STRATEGY_H
#define NEW_DECOMPOSITION_SELECTION_STD_STRATEGY_H

#include "MoochoPack_NewDecompositionSelection_Strategy.hpp"
#include "MoochoPack_DecompositionSystemHandlerSelectNew_Strategy.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"

namespace MoochoPack {

/** \brief Just force the decomposition system object to select a new
 * decomposition and let everyone else fend for themselves.
 */
class NewDecompositionSelectionStd_Strategy
  : public NewDecompositionSelection_Strategy
{
public:

  /// «std comp» members for range/null decomposition handler
  STANDARD_COMPOSITION_MEMBERS( DecompositionSystemHandlerSelectNew_Strategy, decomp_sys_handler );

  /** \brief . */
  NewDecompositionSelectionStd_Strategy(
    const decomp_sys_handler_ptr_t   &decomp_sys_handler
    );

  /** @name Overridden from NewDecompositionSelection_Strategy */
  //@{

  bool new_decomposition(
    NLPAlgo& algo, Algorithm::poss_type step_poss
    ,IterationPack::EDoStepType type, Algorithm::poss_type assoc_step_poss
    );
  /** \brief . */
  void print_new_decomposition(
    const NLPAlgo& algo, Algorithm::poss_type step_poss
    ,IterationPack::EDoStepType type, Algorithm::poss_type assoc_step_poss
    ,std::ostream& out, const std::string& leading_str
    ) const;

  //@}

private:

  // Not defined and not to be called
  NewDecompositionSelectionStd_Strategy();

};	// end class NewDecompositionSelectionStd_Strategy

}	// end namespace MoochoPack 

#endif	// NEW_DECOMPOSITION_SELECTION_STD_STRATEGY_H
