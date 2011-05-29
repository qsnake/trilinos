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

#ifndef LINE_SEARCH_NLE_STEP_HPP
#define LINE_SEARCH_NLE_STEP_HPP

#include "MoochoPack_Types.hpp"
#include "IterationPack_AlgorithmStep.hpp"
#include "ConstrainedOptPack_DirectLineSearch_Strategy.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"

namespace MoochoPack {

/** \brief Delegates the line search to a <tt>DirectLineSearch_Strategy</tt> object.
 */
class LineSearchNLE_Step
  : public IterationPack::AlgorithmStep // doxygen needs full path
{
public:

  /// <<std comp>> members for direct_line_search
  STANDARD_COMPOSITION_MEMBERS(DirectLineSearch_Strategy,direct_line_search);
  /** \brief . */
  LineSearchNLE_Step(
    const direct_line_search_ptr_t& direct_line_search = Teuchos::null
    );

  /** Overridden from AlgorithmStep */
  //@{

  /** \brief . */
  bool do_step(Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
    , poss_type assoc_step_poss);
  /** \brief . */
  void print_step( const Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
    , poss_type assoc_step_poss, std::ostream& out, const std::string& leading_str ) const;

  //@}

};	// end class LineSearchNLE_Step

}	// end namespace MoochoPack 

#endif	// LINE_SEARCH_NLE_STEP_HPP
