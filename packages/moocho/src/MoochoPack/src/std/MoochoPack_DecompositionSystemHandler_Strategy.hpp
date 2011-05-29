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

#ifndef DECOMPOSITION_SYSTEM_HANDLER_STRATEGY_H
#define DECOMPOSITION_SYSTEM_HANDLER_STRATEGY_H

#include "MoochoPack_Types.hpp"
#include "IterationPack_Algorithm.hpp"

namespace MoochoPack {

/** \brief Interface for range/null decomposition handling.
 *
 * ToDo: Finish documentation!
 */
class DecompositionSystemHandler_Strategy {
public:
  
  /** @name Public types */
  //@{

  /** \brief . */
  enum EDecompSysTesting { DST_DEFAULT, DST_TEST, DST_NO_TEST };
  /** \brief . */
  enum EDecompSysPrintLevel { DSPL_USE_GLOBAL, DSPL_LEAVE_DEFAULT };

  //@}

  /** \brief . */
  virtual ~DecompositionSystemHandler_Strategy() {}

  /** \brief Update the decomposition.
   *
   * This method may select a new decomposition (permuting the variables
   * and constriants) and/or take control of the algorithm.
   */
  virtual bool update_decomposition(
    NLPAlgo                                &algo
    ,NLPAlgoState                          &s
    ,NLPFirstOrder                         &nlp
    ,EDecompSysTesting                     decomp_sys_testing
    ,EDecompSysPrintLevel                  decomp_sys_testing_print_level
    ,bool                                  *new_decomp_selected
    ) = 0;

  /** \brief Print the algorithm used for updating the decomposition.
   */
  virtual void print_update_decomposition(
    const NLPAlgo                          &algo
    ,const NLPAlgoState                    &s
    ,std::ostream                          &out
    ,const std::string                     &leading_spaces
    ) const = 0;

}; // end class DecompositionSystemHandler_Strategy

} // end namespace MoochoPack

#endif // DECOMPOSITION_SYSTEM_HANDLER_STRATEGY_H
