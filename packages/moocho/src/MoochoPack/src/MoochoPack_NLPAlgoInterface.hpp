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

#ifndef RSQP_ALGO_INTERFACE_H
#define RSQP_ALGO_INTERFACE_H

#include "MoochoPack_Types.hpp"
#include "MoochoPack_NLPSolverClientInterface.hpp"

namespace MoochoPack {

/** \brief Interface \c NLPAlgoContainer uses to access \c NLPAlgo.
 *
 * This interface helps avoid dangerous usage stategies for an \c NLPAlgo object.
 */
class NLPAlgoInterface {
public:

  /** \brief . */
  virtual ~NLPAlgoInterface() {}

  /** \brief Print the algorithm description.
   */
  virtual void interface_print_algorithm(std::ostream& out) const = 0;

  /** \brief Start the iterations.
    *
    * This function returns true if the solution was found and false
    * if the maximum number of iterations was reached before the
    * solution was found.
    */
  virtual NLPSolverClientInterface::EFindMinReturn dispatch() = 0;

  /** \brief Return the state object.
   */
  virtual const NLPAlgoState& retrieve_state() const = 0;

  /** @name Algorithm timing */
  //@{

  /** \brief . */
  virtual void interface_set_algo_timing( bool algo_timing ) = 0;
  /** \brief . */
  virtual bool interface_algo_timing() const = 0;
  /** \brief . */
  virtual void interface_print_algorithm_times( std::ostream& out ) const = 0;

  //@}

};	// end class NLPAlgoInterface

}	// end namespace MoochoPack

#endif	// RSQP_ALGO_INTERFACE_H
