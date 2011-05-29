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

#ifndef FEASIBILITY_STEP_STRATEGY_H
#define FEASIBILITY_STEP_STRATEGY_H

#include "MoochoPack_Types.hpp"

namespace MoochoPack {

/** \brief Abstract interface for a strategy object that will compute a step that will
 * improve feasibility (at least descent) {abstract}.
 */
class FeasibilityStep_Strategy {
public:

  /** \brief . */
  virtual ~FeasibilityStep_Strategy() {}

  /** \brief Compute a step that improves feasibility (at least locally).
   *
   * This function will compute a step <tt>w</tt> that satisfies:
   *
   * <tt>d_bounds_k.l <= xo + w - x_k <= d_bounds_k.u</tt>
   *
   * and gives descent for <tt>||c(xo + beta*w)||</tt> for at least small <tt>beta > 0</tt>.
   * This norm <tt>||.||</tt> could be any valid norm and the implementation is free to
   * define what descent means any way it would like.  Any information being used
   * in the algorithm can be used to compute this step.
   *
   * @param out     [out] Output stream journal data is written to.
   * @param olevel  [in] Output level for printing to #out#
   * @param algo    [in/out] The NLPAlgo object.  This object can be queryed for
   *                information.
   * @param s       [in/out] NLPAlgoState object.  May be queried or modified if needed.
   * @param xo      [in] Base point vector (size n) xo.
   * @param c_xo    [in] c(xo).
   * @param w       [out] Computed step vector (size n) w.  Must not be NULL.
   *
   * @return Returns true if a step that reduces feasibility subject to the bounds could
   * be found and false otherwise.
   */
   virtual bool compute_feasibility_step(
    std::ostream& out, EJournalOutputLevel olevel, NLPAlgo *algo, NLPAlgoState *s
    ,const Vector& xo, const Vector& c_xo, VectorMutable* w
      ) = 0;

  /** \brief This function will print a description of the computations and logic used.
   */
  virtual void print_step( std::ostream& out, const std::string& leading_str ) const = 0;

}; // end class FeasibilityStep_Strategy

} // end namespace MoochoPack

#endif // FEASIBILITY_STEP_STRATEGY_H
