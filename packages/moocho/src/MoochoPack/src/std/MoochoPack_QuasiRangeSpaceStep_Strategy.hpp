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

#ifndef QUASI_RANGE_SPACE_STEP_STRATEGY_H
#define QUASI_RANGE_SPACE_STEP_STRATEGY_H

#include "MoochoPack_Types.hpp"

namespace MoochoPack {

/** \brief Abstract interface for a strategy object that will compute a step that will
 * approximalty solve a range space subproblem {abstract}.
 */
class QuasiRangeSpaceStep_Strategy {
public:

  /** \brief . */
  virtual ~QuasiRangeSpaceStep_Strategy() {}

  /** \brief Compute a step that will approximatly solve a range-space subproblem.
   *
   * This function will compute a step <tt>v</tt> that will approximatly satisfy:
   *
   * <tt>||Gc_k'*v + c(xo)|| < ||c(xo)||</tt>
   *
   * The above norm ||.|| could be any valid norm and the implementation is free to
   * define what descent means any way it would like.  It is assumed that this
   * step will be computed by using the <tt>Gc_k</tt> but other implementations are
   * possible.  Any information being used in the algorithm can be used
   * to compute this step in a reasonable way.  Note that the 
   * inequalities do not have to (and should not in most cases) be considered
   * in this computation.  Note that whatever means is used to compute <tt>v</tt>
   * that it better give a descent direction for <tt>||c(x)||</tt> but there is no guarantee
   * for this if <tt>||xo - x_k||</tt> is large since <tt>Gc_k</tt> may not accurately approximate
   * <tt>Gc(xo)</tt>.
   *
   * @param out      [out] Output stream journal data is written to.
   * @param olevel   [in] Output level for printing to <tt>out</tt>.
   * @param algo     [in/out] The NLPAlgo object.  This object can be queryed for information.
   * @param s        [in/out] NLPAlgoState object.  May be queried or modified if needed.
   * @param xo       [in] Base point vector (size n) xo.
   * @param c_xo     [out] Constraints residual c(xo).
   * @param v        [out] Computed step vector (size n).  Must not be NULL.
   *
   * @return Returns true if a step could be found and false otherwise.
   */
   virtual bool solve_quasi_range_space_step(
    std::ostream& out, EJournalOutputLevel olevel, NLPAlgo *algo, NLPAlgoState *s
    ,const Vector& xo, const Vector& c_xo, VectorMutable* v
      ) = 0;

  /** \brief This function will print a description of the computations and logic used.
   */
  virtual void print_step( std::ostream& out, const std::string& leading_str ) const = 0;

}; // end class QuasiRangeSpaceStep_Strategy

} // end namespace MoochoPack

#endif // QUASI_RANGE_SPACE_STEP_STRATEGY_H
