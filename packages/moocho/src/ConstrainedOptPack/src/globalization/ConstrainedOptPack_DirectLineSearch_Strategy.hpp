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

#ifndef DIRECT_LINE_SEARCH_STRATEGY_H
#define DIRECT_LINE_SEARCH_STRATEGY_H

#include <stdexcept>
#include <iosfwd>

#include "ConstrainedOptPack_Types.hpp"

namespace ConstrainedOptPack {

/** \brief Abstract strategy interface for 1D line searches {abstract}.
 *
 * This is the interface for strategy objects that perform a line search
 * from an initial point along a search direction given a merit function.
 */
class DirectLineSearch_Strategy {
public:

  /** @name Exceptions */
  //@{

  /// Thrown if the direction vector d_k is not a descent direction for the merit funciton
  class NotDescentDirection : public std::logic_error
  {public: NotDescentDirection(const std::string& what_arg) : std::logic_error(what_arg) {}};

  //@}

  /** \brief . */
  virtual ~DirectLineSearch_Strategy() {}

  /// Set the maximum number of iterations
  virtual void set_max_iter(int max_iter) = 0;

  /// Get the maximum number of iterations
  virtual int max_iter() const = 0;

  /// Get the number of iterations performed
  virtual int num_iterations() const = 0;

  /** \brief Called to perform the linesearch.
   *
   * This operaion  computes the
   * approximate minimum to a merit function along a search direcation.
   * More specifically the following problem is approximatly solved:<br>
   *
   * min  phi(alpha)  s.t. alpha = [0, alpha_upper]<br>
   *
   * Actually, if the initial alpha satisfys an internal descent requirement, then
   * it will be choosen over smaller values of alpha that may result in a 
   * greater reduction in the given merit funciton.
   *
   * If the maximum number of iterations is exceeded then the subclass will return
   * false and will return the values of alpha_k, x_kp1, and phi_kp1 for the 
   * lowest value of phi_kp1 found, and the last call to phi.value(x) will
   * be this best x_kp1.
   *
   * Preconditions: \begin{itemize}
   * \item <tt>phi.deriv(d_k) < 0</tt> (throw NotDescentDirection)
   * \end{itemize}
   *
   * @param  phi    [in] The merit function object that will compute <tt>phi.value(alpha)</tt>
   *                and the descent derivative.
   * @param  phi_k  [in] The value of <tt>phi.value(0)</tt>.  Not computed internally
   *                for the sake of efficency.
   * @param  alpha_k
   *                [in/out] The initial <tt>alpha_k</tt> to try on input (usually 1).
   *                On output <tt>alpha_k</tt> is the accepted value for a successful
   *                line search, or it will be the alpha_k for the minimum phi
   *                found for a line search failure.
   * @param  phi_kp1
   *                [in/out] Merit function at new point.
   *                On input it must be the value of <tt>phi.value(alpha_k)</tt>
   *                and on output is set to <tt>phi.value(alpha_k)</tt>.
   * @param  out    [in/out] If != 0 then output is sent to this stream to record
   *                the progress of the linesearch iterations.  The default
   *                is zero.
   *
   * @return \c true: Successful line search; \c false: Line search failure.
   */
  virtual bool do_line_search(
    const MeritFuncCalc1D   &phi
    ,value_type             phi_k
    ,value_type             *alpha_k
    ,value_type             *phi_kp1
    ,std::ostream           *out      = 0
    ) = 0;

  /** \brief Print the direct line search algorithm.
    *
    * The default does nothing.
    */
  virtual void print_algorithm(std::ostream& out, const std::string& leading_str) const
  {}

};	// end class DirectLineSearch_Strategy

}	// end namespace ConstrainedOptPack

#endif	// DIRECT_LINE_SEARCH_STRATEGY_H
