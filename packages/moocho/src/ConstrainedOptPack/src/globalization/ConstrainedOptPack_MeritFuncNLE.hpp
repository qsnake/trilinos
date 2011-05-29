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

#ifndef MERIT_FUNC_NLE_H
#define MERIT_FUNC_NLE_H

#include <iosfwd>

#include "ConstrainedOptPack_Types.hpp"

namespace ConstrainedOptPack {

/** \brief Base class for all merit functions for systems of NonLinear Equations (NLE) {abstract}.
  */
class MeritFuncNLE {
public:

  /** \brief . */
  class InvalidInitialization : public std::logic_error
  {public: InvalidInitialization(const std::string& what_arg) : std::logic_error(what_arg) {}};

  /** \brief . */
  virtual ~MeritFuncNLE() {}

  /** \brief Return the value of the merit function at <tt>c(x)</tt>.
   * This interface requires the client to compute
   * <tt>c(x)</tt> and pass it to this function to have
   * the value of phi(c) calculated.
   *
   * If the merit function has not been initialized properly
   * then a <tt>InvalidInitialization</tt> exception will be thrown.
   */
  virtual value_type value(const Vector& c) const = 0;

  /** \brief Return the value of the directional derivative of the 
   * merit function w.r.t. <tt>alpha</tt> at <tt>alpha = 0</tt>.  In other words
   * compute return <tt>d(phi(c(x)))/d(alpha_k)</tt> at <tt>alpha_k = 0</tt>
   * where <tt>x = x_k + alpha_k * d_k</tt>.
   *
   * If the merit function has not been initialized properly
   * then a <tt>InvalidInitialization</tt> exception will be thrown.
   */
  virtual value_type deriv() const = 0;

  /** \brief Print the merit funciton
    */
  virtual void print_merit_func(std::ostream& out
    , const std::string& leading_str) const = 0;

};	// end class MeritFuncNLE

}	// end namespace ConstrainedOptPack

#endif	// MERIT_FUNC_NLE_H
