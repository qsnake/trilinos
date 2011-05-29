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

#ifndef MERIT_FUNC_CALC_1D_H
#define MERIT_FUNC_CALC_1D_H

#include "ConstrainedOptPack_Types.hpp"

namespace ConstrainedOptPack {

/** \brief Abstracts a 1D merit function {abstract}.
  *
  * This is the interface that line search algorithms use to compute
  * the value of the merit function at alpha (phi(alpha)) and
  * to retrieve the initial descent derivative of the merit function
  * (using \c deriv()).
  */
class MeritFuncCalc1D {
public:

  /** \brief . */
  virtual ~MeritFuncCalc1D() {}

  /// Return the value of the merit function at alpha.
  virtual value_type operator()( value_type alpha ) const = 0;

  /// Return the derivative of the merit function at alpha = 0
  virtual value_type deriv() const = 0;

  /// Print the particular merit function
  virtual void print_merit_func(
    std::ostream& out, const std::string& leading_str ) const = 0;

};	// end class MeritFuncCalc1D

}	// end namespace ConstrainedOptPack

#endif	// MERIT_FUNC_CALC_1D_H
