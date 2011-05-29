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

#ifndef VARIABLE_BOUNDS_TESTER_H
#define VARIABLE_BOUNDS_TESTER_H

#include "ConstrainedOptPack_Types.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

namespace ConstrainedOptPack {

/** \brief Tests that a set of variables are within their bounds.
  *
  \verbatim
    xL <= x <= xU
  \endverbatim
  *
  * The relative error for each comparison is
  * rel_err(i) = (xL(i)-x(i))/(1+||x||inf)
  * or rel_err(i) = (x(i)-xU(i))/(1+||x||inf).
  * If rel_err(i) >= error_tol, then the tests will be terminated immediately.
  * All of the rel_err(i) >= warning_tol will be printed.
  */
class VariableBoundsTester {
public:

  /** \brief . */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, warning_tol );

  /** \brief . */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, error_tol );

  /** \brief . */
  VariableBoundsTester(
      value_type	warning_tol		= 1e-10
    , value_type	error_tol		= 1e-5
    );

  /** \brief . */
  virtual ~VariableBoundsTester() {}

  /** \brief Check that the variables are within bounds.
   *
   * @param  print_all_warnings
   *              [in] If true, then all errors greater than warning_tol will
   *              be printed.
   * @param  xL   [in] Sparse lower bound vector (xL.size()==x.size())
   * @param  xL_name
   *              [in] The name of the vector xL (null terminated string).
   * @param  xU   [in] Sparse upper bound vector (xU.size()==x.size())
   * @param  xU_name
   *              [in] The name of the vector xU (null terminated string).
   * @param  x    [in] Variable to test that it is in bounds.
   * @param  x_name
   *              [in] The name of the vector x (null terminated string).
   *
   * @return \c true if all of the errors are greater than the error tolerances,
   * otherwise it returns \c false
   */
  virtual bool check_in_bounds(
    std::ostream* out, bool print_all_warnings, bool print_vectors
    ,const Vector& xL, const char xL_name[]
    ,const Vector& xU, const char xU_name[]
    ,const Vector& x,  const char x_name[]
    );

};	// end class VariableBoundsTester

}	// end namespace ConstrainedOptimizationPackTypes

#endif	// VARIABLE_BOUNDS_TESTER_H
