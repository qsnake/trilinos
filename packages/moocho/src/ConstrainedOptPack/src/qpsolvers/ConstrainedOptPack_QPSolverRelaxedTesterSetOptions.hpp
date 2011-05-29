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

#ifndef QP_SOLVER_RELAXED_TESTER_SET_OPTIONS_H
#define QP_SOLVER_RELAXED_TESTER_SET_OPTIONS_H

#include "ConstrainedOptPack_QPSolverRelaxedTester.hpp"
#include "OptionsFromStreamPack_SetOptionsFromStreamNode.hpp"
#include "OptionsFromStreamPack_SetOptionsToTargetBase.hpp"

namespace ConstrainedOptPack {

/** \brief Set options for QPSolverRelaxedTester from an
 * OptionsFromStream object.
 *
 * The default options group name is QPSolverRelaxedTester.
 *
 * The options group is:
 \verbatim
  options_group QPSolverRelaxedTester {
      opt_warning_tol   = 1e-10;  *** Tolerances for optimality conditions
      opt_error_tol     = 1e-5;
      feas_warning_tol  = 1e-10;  *** Tolerances for feasibility
      feas_error_tol    = 1e-5;
      comp_warning_tol  = 1e-10;  *** Tolerances for complementarity
      comp_error_tol    = 1e-5;
  }
  \endverbatim
  */
class QPSolverRelaxedTesterSetOptions
  : public OptionsFromStreamPack::SetOptionsFromStreamNode 
    , public OptionsFromStreamPack::SetOptionsToTargetBase<
      QPSolverRelaxedTester >
{
public:

  /** \brief . */
  QPSolverRelaxedTesterSetOptions(
      QPSolverRelaxedTester* target = 0
    , const char opt_grp_name[] = "QPSolverRelaxedTester" );

protected:

  /// Overridden from SetOptionsFromStreamNode
  void setOption( int option_num, const std::string& option_value );

};	// end class QPSolverRelaxedTesterSetOptions

}	// end namespace ConstrainedOptPack

#endif	// QP_SOLVER_RELAXED_TESTER_SET_OPTIONS_H
