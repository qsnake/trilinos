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

#ifndef QP_SOLVER_RELAXED_QP_SCHUR_SET_OPTIONS_H
#define QP_SOLVER_RELAXED_QP_SCHUR_SET_OPTIONS_H

#include "ConstrainedOptPack_QPSolverRelaxedQPSchur.hpp"
#include "OptionsFromStreamPack_SetOptionsFromStreamNode.hpp"
#include "OptionsFromStreamPack_SetOptionsToTargetBase.hpp"

namespace ConstrainedOptPack {

/** \brief Set options for QPSolverRelaxedQPSchur from an
  * OptionsFromStream object.
  *
  * The default options group name is QPSolverRelaxedQPSchur.
  *
  * The options group is:
  *
  \begin{verbatim}
  options_group QPSolverRelaxedQPSchur {
  *    max_qp_iter_frac  = 10.0;   *** (+dbl) max_qp_itr = max_qp_itr_frac * (# variables)
  *    max_real_runtime  = 1e+20;  *** (+dbl) maximum runtime in minutes
  *    inequality_pick_policy = ADD_BOUNDS_THEN_MOST_VIOLATED_INEQUALITY;
  *    inequality_pick_policy = ADD_BOUNDS_THEN_FIRST_VIOLATED_INEQUALITY; *** not supported yet!
  *    inequality_pick_policy = ADD_MOST_VIOLATED_BOUNDS_AND_INEQUALITY;
  *    bounds_tol        = 1e-10;  *** (+dbl) feasibility tolerance for bound constriants
  *    inequality_tol    = 1e-10;  *** (+dbl) feasibility tolerance for general inequality constriants
  *    equality_tol      = 1e-10;  *** (+dbl) feasibility tolerance for general equality constriants
  *    loose_feas_tol    = 1e-9;   *** (+dbl) (Expert use only)
  *    dual_infeas_tol   = 1e-12;  *** (+dbl) allowable dual infeasiblity before error
  *    huge_primal_step  = 1e+20;  *** (+dbl) value of a near infinite primal step
  *    huge_dual_step    = 1e+20;  *** (+dbl) value of a near infinite dual step
  *    bigM              = 1e+10;  *** (+dbl) value or relaxation penalty in objective
  *    warning_tol   = 1e-10;  *** Testing warning tolerance
  *    error_tol     = 1e-5;   *** Testing error tolerance
  *    iter_refine_min_iter = 1;  *** (+int) Minimum number of iterative refinement iterations
  *    iter_refine_max_iter = 3;  *** (+int) Maximum number of iterative refinement iterations
  *    iter_refine_opt_tol  = 1e-12; *** (+dbl) Optimality tolarance for iterative refinement
  *    iter_refine_feas_tol = 1e-12; *** (+dbl) Feasibility tolerance for iterative refinement
  *    iter_refine_at_solution = true; *** (+dbl) If true then iterative refinement will always be used
  *    pivot_warning_tol        = 1e-6;  *** (+dbl) Relative warning tolerance for a pivot element in the schur complement
  *    pivot_singular_tol       = 1e-8;  *** (+dbl) Relative singularity tolerance for a pivot element in the schur complement
  *    pivot_wrong_iniertia_tol = 1e-10; *** (+dbl) Relative tolerance for a pivot element in the schur complement for wrong inertia
  *    print_level = USE_INPUT_ARG;  *** Use the input argument to solve_qp(...)
  *    print_level = NO_OUTPUT;
  *    print_level = OUTPUT_BASIC_INFO;
  *    print_level = OUTPUT_ITER_SUMMARY;
  *    print_level = OUTPUT_ITER_STEPS;
  *    print_level = OUTPUT_ACT_SET;
  *    print_level = OUTPUT_ITER_QUANTITIES;
  }
  \end{verbatim}
  */
class QPSolverRelaxedQPSchurSetOptions
  : public OptionsFromStreamPack::SetOptionsFromStreamNode 
    , public OptionsFromStreamPack::SetOptionsToTargetBase<
      QPSolverRelaxedQPSchur >
{
public:

  /** \brief . */
  QPSolverRelaxedQPSchurSetOptions(
      QPSolverRelaxedQPSchur* target = 0
    , const char opt_grp_name[] = "QPSolverRelaxedQPSchur" );

protected:

  /// Overridden from SetOptionsFromStreamNode
  void setOption( int option_num, const std::string& option_value );

};	// end class QPSolverRelaxedQPSchurSetOptions

}	// end namespace ConstrainedOptPack

#endif	// QP_SOLVER_RELAXED_QP_SCHUR_SET_OPTIONS_H
