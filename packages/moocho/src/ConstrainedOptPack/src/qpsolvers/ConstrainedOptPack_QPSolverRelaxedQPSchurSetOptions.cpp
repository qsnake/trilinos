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

#include <assert.h>
#include <math.h>

#include "ConstrainedOptPack_QPSolverRelaxedQPSchurSetOptions.hpp"
#include "OptionsFromStreamPack_StringToBool.hpp"

// Define the options
namespace {

  const int local_num_options = 22;

  enum local_EOptions {
    MAX_QP_ITER_FRAC
    ,MAX_REAL_RUNTIME
    ,INEQUALITY_PICK_POLICY
    ,BOUND_TOL
    ,INEQUALITY_TOL
    ,EQUALITY_TOL
    ,LOOSE_FEAS_TOL
    ,DUAL_INFEAS_TOL
    ,HUGE_PRIMAL_STEP
    ,HUGE_DUAL_STEP
    ,BIGM
    ,WARNING_TOL
    ,ERROR_TOL
    ,ITER_REFINE_MIN_ITER
    ,ITER_REFINE_MAX_ITER
    ,ITER_REFINE_OPT_TOL
    ,ITER_REFINE_FEAS_TOL
    ,ITER_REFINE_AT_SOLUTION
    ,PIVOT_WARNING_TOL
    ,PIVOT_SINGULAR_TOL
    ,PIVOT_WRONG_INERTIA_TOL
    ,PRINT_LEVEL
  };

  const char* local_SOptions[local_num_options]	= {
    "max_qp_iter_frac"
    ,"max_real_runtime"
    ,"inequality_pick_policy"
    ,"bounds_tol"
    ,"inequality_tol"
    ,"equality_tol"
    ,"loose_feas_tol"
    ,"dual_infeas_tol"
    ,"huge_primal_step"
    ,"huge_dual_step"
    ,"bigM"
    ,"warning_tol"
    ,"error_tol"
    ,"iter_refine_min_iter"
    ,"iter_refine_max_iter"
    ,"iter_refine_opt_tol"
    ,"iter_refine_feas_tol"
    ,"iter_refine_at_solution"
    ,"pivot_warning_tol"
    ,"pivot_singular_tol"
    ,"pivot_wrong_inertia_tol"
    ,"print_level"
  };

}

namespace ConstrainedOptPack {

QPSolverRelaxedQPSchurSetOptions::QPSolverRelaxedQPSchurSetOptions(
        QPSolverRelaxedQPSchur* target
      , const char opt_grp_name[] )
  :	OptionsFromStreamPack::SetOptionsFromStreamNode(
        opt_grp_name, local_num_options, local_SOptions )
    , OptionsFromStreamPack::SetOptionsToTargetBase<
      QPSolverRelaxedQPSchur >( target )
{}

void QPSolverRelaxedQPSchurSetOptions::setOption(
  int option_num, const std::string& option_value )
{
  using OptionsFromStreamPack::StringToBool;
  typedef QPSolverRelaxedQPSchur target_t;
  typedef QPSchurPack::ConstraintsRelaxedStd constr_t;
  switch( (local_EOptions)option_num ) {
    case MAX_QP_ITER_FRAC:
      target().max_qp_iter_frac(std::fabs(std::atof(option_value.c_str())));
      break;
    case MAX_REAL_RUNTIME:
      target().max_real_runtime(std::fabs(std::atof(option_value.c_str())));
      break;
    case INEQUALITY_PICK_POLICY:
      if(			option_value == "ADD_BOUNDS_THEN_MOST_VIOLATED_INEQUALITY" )
        target().inequality_pick_policy( constr_t::ADD_BOUNDS_THEN_MOST_VIOLATED_INEQUALITY );
      else if(	option_value == "ADD_BOUNDS_THEN_FIRST_VIOLATED_INEQUALITY" )
        target().inequality_pick_policy( constr_t::ADD_BOUNDS_THEN_FIRST_VIOLATED_INEQUALITY );
      else if(	option_value == "ADD_MOST_VIOLATED_BOUNDS_AND_INEQUALITY" )
        target().inequality_pick_policy( constr_t::ADD_MOST_VIOLATED_BOUNDS_AND_INEQUALITY );
      else
        throw std::invalid_argument( "QPSolverRelaxedQPSchurSetOptions::"
          "setOption(...) : Error, only the values of\n"
          " ADD_BOUNDS_THEN_MOST_VIOLATED_INEQUALITY\n"
          ", ADD_BOUNDS_THEN_FIRST_VIOLATED_INEQUALITY and"
          " ADD_MOST_VIOLATED_BOUNDS_AND_INEQUALITY \nare valid for the option"
          " \"QPSolverRelaxedQPSchur::inequality_pick_policy\"" );
      break;
    case BOUND_TOL:
      target().bounds_tol(std::fabs(std::atof(option_value.c_str())));
      break;
    case INEQUALITY_TOL:
      target().inequality_tol(std::fabs(std::atof(option_value.c_str())));
      break;
    case EQUALITY_TOL:
      target().equality_tol(std::fabs(std::atof(option_value.c_str())));
      break;
    case LOOSE_FEAS_TOL:
      target().loose_feas_tol(std::fabs(std::atof(option_value.c_str())));
      break;
    case DUAL_INFEAS_TOL:
      target().dual_infeas_tol(std::fabs(std::atof(option_value.c_str())));
      break;
    case HUGE_PRIMAL_STEP:
      target().huge_primal_step(std::fabs(std::atof(option_value.c_str())));
      break;
    case HUGE_DUAL_STEP:
      target().huge_dual_step(std::fabs(std::atof(option_value.c_str())));
      break;
    case BIGM:
      target().bigM(std::fabs(std::atof(option_value.c_str())));
      break;
    case WARNING_TOL:
      target().warning_tol(std::fabs(std::atof(option_value.c_str())));
      break;
    case ERROR_TOL:
      target().error_tol(std::fabs(std::atof(option_value.c_str())));
      break;

    case ITER_REFINE_MIN_ITER:
      target().iter_refine_min_iter(std::abs(std::atoi(option_value.c_str())));
      break;
    case ITER_REFINE_MAX_ITER:
      target().iter_refine_max_iter(std::abs(std::atoi(option_value.c_str())));
      break;
    case ITER_REFINE_OPT_TOL:
      target().iter_refine_opt_tol(std::fabs(std::atof(option_value.c_str())));
      break;
    case ITER_REFINE_FEAS_TOL:
      target().iter_refine_feas_tol(std::fabs(std::atof(option_value.c_str())));
      break;
    case ITER_REFINE_AT_SOLUTION:
      target().iter_refine_at_solution(StringToBool( "iter_refine_at_solution", option_value.c_str() ));
      break;
    case PIVOT_WARNING_TOL:
      target().pivot_warning_tol(std::fabs(std::atof(option_value.c_str())));
      break;
    case PIVOT_SINGULAR_TOL:
      target().pivot_singular_tol(std::fabs(std::atof(option_value.c_str())));
      break;
      case PIVOT_WRONG_INERTIA_TOL:
      target().pivot_wrong_inertia_tol(std::fabs(std::atof(option_value.c_str())));
      break;
    case PRINT_LEVEL:
      if(			option_value == "USE_INPUT_ARG" )
        target().print_level( target_t::USE_INPUT_ARG );
      else if(	option_value == "NO_OUTPUT" )
        target().print_level( target_t::NO_OUTPUT );
      else if(	option_value == "OUTPUT_BASIC_INFO" )
        target().print_level( target_t::OUTPUT_BASIC_INFO );
      else if(	option_value == "OUTPUT_ITER_SUMMARY" )
        target().print_level( target_t::OUTPUT_ITER_SUMMARY );
      else if(	option_value == "OUTPUT_ITER_STEPS" )
        target().print_level( target_t::OUTPUT_ITER_STEPS );
      else if(	option_value == "OUTPUT_ACT_SET" )
        target().print_level( target_t::OUTPUT_ACT_SET );
      else if(	option_value == "OUTPUT_ITER_QUANTITIES" )
        target().print_level( target_t::OUTPUT_ITER_QUANTITIES );
      else
        throw std::invalid_argument( "QPSolverRelaxedQPSchurSetOptions::"
          "setOption(...) : Error, only the values of USE_INPUT_ARG, NO_OUTPUT"
          ", OUTPUT_BASIC_INFO, OUTPUT_ITER_SUMMARY\n"
          ", OUTPUT_ITER_STEPS, OUTPUT_ACT_SET and"
          " OUTPUT_ITER_QUANTITIES \nare valid for the option"
          " \"QPSolverRelaxedQPSchur::print_level\"" );
      break;
    default:
      TEST_FOR_EXCEPT(true);	// Local error only?
  }
}

}	// end namespace ConstrainedOptPack
