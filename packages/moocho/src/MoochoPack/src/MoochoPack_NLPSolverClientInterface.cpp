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

#include "MoochoPack_NLPSolverClientInterface.hpp"
#include "IterationPack_AlgorithmTracker.hpp"

MoochoPack::NLPSolverClientInterface::NLPSolverClientInterface(
  int                      max_iter
  ,double                  max_run_time
  ,value_type              opt_tol
  ,value_type              feas_tol
  ,value_type              comp_tol
  ,value_type              step_tol
  ,EJournalOutputLevel     journal_output_level
  ,EJournalOutputLevel     null_space_journal_output_level
  ,int                     journal_print_digits
  ,bool                    check_results
  ,bool                    calc_conditioning
  ,bool                    calc_matrix_norms
  ,bool                    calc_matrix_info_null_space_only
  )
  :max_iter_(max_iter)
  ,max_run_time_(max_run_time)
  ,opt_tol_(opt_tol)
  ,feas_tol_(feas_tol)
  ,comp_tol_(comp_tol)
  ,step_tol_(step_tol)
  ,journal_output_level_(journal_output_level)
  ,null_space_journal_output_level_(null_space_journal_output_level)
  ,journal_print_digits_(journal_print_digits)
  ,check_results_(check_results)
  ,calc_conditioning_(calc_conditioning)
  ,calc_matrix_norms_(calc_matrix_norms)
  ,calc_matrix_info_null_space_only_(calc_matrix_info_null_space_only)
{}
