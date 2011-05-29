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
//

#ifndef EXAMPLE_NLP_FIRST_ORDER_DIRECT_RUN_H
#define EXAMPLE_NLP_FIRST_ORDER_DIRECT_RUN_H

#include <iosfwd>

#include "NLPInterfacePack_Types.hpp"
#include "MoochoPack_MoochoSolver.hpp"

namespace NLPInterfacePack {

/** \defgroup ExampleNLPDirectRun_grp Helper function for ExampleNLPDirect */
//@{

/** \brief Function accepts a VectorSpace object and then uses it to define
 * an example NLP and run <tt>MoochoPack::MoochoSolver</tt> on it.
 *
 * @param  vec_space   [in] The vector space object used to create all of the
 *                     needed vector spaces and vectors.  This vector space and
 *                     the vectors it creates will get a though testing.
 * @param  xo          [in] The initial starting point for unknown variables (before
 *                     they are forced in bounds).
 * @param  has_bounds  [in] If true, then the NLP will have bounds on the variables.
 * @param  dep_bounded [in] (valid only if has_bounds == true) If true, then
 *                     the dependent variables will be bounded, if false the
 *                     independent variables will be bounded.
 * @param  console_out [in/out] If != NULL then *console_out gets the output.
 * @param  error_out   [in/out] If != NULL then *eout gets minimal summary output.
 * @param  throw_solve_exception
 *                     [in] If true then solver will not throw exception (but other code may).
 * @param  algo_out    [in/out] If != NULL then it gets algo outptut, otherwise goes to 'MoochoAlgo.out'
 * @param  summary_out [in/out] If != NULL then it gets summary outptut, otherwise goes to 'MoochoSummary.out'
 * @param  journal_out [in/out] If != NULL then it gets journal outptut, otherwise goes to 'MoochoJournal.out'
 *
 * @returns Returns the return value from <tt>MoochoPack::rsqp_mama_jama_solve()</tt>
 * (see this function for most of the documentation).
 */
MoochoPack::MoochoSolver::ESolutionStatus
ExampleNLPDirectRun(
  const VectorSpace&   vec_space
  ,value_type          xo
  ,bool                has_bounds
  ,bool                dep_bounded
  ,std::ostream*       console_out
  ,std::ostream*       error_out
  ,bool                throw_solve_exception = false
  ,std::ostream*       algo_out              = NULL
  ,std::ostream*       summary_out           = NULL
  ,std::ostream*       journal_out           = NULL
  );

//@}

} // end namespace NLPInterfacePack

#endif // EXAMPLE_NLP_FIRST_ORDER_DIRECT_RUN_H


