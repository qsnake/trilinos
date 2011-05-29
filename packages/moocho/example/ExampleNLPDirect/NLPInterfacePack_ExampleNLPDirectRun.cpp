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

#include <assert.h>

#include <fstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <typeinfo>

#include "NLPInterfacePack_ExampleNLPDirectRun.hpp"
#include "NLPInterfacePack_ExampleNLPDirect.hpp"
#include "MoochoPack_NLPAlgoConfigMamaJama.hpp"
#include "IterationPack_AlgorithmTracker.hpp"
#include "AbstractLinAlgPack_VectorSpace.hpp"
#include "AbstractLinAlgPack_BasisSystem.hpp"
#include "OptionsFromStreamPack_OptionsFromStream.hpp"

MoochoPack::MoochoSolver::ESolutionStatus
NLPInterfacePack::ExampleNLPDirectRun(
  const VectorSpace&   vec_space
  ,value_type          xo
  ,bool                has_bounds
  ,bool                dep_bounded
  ,std::ostream*       console_out
  ,std::ostream*       error_out
  ,bool                throw_solve_exception
  ,std::ostream*       algo_out
  ,std::ostream*       summary_out
  ,std::ostream*       journal_out
  )
{
  using std::endl;
  using std::setw;
  namespace rcp = MemMngPack;
  using Teuchos::RCP;
  namespace ofsp = OptionsFromStreamPack;
  using ofsp::OptionsFromStream;
  namespace rsqp = MoochoPack;
  using rsqp::MoochoSolver;
  using rsqp::NLPAlgoConfigMamaJama;

  MoochoSolver::ESolutionStatus
    solve_return = MoochoSolver::SOLVE_RETURN_EXCEPTION;

  int err = 0;
  
  int w = 15;
  int prec = 8;

  if(console_out)
    *console_out
      << std::setprecision(prec)
      << std::scientific
      << "***************************************************\n"
      << "*** Running Tests on ExampleNLPDirect ***\n"
      << "***************************************************\n"
      << "\nUsing a vector space of type \'" << typeName(vec_space) << "\'"
      << "\nwith a dimension of vec_space.dim() = " << vec_space.dim()
      << std::endl;

  // Create the nlp
  ExampleNLPDirect
    nlp(VectorSpace::space_ptr_t(&vec_space,false),xo,has_bounds,dep_bounded);

  // Create the solver object and set it up
  MoochoSolver solver;
  solver.set_nlp(Teuchos::rcp(&nlp,false));                  // Set the NLP!
  solver.set_error_handling(                             // set up outputting
    throw_solve_exception
    ,Teuchos::rcp(error_out,false)
    );
  solver.set_console_out(Teuchos::rcp(console_out,false));
  solver.set_summary_out(Teuchos::rcp(summary_out,false));
  solver.set_journal_out(Teuchos::rcp(journal_out,false));
  solver.set_algo_out(   Teuchos::rcp(algo_out,false)   );

  // Run MOOCHO using the MamaJama configuration
  solve_return = solver.solve_nlp();

  return solve_return;
}
