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

#include <iostream>

#include "NLPInterfacePack_NLPWBCounterExample.hpp"
#include "MoochoPack_MoochoSolver.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

int main( int argc, char* argv[] )
{
  using MoochoPack::MoochoSolver;
  using NLPInterfacePack::NLPWBCounterExample;
  using Teuchos::CommandLineProcessor;
  bool success = true;

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  try {

    MoochoSolver  solver;
  
    //
    // Get options from the command line
    //
    
    double   xinit[3]          = { 0.0, 0.0, 0.0 };
    double   a                 = 0.0;
    double   b                 = 1.0;
    bool     nlp_selects_basis = true;
    bool     linear_obj        = true;

    CommandLineProcessor  clp(false); // don't throw exceptions
    solver.setup_commandline_processor(&clp);
    clp.setOption( "x1-init",  &xinit[0], "Initail guess for x(1)" );
    clp.setOption( "x2-init",  &xinit[1], "Initail guess for x(2)" );
    clp.setOption( "x3-init",  &xinit[2], "Initail guess for x(3)" );
    clp.setOption( "a",  &a, "Constant for c(1)" );
    clp.setOption( "b",  &b, "Constant for c(2)" );
    clp.setOption(
      "nlp-selects-basis", "no-nlp-selects-basis", &nlp_selects_basis
      ,"Determine if NLP will select basis" );
    clp.setOption(
      "linear-obj", "nonlinear-obj", &linear_obj
      ,"Determine if objective is linear" );
    CommandLineProcessor::EParseCommandLineReturn
      parse_return = clp.parse(argc,argv,&std::cerr);
    if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL )
      return parse_return;
    
    //
    // Create the NLP and solve it
    //

    // Create the NLP
    NLPWBCounterExample
      nlp(xinit,a,b,nlp_selects_basis,linear_obj);
    // Set the NLP
    solver.set_nlp( Teuchos::rcp(&nlp,false) );
    // Solve the NLP
    const MoochoSolver::ESolutionStatus
      solution_status = solver.solve_nlp();
    
    //
    // Return the solution status (0 if sucessfull)
    //

    return solution_status;

  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true,std::cout,success)

  return MoochoSolver::SOLVE_RETURN_EXCEPTION;

}
