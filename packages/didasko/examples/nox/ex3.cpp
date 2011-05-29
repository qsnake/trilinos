
// @HEADER
// ***********************************************************************
// 
//                      Didasko Tutorial Package
//                 Copyright (2005) Sandia Corporation
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
//
// Questions about Didasko? Contact Marzio Sala (marzio.sala _AT_ gmail.com)
// 
// ***********************************************************************
// @HEADER

// -----------------
// The files in this directory solve the Chandrasekhar H-equation using NOX. 
//
// The Chandrasekhar H-equation is defined as:
//   F(H)(u) = H(u) - [ 1 - c/2 * Int_0^1 { u*H(v)dv / (u+v) } ]^{-1}
// The discretized H-equation is:
//   F(x)_i = x_i - [ 1 - c/(2N) * Sum_{j=1}^N { u_i*x_j / (u_i+u_j) } ]^{-1}
// where u_i = (i - 0.5)/N for 1 <= i <= N.
// 
// The H-equation has two solutions for c in (0,1).  It has one solution for c=1,
//   which is also the turning point for the solution graph.
// 
//
// This file will solve the nonlinear system for a specific value of the parameter
// "c" to find the corresponding solution vector.  It will output the 2-norm of the 
// solution vector as well as the solution vector itself.
//

// include all relevant header files

// need to include all didasko (example) info
#include "Didasko_ConfigDefs.h"
#if defined(HAVE_DIDASKO_EPETRA) && defined(HAVE_DIDASKO_NOX) && defined(HAVE_NOX_EPETRA)

// include all header files that are needed for Loca continuation
#include "Epetra_LinearProblem.h"
#include "NOX_Epetra_MatrixFree.H"
// header file with the specific problem info
#include "ProblemInterface.H"
// for output file
#include <fstream>

// Required for reading and writing parameter lists from xml format
// Configure Trilinos with --enable-teuchos-extended
#ifdef HAVE_TEUCHOS_EXTENDED
#include "Teuchos_XMLParameterListHelpers.hpp"
#endif

// =========== //
// main driver //
// =========== //

int main( int argc, char **argv )
{

// check for parallel computation
#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif


// define main parameters

  double c = 0.5;             // continuation parameter
  int N = 50;                  // number of grid points
  int maxNewtonIters = 20;     // max number of Newton iterations
  double n2;                   // temporary variable to hold the norm of the soln vector
  int NumProc, MyPID;          // NumProc = total no. of processors used, 
                               // MyPID = ID of the current processor
  int NumEle;                  // No. of elements on the current processor
  int i, j;                    // temporary variables used for loops

// Set flag for whether the computations will be Matrix-free (true) or will use a computed
//   Jacobian (false)
  bool doMatFree = false;      
      
// Create output file to save solutions
  ofstream outFile("Heq3.dat",ios::out);
  outFile.close();
  outFile.open("Heq3.dat",ios::app);
  outFile.setf(ios::scientific, ios::floatfield);
  outFile.precision(10);

// Define the problem class
  HeqProblem Problem(N,&Comm,outFile);
  
// Create the initial guess vector and set it to all ones
  Epetra_Vector InitialGuess(Problem.GetMap());
  InitialGuess.PutScalar(1.0);     

// Create the top level parameter list
  Teuchos::RCP<Teuchos::ParameterList> ParamList = Teuchos::rcp(new Teuchos::ParameterList);

  // Set up NOX info
  Teuchos::ParameterList& nlParams = ParamList->sublist("NOX");

  // Set the nonlinear solver method
  nlParams.set("Nonlinear Solver", "Line Search Based");

  // Set the printing parameters in the "Printing" sublist
  Teuchos::ParameterList& printParams = nlParams.sublist("Printing");
  printParams.set("MyPID", Comm.MyPID()); 
  printParams.set("Output Precision", 5);
  printParams.set("Output Processor", 0);
  printParams.set("Output Information", 
			NOX::Utils::OuterIteration + 
			NOX::Utils::OuterIterationStatusTest + 
			NOX::Utils::InnerIteration +
			NOX::Utils::LinearSolverDetails +
			NOX::Utils::Parameters + 
			NOX::Utils::Details + 
			NOX::Utils::Warning);

  // NOX parameters - Sublist for line search 
  Teuchos::ParameterList& searchParams = nlParams.sublist("Line Search");
  searchParams.set("Method", "Full Step");
//  searchParams.set("Method", "Backtrack");
//  searchParams.set("Method", "NonlinearCG");

  // Sublist for direction
  Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
  dirParams.set("Method", "Newton");
//  dirParams.set("Method", "NonlinearCG");

  Teuchos::ParameterList& newtonParams = dirParams.sublist("Newton");
  newtonParams.set("Forcing Term Method", "Constant");

  // Sublist for linear solver for the Newton method
  Teuchos::ParameterList& lsParams = newtonParams.sublist("Linear Solver");
  lsParams.set("Aztec Solver", "GMRES");  
  lsParams.set("Max Iterations", 800);  
  lsParams.set("Tolerance", 1e-8);
  lsParams.set("Output Frequency", 1);    
  lsParams.set("Preconditioner", "None");
//  lsParams.set("Preconditioner", "AztecOO");
//  lsParams.set("Aztec Preconditioner", "ilu"); 
//  lsParams.set("Preconditioner", "Ifpack");
//  lsParams.set("Ifpack Preconditioner", "ILU");
//  lsParams.set("Preconditioner", "New Ifpack");
//  Teuchos::ParameterList& ifpackParams = lsParams.sublist("Ifpack");
//    ifpackParams.set("fact: level-of-fill", 1);


// Set up the problem interface
  Teuchos::RCP<SimpleProblemInterface> interface = 
    Teuchos::rcp(new SimpleProblemInterface(&Problem,c) );

  Teuchos::RCP<LOCA::Epetra::Interface::Required> iReq = interface;

// Create the operator to hold either the Jacobian matrix or the Matrix-free operator
  Teuchos::RCP<Epetra_Operator> A;
  Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac;

  // Need a NOX::Epetra::Vector for constructor
  // This becomes the initial guess vector that is used for the nonlinear solves
  NOX::Epetra::Vector noxInitGuess(InitialGuess, NOX::DeepCopy);   

  if (doMatFree) {
    // Matrix Free application (Epetra Operator):
    Teuchos::RCP<NOX::Epetra::MatrixFree> MF = 
      Teuchos::rcp(new NOX::Epetra::MatrixFree(printParams, interface, noxInitGuess)); 
    A = MF;
    iJac = MF;
  }
  else  {  // Computed Jacobian application
    A = Teuchos::rcp( Problem.GetMatrix(), false );
    iJac = interface;
  }
 
// Build the linear system solver
  Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linSys = 
    Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams, lsParams,
						      iReq,
						      iJac, A, 
						      noxInitGuess));

// Create the Group - must be NOX group
  Teuchos::RCP<NOX::Epetra::Group> grpPtr = 
    Teuchos::rcp(new NOX::Epetra::Group(printParams, 
					iReq, noxInitGuess, 
					linSys)); 

  // Calculate the first F(x0) as a starting point.  This is only needed for
  // certain status tests, to ensure that an initial residual (|r0|) is calculated
  grpPtr->computeF();

// Set up the status tests to check for convergence
  // Determines the error tolerance for the Newton solves 
  Teuchos::RCP<NOX::StatusTest::NormF> testNormF = 
    Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-4));
  // Checks for the max number of nonlinear (Newton) iterations to be taken 
  Teuchos::RCP<NOX::StatusTest::MaxIters> testMaxIters = 
    Teuchos::rcp(new NOX::StatusTest::MaxIters(maxNewtonIters));
// This combination of tests will be used by NOX to determine whether the step converged
  Teuchos::RCP<NOX::StatusTest::Combo> combo = 
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR, 
					    testNormF, testMaxIters));

// This is sample code to write and read parameters to/from a file.  Currently not activated! 
// To use, change the 'XXXHAVE_TEUCHOS_EXTENDED' TO 'HAVE_TEUCHOS_EXTENDED'
#ifdef XXXHAVE_TEUCHOS_EXTENDED
  // Write the parameter list to a file
  cout << "Writing parameter list to \"input.xml\"" << endl;
  Teuchos::writeParameterListToXmlFile(*ParamList, "input.xml");

  // Read in the parameter list from a file
  cout << "Reading parameter list from \"input.xml\"" << endl;
  Teuchos::RCP<Teuchos::ParameterList> paramList2 = 
    Teuchos::rcp(new Teuchos::ParameterList);
  Teuchos::updateParametersFromXmlFile("input.xml", paramList2.get());
  ParamList = paramList2;
#endif

// Create the solver
  Teuchos::RCP<NOX::Solver::Generic> solver = 
    NOX::Solver::buildSolver(grpPtr, combo, ParamList);

// Solve the nonlinear system
  NOX::StatusTest::StatusType status = solver->solve();

// Output whether the nonlinear solver converged
  if( NOX::StatusTest::Converged  == status )
    cout << "\n" << "-- NOX solver converged --" << "\n";
  else
    cout << "\n" << "-- NOX solver did not converge --" << "\n";

// Output the NOX parameter info
  if( Comm.MyPID() == 0 ) {
    cout << "\n" << "-- Parameter List From Solver --" << "\n";
    solver->getList().print(cout);
  }

// Get the Epetra_Vector with the final solution from the solver
  const NOX::Epetra::Group & finalGroup = 
    dynamic_cast<const NOX::Epetra::Group&>(solver->getSolutionGroup());
  const Epetra_Vector & finalSolution = 
      (dynamic_cast<const NOX::Epetra::Vector&>(finalGroup.getX())).getEpetraVector();

// Compute the 2-norm of the solution vector
  finalSolution.Norm2(&n2);                   
  
  if (Comm.MyPID()==0) {
    outFile << "\n The two-norm of the solution vector is: " << setprecision(10) << n2;
    outFile << "\n \n";
    outFile << "Computed solution : " << endl;
  }

// Write out the solution in an acceptable format

// Get the number of processors and the ID of the current processor
  NumProc = Comm.NumProc();
  MyPID = Comm.MyPID();
  NumEle = finalSolution.Map().NumMyElements();

// Loop through all processors in turn to print solution in sequential order
  for (i=0; i<NumProc; i++)  {
    if (i == MyPID) {
      for (j=0; j<NumEle; j++) 
        outFile << setprecision(10) << finalSolution[j] << endl;
    }
     
    else Comm.Barrier();
  }

// Make sure all processors are done and close the output file
Comm.Barrier();
outFile.close();


#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  return(EXIT_SUCCESS);
}  // DONE!!

// If initial includes did not have Didasko, error messages result
#else

#include <stdlib.h>
#include <stdio.h>
#ifdef HAVE_MPI
#include "mpi.h"
#endif

int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif
  puts("Please configure Didasko with:");
  puts("--enable-epetra");
  puts("--enable-nox-epetra");
  puts("--enable-ifpack");
  puts("--enable-aztecoo");
  puts("--enable-nox");

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  return(EXIT_SUCCESS);
}
#endif
