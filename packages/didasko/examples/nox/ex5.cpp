
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
// The files in this directory solve the Chandrasekhar H-equation using NOX & LOCA. 
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
// This file solve for the parameter "c" to find the turning point,
// which is at c=1.  Turning point continuation requires a second parameter,
// but in this case we aren't interested in the 2nd parameter continuation, so 
// the 2nd parameter is called "dummy" and is set to a fixed value.
// The eigenvalues are also computed.
//
// The output file (currently defined as "Heq5.dat") will output the value of 
// the continuation parameter and the 1-norm of the solution vector x.  
// However, since we have two continuation parameters defined ("c" and "dummy"),
// the only output we are interested in is the third piece of data, which is
// the value of "c" at the turning point.  (The first two values are those 
// associated with the "dummy" parameter.)     

// include all relevant header files

// need to include all didasko (example) info
#include "Didasko_ConfigDefs.h"
#if defined(HAVE_DIDASKO_EPETRA) && defined(HAVE_DIDASKO_NOX) && defined(HAVE_NOX_EPETRA)

// include all header files that are needed for Loca continuation
#include "Epetra_LinearProblem.h"
#include "AztecOO.h"
#include "LOCA.H"
#include "LOCA_GlobalData.H"
#include "LOCA_Epetra.H"
#include "NOX_Epetra_MatrixFree.H"
#include "LOCA_Epetra_Interface_TimeDependentMatrixFree.H"
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

  double c = 0.9999;           // continuation parameter
  int N = 50;                  // number of grid points
  int maxNewtonIters = 20;     // max number of Newton iterations
  int maxSteps = 50;           // max number of continuation steps taken
  int ilocal, iglobal;         // counter variables used for loops:
                               //   ilocal = counter for local elements on this processor;
                               //   iglobal = counter to signify global position across all procs 
  int Myele;                   // holds the number of elements on the processor

// Set flag for whether the computations will be Matrix-free (true) or will use a computed
//   Jacobian (false)
  bool doMatFree = false;      
   
// Create output file to save solutions
  ofstream outFile("Heq5.dat");
  outFile.setf(ios::scientific, ios::floatfield);
  outFile.precision(10);

// Define the problem class
  HeqProblem Problem(N,&Comm,outFile);
  
// Build initial guess.  The initial guess should be a solution vector x close to the 
//   bifurcation point.

  // Create the initial guess vector
  Epetra_Vector InitialGuess(Problem.GetMap());

  // Get the number of elements on this processor  
  Myele = Problem.GetMap().NumMyElements();

  // Compute the initial guess.  For this example, it is a line from (0,1) to (1,8/3)
  for (ilocal=0; ilocal<Myele; ilocal++) {
     iglobal=Problem.GetMap().GID(ilocal);
     InitialGuess[ilocal]= 1.0 + (5.0*iglobal)/(3.0*(N-1));
  }

// Create the null vector for the Jacobian (ie, J*v=0, used to solve the system of equations
//   f(x,p)=0; J*v=0; v0*v=1.  The solution of the system is [x*,v*,p*]. )
  Teuchos::RCP<NOX::Abstract::Vector> nullVec =
    Teuchos::rcp(new NOX::Epetra::Vector(InitialGuess));

  // Initialize to all ones
  nullVec->init(1.0);  
    // NOTE:  init is a function within the NOX::Abstract:Vector class which initializes every
    // value of the vector to the value within the parentheses (must be in 'double' format) 

// Create the top level parameter list
  Teuchos::RCP<Teuchos::ParameterList> ParamList = Teuchos::rcp(new Teuchos::ParameterList);

  // Create LOCA sublist
  Teuchos::ParameterList& locaParamsList = ParamList->sublist("LOCA");

  // Create the sublist for continuation and set the stepper parameters
  Teuchos::ParameterList& stepperList = locaParamsList.sublist("Stepper");
    //stepperList.set("Continuation Method", "Arc Length");// Default
    stepperList.set("Continuation Method", "Natural");
    stepperList.set("Continuation Parameter", "dummy");  // Must set
    stepperList.set("Initial Value", 999.0);             // Must set
    stepperList.set("Max Value", 50.0e4);             // Must set
    stepperList.set("Min Value", 0.0);             // Must set
    stepperList.set("Max Steps", maxSteps);                    // Should set
    stepperList.set("Max Nonlinear Iterations", maxNewtonIters); // Should set
    stepperList.set("Bordered Solver Method", "Bordering");

  //  Teuchos::ParameterList& nestedList = 
  //    stepperList.sublist("Nested Bordered Solver");
  //  nestedList.set("Bordered Solver Method", "Householder");
  //  nestedList.set("Include UV In Preconditioner", true);
  //  //nestedList.set("Use P For Preconditioner", true);
  //  nestedList.set("Preconditioner Method", "SMW");

// Set up parameters to compute Eigenvalues
#ifdef HAVE_LOCA_ANASAZI
  // Create Anasazi Eigensolver sublist (needs --with-loca-anasazi)
  stepperList.set("Compute Eigenvalues",true);
  Teuchos::ParameterList& aList = stepperList.sublist("Eigensolver");
  aList.set("Method", "Anasazi");
  aList.set("Block Size", 1);        // Size of blocks
  aList.set("Num Blocks", 20);       // Size of Arnoldi factorization
  aList.set("Num Eigenvalues", 5);   // Number of eigenvalues
  //  aList.set("Sorting Order", "SR");
  aList.set("Convergence Tolerance", 2.0e-7);          // Tolerance
  aList.set("Step Size", 1);         // How often to check convergence
  aList.set("Maximum Restarts",2);   // Maximum number of restarts
  aList.set("Verbosity",  
	    Anasazi::Errors + 
	    Anasazi::Warnings +
	    Anasazi::FinalSummary);        // Verbosity
#else
    stepperList.set("Compute Eigenvalues",false);
#endif
  
  // Create bifurcation sublist.  Note that for turning point continuation, the "type"
  //   is set to "Turning Point".  If not doing TP, type should be "None".
  Teuchos::ParameterList& bifurcationList = locaParamsList.sublist("Bifurcation");
  bifurcationList.set("Type", "Turning Point");
  bifurcationList.set("Bifurcation Parameter", "c");
  //  bifurcationList.set("Formulation", "Minimally Augmented");
  bifurcationList.set("Symmetric Jacobian", false); 
  bifurcationList.set("Update Null Vectors Every Continuation Step", true);
  bifurcationList.set("Update Null Vectors Every Nonlinear Iteration", false);
  bifurcationList.set("Transpose Solver Method","Explicit Transpose");
  //  bifurcationList.set("Transpose Solver Method","Transpose Preconditioner");
  //  bifurcationList.set("Transpose Solver Method","Left Preconditioning");
  bifurcationList.set("Initial Null Vector Computation", "Solve df/dp");
  //  bifurcationList.set("Initial A Vector", nullVec);      // minimally augmented
  //  bifurcationList.set("Initial B Vector", nullVec);      //minimally augmented
  
  //  bifurcationList.set("Bordered Solver Method", "Householder");
  //  bifurcationList.set("Include UV In Preconditioner", true);
  //  //bifurcationList.set("Use P For Preconditioner", true);
  //  bifurcationList.set("Preconditioner Method", "SMW");

  bifurcationList.set("Formulation", "Moore-Spence");
  bifurcationList.set("Solver Method", "Phipps Bordering"); // better for nearly singular matrices
  //  bifurcationList.set("Solver Method", "Salinger Bordering");   
  bifurcationList.set("Initial Null Vector", nullVec);
  bifurcationList.set("Length Normalization Vector", nullVec);

    // Create the sublist for the predictor
    Teuchos::ParameterList& predictorList = locaParamsList.sublist("Predictor");
    predictorList.set("Method", "Secant");         // Default
    // predictorList.set("Method", "Constant");     // Other options
    // predictorList.set("Method", "Tangent");      // Other options

    // Create step size sublist
    Teuchos::ParameterList& stepSizeList = locaParamsList.sublist("Step Size");
    stepSizeList.set("Method", "Adaptive");             // Default
    stepSizeList.set("Initial Step Size", 0.1);   // Should set
    stepSizeList.set("Min Step Size", 1.0e-6);    // Should set
    stepSizeList.set("Max Step Size", 1.0);      // Should set
    stepSizeList.set("Aggressiveness", 0.1);

// Set up NOX info
  Teuchos::ParameterList& nlParams = ParamList->sublist("NOX");

// Set the nonlinear solver method
  nlParams.set("Nonlinear Solver", "Line Search Based");

// Set the printing parameters in the "Printing" sublist.  This list determines how much
//   of the NOX information is output
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
			NOX::Utils::Warning +
         NOX::Utils::StepperIteration +
         NOX::Utils::StepperDetails +
         NOX::Utils::StepperParameters);

  // NOX parameters - Sublist for line search 
  Teuchos::ParameterList& searchParams = nlParams.sublist("Line Search");
  searchParams.set("Method", "Backtrack");
  //  searchParams.set("Method", "Full Step");

  // Sublist for direction
  Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
  dirParams.set("Method", "Newton");

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
  //  lsParams.set("Scaling", "None");
  //  lsParams.set("Scaling", "Row Sum");  
  lsParams.set("Compute Scaling Manually", false);
  //  lsParams.set("Preconditioner", "Ifpack");
  //  lsParams.set("Ifpack Preconditioner", "ILU");
  //  lsParams.set("Preconditioner", "New Ifpack");
  //  Teuchos::ParameterList& ifpackParams = lsParams.sublist("Ifpack");
  //  ifpackParams.set("fact: level-of-fill", 1);

// Set up the continuation parameter vector
  LOCA::ParameterVector p;
  p.addParameter("c",c);  
  p.addParameter("dummy",999.0);  

// Create the problem interface
  Teuchos::RCP<SimpleProblemInterface> interface = 
    Teuchos::rcp(new SimpleProblemInterface(&Problem,c) );

  Teuchos::RCP<LOCA::Epetra::Interface::Required> iReq = interface;

// Create the operator to hold either the Jacobian matrix or the Matrix-free operator
  Teuchos::RCP<Epetra_Operator> A;
  //  Teuchos::RCP<Epetra_RowMatrix> A;
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
 
  // Create scaling object
  Teuchos::RCP<NOX::Epetra::Scaling> scaling = Teuchos::null;
  //   scaling = Teuchos::rcp(new NOX::Epetra::Scaling);
  //   Teuchos::RCP<Epetra_Vector> scalingVector = 
  //     Teuchos::rcp(new Epetra_Vector(soln.Map()));
  //   //scaling->addRowSumScaling(NOX::Epetra::Scaling::Left, scalingVector);
  //   scaling->addColSumScaling(NOX::Epetra::Scaling::Right, scalingVector);

  // Create transpose scaling object
  Teuchos::RCP<NOX::Epetra::Scaling> trans_scaling = Teuchos::null;
  //   trans_scaling = Teuchos::rcp(new NOX::Epetra::Scaling);
  //   Teuchos::RCP<Epetra_Vector> transScalingVector = 
  //     Teuchos::rcp(new Epetra_Vector(soln.Map()));
  //   trans_scaling->addRowSumScaling(NOX::Epetra::Scaling::Right, 
  // 				  transScalingVector);
  //   trans_scaling->addColSumScaling(NOX::Epetra::Scaling::Left, 
  // 				  transScalingVector);
    //bifurcationList.set("Transpose Scaling", trans_scaling);

// Build the linear system solver
  Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linSys = 
    Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams, lsParams,
						      iReq,
						      iJac, A, 
                        noxInitGuess, scaling));            // use if scaling
//                        noxInitGuess));                     // use if no scaling

// Create the Loca (continuation) vector
  NOX::Epetra::Vector locaSoln(noxInitGuess);
  
  // Create Epetra Factory
  Teuchos::RCP<LOCA::Abstract::Factory> epetraFactory = Teuchos::rcp(new LOCA::Epetra::Factory);

  // Create global data object
  Teuchos::RCP<LOCA::GlobalData> globalData = LOCA::createGlobalData(ParamList, epetraFactory);
 
  // Create the Group - must be LOCA group
  Teuchos::RCP<LOCA::Epetra::Group> grpPtr = 
    Teuchos::rcp(new LOCA::Epetra::Group(globalData, printParams, 
					iReq, locaSoln, 
					linSys, p)); 

  // Calculate the first F(x0) as a starting point.  This is only needed for
  // certain status tests, to ensure that an initial residual (|r0|) is calculated
  grpPtr->computeF();

// Set up the status tests to check for convergence
  // Determines the error tolerance for the Newton solves 
  Teuchos::RCP<NOX::StatusTest::NormF> testNormF = 
    Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-4));
  // Sets the max number of nonlinear (Newton) iterations that will be taken.  If this is not
  //   already set, it will default to the '20' given 
  Teuchos::RCP<NOX::StatusTest::MaxIters> testMaxIters = 
    Teuchos::rcp(new NOX::StatusTest::MaxIters(stepperList.get("Max Nonlinear Iterations", 20)));
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


// Create the stepper
  LOCA::Stepper stepper(globalData, grpPtr, combo, ParamList);
  LOCA::Abstract::Iterator::IteratorStatus status = stepper.run();
  
  // Check if the stepper completed
  if  (status == LOCA::Abstract::Iterator::Finished)
    globalData->locaUtils->out() << "\nAll tests passed!" << endl;
  else 
    if (globalData->locaUtils->isPrintType(NOX::Utils::Error))
      globalData->locaUtils->out() << "\nStepper failed to converge!"  << endl;

// Output the stepper parameter list info
  if (globalData->locaUtils->isPrintType(NOX::Utils::StepperParameters)) {
    globalData->locaUtils->out() << endl << "Final Parameters" << endl
    << "*******************" << endl;
    stepper.getList()->print(globalData->locaUtils->out());
    globalData->locaUtils->out() << endl;
  }

// Make sure all processors are done and close the output file
Comm.Barrier();
outFile.close();

// Deallocate memory
LOCA::destroyGlobalData(globalData);

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
