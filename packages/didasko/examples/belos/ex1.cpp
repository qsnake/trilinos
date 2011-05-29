
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
// Questions about Didasko? Contact Marzio Sala (marzio.sala _AT_ // gmail.com)
// 
// ***********************************************************************
// @HEADER

// Compute the solutions to an SPD linear system using Belos::BlockCG

#include "Didasko_ConfigDefs.h"
#if defined(HAVE_DIDASKO_EPETRA) && defined(HAVE_DIDASKO_BELOS) && defined(HAVE_DIDASKO_TEUCHOS) && defined(HAVE_DIDASKO_TRIUTILS) 

#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "Epetra_MultiVector.h"
#include "Trilinos_Util_CrsMatrixGallery.h"
#include "BelosLinearProblem.hpp"
#include "BelosEpetraAdapter.hpp"
#include "BelosBlockCGSolMgr.hpp"

int main(int argc, char *argv[]) {
    
#ifdef HAVE_MPI
  // Initialize MPI and setup an Epetra communicator
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  // If we aren't using MPI, then setup a serial communicator.
  Epetra_SerialComm Comm;
#endif

  // Some typedefs for oft-used data types
  typedef Epetra_MultiVector MV;
  typedef Epetra_Operator OP;
  typedef Belos::MultiVecTraits<double, Epetra_MultiVector> MVT;

  bool ierr;

  // Get our processor ID (0 if running serial)
  int MyPID = Comm.MyPID();

  // Verbose flag: only processor 0 should print
  bool verbose = (MyPID==0);

  // Initialize a Gallery object, from which we will select a matrix
  // Select a 2-D laplacian, of order 100
  Trilinos_Util::CrsMatrixGallery Gallery("laplace_2d", Comm);
  Gallery.Set("problem_size", 100);

  // Say hello and print some information about the gallery matrix
  if (verbose) {
    cout << "Belos Example: Block CG" << endl;
    cout << "Problem info:" << endl;
    cout << Gallery;
    cout << endl;
  }

  // Setup some more problem/solver parameters:

  // Block size
  int blocksize = 4;

  // Get a pointer to the system matrix, inside of a Teuchos::RCP
  // The Teuchos::RCP is a reference counting pointer that handles
  // garbage collection for us, so that we can perform memory allocation without
  // having to worry about freeing memory manually.
  Teuchos::RCP<Epetra_CrsMatrix> A = Teuchos::rcp( Gallery.GetMatrix(), false );

  // Create an Belos MultiVector, based on Epetra MultiVector
  const Epetra_Map * Map = &(A->RowMap());
  Teuchos::RCP<Epetra_MultiVector> B = 
    Teuchos::rcp( new Epetra_MultiVector(*Map,blocksize) );
  Teuchos::RCP<Epetra_MultiVector> X = 
    Teuchos::rcp( new Epetra_MultiVector(*Map,blocksize) );

  // Initialize the solution with zero and right-hand side with random entries
  X->PutScalar( 0.0 );
  B->Random();

  // Setup the linear problem, with the matrix A and the vectors X and B
  Teuchos::RCP< Belos::LinearProblem<double,MV,OP> > myProblem = 
    Teuchos::rcp( new Belos::LinearProblem<double,MV,OP>(A, X, B) );

  // The 2-D laplacian is symmetric. Specify this in the linear problem.
  myProblem->setHermitian();

  // Signal that we are done setting up the linear problem
  ierr = myProblem->setProblem();

  // Check the return from setProblem(). If this is true, there was an
  // error. This probably means we did not specify enough information for
  // the eigenproblem.
  assert(ierr == true);

  // Specify the verbosity level. Options include:
  // Belos::Errors 
  //   This option is always set
  // Belos::Warnings 
  //   Warnings (less severe than errors)
  // Belos::IterationDetails 
  //   Details at each iteration, such as the current eigenvalues
  // Belos::OrthoDetails 
  //   Details about orthogonality
  // Belos::TimingDetails
  //   A summary of the timing info for the solve() routine
  // Belos::FinalSummary 
  //   A final summary 
  // Belos::Debug 
  //   Debugging information
  int verbosity = Belos::Warnings + Belos::Errors + Belos::FinalSummary + Belos::TimingDetails;

  // Create the parameter list for the eigensolver
  Teuchos::RCP<Teuchos::ParameterList> myPL = Teuchos::rcp( new Teuchos::ParameterList() );
  myPL->set( "Verbosity", verbosity );
  myPL->set( "Block Size", blocksize );
  myPL->set( "Maximum Iterations", 100 );
  myPL->set( "Convergence Tolerance", 1.0e-8 );

  // Create the Block CG solver
  // This takes as inputs the linear problem and the solver parameters
  Belos::BlockCGSolMgr<double,MV,OP> mySolver(myProblem, myPL);

  // Solve the linear problem, and save the return code
  Belos::ReturnType solverRet = mySolver.solve();

  // Check return code of the solver: Unconverged, Failed, or OK
  switch (solverRet) {

  // UNCONVERGED
  case Belos::Unconverged:
    if (verbose) 
      cout << "Belos::BlockCGSolMgr::solve() did not converge!" << endl;
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return 0;
    break;

  // CONVERGED
  case Belos::Converged:
    if (verbose) 
      cout << "Belos::BlockCGSolMgr::solve() converged!" << endl;
    break;
  }

  // Test residuals
  Epetra_MultiVector R( B->Map(), blocksize );

  // R = A*X
  A->Apply( *X, R );

  // R -= B 
  MVT::MvAddMv( -1.0, *B, 1.0, R, R );

  // Compute the 2-norm of each vector in the MultiVector
  // and store them to a std::vector<double>
  std::vector<double> normR(blocksize), normB(blocksize);
  MVT::MvNorm( R, normR );
  MVT::MvNorm( *B, normB );

  // Output results to screen
  if(verbose) {
    cout << scientific << setprecision(6) << showpoint;
    cout << "******************************************************\n"
         << "           Results (outside of linear solver)           \n" 
         << "------------------------------------------------------\n"
         << "  Linear System\t\tRelative Residual\n"
         << "------------------------------------------------------\n";
    for( int i=0 ; i<blocksize ; ++i ) {
      cout << "  " << i+1 << "\t\t\t" << normR[i]/normB[i] << endl;
    }
    cout << "******************************************************\n" << endl;
  }


#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return(EXIT_SUCCESS);
}

#else

#include <stdlib.h>
#include <stdio.h>
#ifdef HAVE_MPI
#include "mpi.h"
#endif

int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
#endif

  puts("Please configure Didasko with:");
  puts("--enable-epetra");
  puts("--enable-teuchos");
  puts("--enable-triutils");
  puts("--enable-anasazi");

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  return 0;
}

#endif

