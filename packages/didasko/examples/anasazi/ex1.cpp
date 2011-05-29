
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

// Compute the smallest magnitude eigenvalues and corresponding eigenvectors using
// Anasazi::BlockKrylovSchur

#include "Didasko_ConfigDefs.h"
#if defined(HAVE_DIDASKO_EPETRA) && defined(HAVE_DIDASKO_ANASAZI) && defined(HAVE_DIDASKO_TEUCHOS) && defined(HAVE_DIDASKO_TRIUTILS) 

#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "Epetra_MultiVector.h"
#include "Trilinos_Util_CrsMatrixGallery.h"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziEpetraAdapter.hpp"
#include "AnasaziBlockKrylovSchurSolMgr.hpp"

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
  typedef Anasazi::MultiVecTraits<double, Epetra_MultiVector> MVT;

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
    cout << "Anasazi Example: Block Kyrlov Schur" << endl;
    cout << "Problem info:" << endl;
    cout << Gallery;
    cout << endl;
  }

  // Setup some more problem/solver parameters:

  // Number of desired eigenvalues
  int nev = 4;

  // Block size
  int blocksize = 4;


  // Get a pointer to the system matrix, inside of a Teuchos::RCP
  // The Teuchos::RCP is a reference counting pointer that handles
  // garbage collection for us, so that we can perform memory allocation without
  // having to worry about freeing memory manually.
  Teuchos::RCP<Epetra_CrsMatrix> A = Teuchos::rcp( Gallery.GetMatrix(), false );

  // Create an Anasazi MultiVector, based on Epetra MultiVector
  const Epetra_Map * Map = &(A->RowMap());
  Teuchos::RCP<Epetra_MultiVector> ivec = 
    Teuchos::rcp( new Epetra_MultiVector(*Map,blocksize) );

  // Fill it with random numbers
  ivec->Random();


  // Setup the eigenproblem, with the matrix A and the initial vectors ivec
  Teuchos::RCP< Anasazi::BasicEigenproblem<double,MV,OP> > MyProblem = 
    Teuchos::rcp( new Anasazi::BasicEigenproblem<double,MV,OP>(A, ivec) );

  // The 2-D laplacian is symmetric. Specify this in the eigenproblem.
  MyProblem->setHermitian(true);

  // Specify the desired number of eigenvalues
  MyProblem->setNEV( nev );

  // Signal that we are done setting up the eigenvalue problem
  ierr = MyProblem->setProblem();

  // Check the return from setProblem(). If this is true, there was an
  // error. This probably means we did not specify enough information for
  // the eigenproblem.
  assert(ierr == true);

  // Specify the verbosity level. Options include:
  // Anasazi::Errors 
  //   This option is always set
  // Anasazi::Warnings 
  //   Warnings (less severe than errors)
  // Anasazi::IterationDetails 
  //   Details at each iteration, such as the current eigenvalues
  // Anasazi::OrthoDetails 
  //   Details about orthogonality
  // Anasazi::TimingDetails
  //   A summary of the timing info for the solve() routine
  // Anasazi::FinalSummary 
  //   A final summary 
  // Anasazi::Debug 
  //   Debugging information
  int verbosity = Anasazi::Warnings + Anasazi::Errors + Anasazi::FinalSummary + Anasazi::TimingDetails;

  // Choose which eigenvalues to compute
  // Choices are:
  // LM - target the largest magnitude  [default]
  // SM - target the smallest magnitude 
  // LR - target the largest real 
  // SR - target the smallest real 
  // LI - target the largest imaginary
  // SI - target the smallest imaginary
  std::string which("SM");

  // Create the parameter list for the eigensolver
  Teuchos::ParameterList MyPL;
  MyPL.set( "Verbosity", verbosity );
  MyPL.set( "Which", which );
  MyPL.set( "Block Size", blocksize );
  MyPL.set( "Num Blocks", 20 );
  MyPL.set( "Maximum Restarts", 100 );
  MyPL.set( "Convergence Tolerance", 1.0e-8 );

  // Create the Block Krylov Schur solver
  // This takes as inputs the eigenvalue problem and the solver parameters
  Anasazi::BlockKrylovSchurSolMgr<double,MV,OP> 
    MyBlockKrylovSchur(MyProblem, MyPL );

  // Solve the eigenvalue problem, and save the return code
  Anasazi::ReturnType solverreturn = MyBlockKrylovSchur.solve();

  // Check return code of the solver: Unconverged, Failed, or OK
  switch (solverreturn) {

  // UNCONVERGED
  case Anasazi::Unconverged:
    if (verbose) 
      cout << "Anasazi::BlockKrylovSchur::solve() did not converge!" << endl;
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return 0;
    break;

  // CONVERGED
  case Anasazi::Converged:
    if (verbose) 
      cout << "Anasazi::BlockKrylovSchur::solve() converged!" << endl;
    break;
  }

  // Get eigensolution struct
  Anasazi::Eigensolution<double, Epetra_MultiVector> sol = MyProblem->getSolution();

  // Get the number of eigenpairs returned
  int numev = sol.numVecs;

  // Get eigenvectors
  Teuchos::RCP<Epetra_MultiVector> evecs = sol.Evecs;

  // Get eigenvalues
  std::vector<Anasazi::Value<double> > evals = sol.Evals;

  // Test residuals
  // Generate a (numev x numev) dense matrix for the eigenvalues...
  // This matrix is automatically initialized to zero
  Teuchos::SerialDenseMatrix<int, double> D(numev, numev);

  // Add the eigenvalues on the diagonals (only the real part since problem is Hermitian)
  for (int i=0; i<numev; i++) {
    D(i,i) = evals[i].realpart;
  }

  // Generate a multivector for the product of the matrix and the eigenvectors
  Epetra_MultiVector R(*Map, numev); 

  // R = A*evecs
  A->Apply( *evecs, R );

  // R -= evecs*D 
  //    = A*evecs - evecs*D
  MVT::MvTimesMatAddMv( -1.0, *evecs, D, 1.0, R );

  // Compute the 2-norm of each vector in the MultiVector
  // and store them to a std::vector<double>
  std::vector<double> normR(numev);
  MVT::MvNorm( R, normR );

  // Output results to screen
  if(verbose) {
    cout << scientific << showpos << setprecision(6) << showpoint;
    cout << "******************************************************\n"
         << "           Results (outside of eigensolver)           \n" 
         << "------------------------------------------------------\n"
         << "  evals\t\t\tnormR\n"
         << "------------------------------------------------------\n";
    for( int i=0 ; i<numev ; ++i ) {
      cout << "  " << evals[i].realpart << "\t\t" << normR[i] << endl;
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

