
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

// Compute the lowest eigenvalues and corresponding eigenvectors using
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
  typedef Anasazi::MultiVecTraits<double, MV> MVT;
  typedef Anasazi::OperatorTraits<double, MV, OP> OPT;

  bool ierr;

  // Get our processor ID (0 if running serial)
  int MyPID = Comm.MyPID();

  // Verbose flag: only processor 0 should print
  bool verbose = (MyPID==0);

  // Initialize a Gallery object, from which we will select a matrix
  // Select recirc_2d, of order 100
  Trilinos_Util::CrsMatrixGallery Gallery("recirc_2d", Comm);
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

  // The problem is non-symmetric. Specify this in the eigenproblem.
  MyProblem->setHermitian(false);

  // Specify the desired number of eigenvalues
  MyProblem->setNEV( nev );

  // Signal that we are done setting up the eigenvalue problem
  ierr = MyProblem->setProblem();

  // Check the return from setProblem(). If this is true, there was an
  // error. This probably means we did not specify enough information for
  // the eigenproblem.
  assert(ierr);

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
  int verbosity = Anasazi::Warnings + Anasazi::Errors + Anasazi::FinalSummary;

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
  MyPL.set( "Num Blocks", 5 );
  MyPL.set( "Maximum Restarts", 300 );
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

  // Get the index vector for the eigenvector storage
  std::vector<int> index = sol.index;
  
  // Get eigenvectors
  Teuchos::RCP<Epetra_MultiVector> evecs = sol.Evecs;
  
  // Get eigenvalues
  // Because the matrix is nonsymmetric (and the eigenvalues are potentially 
  // complex), Evals is a vector of length numev where each entry is a pair
  // with a real and imaginary part.
  std::vector<Anasazi::Value<double> > evals = sol.Evals;
  
  // Compute residuals.
  Teuchos::LAPACK<int,double> lapack;
  std::vector<double> normA(numev);
  
  int i=0;
  std::vector<int> curind(1);
  std::vector<double> resnorm(1), tempnrm(1);
  Teuchos::RCP<MV> evecr, eveci, tempAevec;
  Epetra_MultiVector Aevec(*Map,numev);
  
  // Compute A*evecs
  OPT::Apply( *A, *evecs, Aevec );
  
  Teuchos::SerialDenseMatrix<int,double> Breal(1,1), Bimag(1,1);
  while (i<numev) {
    if (index[i]==0) {
      // Get a view of the current eigenvector (evecr)
      curind[0] = i;
      evecr = MVT::CloneViewNonConst( *evecs, curind );
      
      // Get a copy of A*evecr
      tempAevec = MVT::CloneCopy( Aevec, curind );
      
      // Compute A*evecr - lambda*evecr
      Breal(0,0) = evals[i].realpart;
      MVT::MvTimesMatAddMv( -1.0, *evecr, Breal, 1.0, *tempAevec );
      
      // Compute the norm of the residual and increment counter
      MVT::MvNorm( *tempAevec, resnorm );
      normA[i] = resnorm[0]/Teuchos::ScalarTraits<double>::magnitude( evals[i].realpart );
      i++;
    } else {
      // Get a view of the real part of the eigenvector (evecr)
      curind[0] = i;
      evecr = MVT::CloneViewNonConst( *evecs, curind );
      
      // Get a copy of A*evecr
      tempAevec = MVT::CloneCopy( Aevec, curind );
      
      // Get a view of the imaginary part of the eigenvector (eveci)
      curind[0] = i+1;
      eveci = MVT::CloneViewNonConst( *evecs, curind );
      
      // Set the eigenvalue into Breal and Bimag
      Breal(0,0) = evals[i].realpart;
      Bimag(0,0) = evals[i].imagpart;
      
      // Compute A*evecr - evecr*lambdar + eveci*lambdai
      MVT::MvTimesMatAddMv( -1.0, *evecr, Breal, 1.0, *tempAevec );
      MVT::MvTimesMatAddMv( 1.0, *eveci, Bimag, 1.0, *tempAevec );
      MVT::MvNorm( *tempAevec, tempnrm );
      
      // Get a copy of A*eveci
      tempAevec = MVT::CloneCopy( Aevec, curind );
      
      // Compute A*eveci - eveci*lambdar - evecr*lambdai
      MVT::MvTimesMatAddMv( -1.0, *evecr, Bimag, 1.0, *tempAevec );
      MVT::MvTimesMatAddMv( -1.0, *eveci, Breal, 1.0, *tempAevec );
      MVT::MvNorm( *tempAevec, resnorm );
      
      // Compute the norms and scale by magnitude of eigenvalue
      normA[i] = lapack.LAPY2( tempnrm[i], resnorm[i] ) /
	lapack.LAPY2( evals[i].realpart, evals[i].imagpart );
      normA[i+1] = normA[i];
      
      i=i+2;
    }
  }
  
  // Output results to screen
  if(verbose) {
    cout << scientific << showpos << setprecision(6) << showpoint;
    cout << "******************************************************\n"
         << "           Results (outside of eigensolver)           \n" 
         << "------------------------------------------------------\n"
         << "  real(evals)\t\timag(evals)\tnormR\n"
         << "------------------------------------------------------\n";
    for(int i=0; i<numev; ++i) {
      cout << "  " << evals[i].realpart << "\t\t" << evals[i].imagpart
           << "\t" << normA[i] << endl;
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

