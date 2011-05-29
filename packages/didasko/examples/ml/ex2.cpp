
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

// Use of ML as a black-box smoothed aggregation preconditioner

#include "Didasko_ConfigDefs.h"
#if defined(HAVE_DIDASKO_EPETRA) && defined(HAVE_DIDASKO_ML) && defined(HAVE_DIDASKO_TRIUTILS) && defined(HAVE_DIDASKO_AZTECOO) && defined(HAVE_DIDASKO_TEUCHOS)

#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Time.h"
#include "AztecOO.h"

// includes required by ML
#include "ml_epetra_preconditioner.h"

#include "Trilinos_Util_CrsMatrixGallery.h"

using namespace Teuchos;
using namespace Trilinos_Util;

#include <iostream>

int main(int argc, char *argv[])
{
  
#ifdef EPETRA_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif
  
  Epetra_Time Time(Comm);

  // initialize an Gallery object
  CrsMatrixGallery Gallery("laplace_3d", Comm);
  Gallery.Set("problem_size", 1000);

  // retrive pointers to matrix and linear problem
  Epetra_RowMatrix * A = Gallery.GetMatrix();
  Epetra_LinearProblem * Problem = Gallery.GetLinearProblem();
  
  // Construct a solver object for this problem
  AztecOO solver(*Problem);
  
  // create the preconditioner object and compute hierarchy
  ML_Epetra::MultiLevelPreconditioner * MLPrec = 
    new ML_Epetra::MultiLevelPreconditioner(*A, true);

  // tell AztecOO to use this preconditioner, then solve
  solver.SetPrecOperator(MLPrec);
 
  solver.SetAztecOption(AZ_solver, AZ_gmres_condnum);
  solver.SetAztecOption(AZ_output, 32);
  solver.SetAztecOption(AZ_kspace, 160);
   
  int Niters = 500;
  solver.Iterate(Niters, 1e-12);

  // print out some information about the preconditioner
  if( Comm.MyPID() == 0 ) cout << MLPrec->GetOutputList();
  
  delete MLPrec;
  
  // compute the real residual

  double residual, diff;
  Gallery.ComputeResidual(&residual);
  Gallery.ComputeDiffBetweenStartingAndExactSolutions(&diff);
  
  if( Comm.MyPID()==0 ) {
    cout << "||b-Ax||_2 = " << residual << endl;
    cout << "||x_exact - x||_2 = " << diff << endl;
    cout << "Total Time = " << Time.ElapsedTime() << endl;
  }

  if (residual > 1e-5)
    exit(EXIT_FAILURE);
#ifdef EPETRA_MPI
  MPI_Finalize() ;
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
  MPI_Init(&argc,&argv);
#endif
  
  puts("Please configure Didasko with:\n"
       "--enable-epetra\n"
       "--enable-teuchos\n"
       "--enable-triutils\n"
       "--enable-aztecoo\n"
       "--enable-ml");

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  
  return(EXIT_SUCCESS);
}
#endif
