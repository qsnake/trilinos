
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

// Solve a linear system with AztecOO. 
// The linear system is created using MatrixGallery

#include "Didasko_ConfigDefs.h"
#if defined(HAVE_DIDASKO_EPETRA) && defined(HAVE_DIDASKO_TRIUTILS)

#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "AztecOO.h"
#include "Epetra_CrsMatrix.h"

#include "Trilinos_Util_CrsMatrixGallery.h"

using namespace Trilinos_Util;

int main(int argc, char *argv[])
{
    
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  // initialize an Gallery object
  CrsMatrixGallery Gallery("laplace_2d", Comm);
  Gallery.Set("problem_size",900);

  // Get pointers to the linear problem (containing matrix, LHS and RHS).
  Epetra_LinearProblem * Problem= Gallery.GetLinearProblem();

  // initialize the AztecOO solve object, based on current linear problem
  AztecOO solver(*Problem);

  // here set some AztecOO options:
  // - symmetric problem;
  // - domain decomposition preconditioner
  // - ICC factorization on each subdomain
  solver.SetAztecOption(AZ_solver, AZ_cg_condnum);
  solver.SetAztecOption(AZ_precond, AZ_dom_decomp);
  solver.SetAztecOption(AZ_overlap,0);
  solver.SetAztecOption(AZ_subdomain_solve, AZ_icc);

  // solve the linear system
  solver.Iterate(1550, 1e-12);

  // AztecOO defined a certain number of output parameters, and store them
  // in a double vector called status. 
  double status[AZ_STATUS_SIZE];
  solver.GetAllAztecStatus(status);
  
  // verify that linear system has been solved as required
  double residual, diff;

  Gallery.ComputeResidual(&residual);
  Gallery.ComputeDiffBetweenStartingAndExactSolutions(&diff);
  
  if( Comm.MyPID()==0 ) {
    cout << "||b-Ax||_2 = " << residual << endl;
    cout << "||x_exact - x||_2 = " << diff << endl;
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
  MPI_Init(&argc,&argv);
#endif

  puts("Please configure Didasko with:\n"
       "--enable-epetra\n"
       "--enable-triutils\n"
       "--enable-aztecoo\n");

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  return 0;
}

#endif


