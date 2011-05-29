/*! \@HEADER */
/*
************************************************************************

                CTrilinos:  C interface to Trilinos
                Copyright (2009) Sandia Corporation

Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
license for use of this work by or on behalf of the U.S. Government.

This library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation; either version 2.1 of the
License, or (at your option) any later version.

This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
USA
Questions? Contact M. Nicole Lemaster (mnlemas\@sandia.gov)

************************************************************************
*/
/*! \@HEADER */


#include "CTrilinos_config.h"

#include <stdlib.h>
#include <stdio.h>

#include "CTrilinos_enums.h"
#include "CTrilinos_flex_enums.h"
#ifdef HAVE_MPI
#include "CEpetra_MpiComm.h"
#include "mpi.h"
#else
#include "CEpetra_SerialComm.h"
#endif
#include "CEpetra_Comm.h"
#include "CEpetra_Map.h"
#include "CEpetra_BlockMap.h"
#include "CEpetra_Vector.h"
#include "CEpetra_MultiVector.h"
#include "CEpetra_CrsMatrix.h"
#include "CEpetra_RowMatrix.h"
#include "CEpetra_LinearProblem.h"
#include "CTeuchos_ParameterList.h"
#include "CGaleri_Maps.h"
#include "CGaleri_CrsMatrices.h"
#include "CGaleri_Utils.h"

/* ===========
 * main driver
 * =========== */

int main(int argc, char* argv[])
{
  int NumProc, MyPID;
  double ResidualNorm;

  CT_Epetra_Map_ID_Flex_t Map;
  CT_Epetra_CrsMatrix_ID_Flex_t Matrix;
  CT_Teuchos_ParameterList_ID_t GaleriList;
  CT_Epetra_Vector_ID_Flex_t ExactSolution, LHS, RHS;
  CT_Epetra_LinearProblem_ID_t Problem;

#ifdef HAVE_MPI
  CT_Epetra_MpiComm_ID_Flex_t Comm;
  MPI_Init(&argc, &argv);
  Comm.Epetra_MpiComm = Epetra_MpiComm_Create(MPI_COMM_WORLD);
#else
  CT_Epetra_SerialComm_ID_Flex_t Comm;
  Comm.Epetra_SerialComm = Epetra_SerialComm_Create();
#endif
  NumProc = Epetra_Comm_NumProc(Comm.Epetra_Comm);
  MyPID = Epetra_Comm_MyPID(Comm.Epetra_Comm);

  /* Here we create the linear problem
   *
   *   Matrix * LHS = RHS
   *
   * with Matrix arising from a 5-point formula discretization. */

  GaleriList = Teuchos_ParameterList_Create();
  
  /* dimension of the problem is nx x ny */
  Teuchos_ParameterList_set_int(GaleriList, "nx", 10 * NumProc, "");
  Teuchos_ParameterList_set_int(GaleriList, "ny", 10, "");
  /* total number of processors is mx x my */
  Teuchos_ParameterList_set_int(GaleriList, "mx", NumProc, "");
  Teuchos_ParameterList_set_int(GaleriList, "my", 1, "");

  Map.Epetra_Map = Galeri_Maps_CreateMap("Cartesian2D", Comm.Epetra_Comm, GaleriList);

  Matrix.Epetra_CrsMatrix = Galeri_CrsMatrices_CreateCrsMatrix("Laplace2D",
      Map.Epetra_Map, GaleriList);

  ExactSolution.Epetra_Vector = Epetra_Vector_Create(Map.Epetra_BlockMap, TRUE);
  Epetra_MultiVector_Random(ExactSolution.Epetra_MultiVector);

  LHS.Epetra_Vector = Epetra_Vector_Create(Map.Epetra_BlockMap, TRUE);
  RHS.Epetra_Vector = Epetra_Vector_Create(Map.Epetra_BlockMap, TRUE);

  Epetra_CrsMatrix_Multiply_Vector(Matrix.Epetra_CrsMatrix, FALSE,
      ExactSolution.Epetra_Vector, RHS.Epetra_Vector);

  Problem = Epetra_LinearProblem_Create_FromMatrix(Matrix.Epetra_RowMatrix,
      LHS.Epetra_MultiVector, RHS.Epetra_MultiVector);

  /* at this point any object that understand Epetra_LinearProblem can be
   * used, for example AztecOO, Amesos. IFPACK and ML can be used to define a
   * preconditioner for Matrix. Here we use a simple solver, based on
   * LAPACK, that is meant for simple testing only. */
    
  Galeri_Utils_Solve_LinearProblem(Problem);

  /* and we compute the norm of the true residual. */
  ResidualNorm = Galeri_Utils_ComputeNorm_Matrix(Matrix.Epetra_RowMatrix,
      LHS.Epetra_MultiVector, RHS.Epetra_MultiVector);

  if (MyPID == 0)
    printf("%g\n", ResidualNorm);

  Epetra_BlockMap_Destroy(&Map.Epetra_BlockMap);
  Epetra_CrsMatrix_Destroy(&Matrix.Epetra_CrsMatrix);
  Epetra_Vector_Destroy(&LHS.Epetra_Vector);
  Epetra_Vector_Destroy(&RHS.Epetra_Vector);
  Epetra_Vector_Destroy(&ExactSolution.Epetra_Vector);

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return(EXIT_SUCCESS);
}
