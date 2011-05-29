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

#include <stdio.h>

#include "az_aztec_defs.h"
#include "CAztecOO.h"

#include "CTrilinos_enums.h"
#include "CTrilinos_flex_enums.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "CEpetra_MpiComm.h"
#else
#include "CEpetra_SerialComm.h"
#endif
#include "CEpetra_Comm.h"

#include "CEpetra_Map.h"
#include "CEpetra_BlockMap.h"
#include "CEpetra_Vector.h"
#include "CEpetra_MultiVector.h"
#include "CEpetra_CrsMatrix.h"
#include "CEpetra_LinearProblem.h"

int main(int argc, char *argv[])
{
  CT_Epetra_Map_ID_Flex_t Map;
  CT_Epetra_CrsMatrix_ID_Flex_t A;
  CT_Epetra_Vector_ID_Flex_t x, b;
  CT_Epetra_LinearProblem_ID_t problem;
  CT_AztecOO_ID_t solver;

  int NumMyElements, NumGlobalElements, i, GlobalRow, RowLess1, RowPlus1;
  double negOne, posTwo;

#ifdef HAVE_MPI
  /* Initialize MPI */
  CT_Epetra_MpiComm_ID_Flex_t Comm;
  MPI_Init(&argc, &argv);
  Comm.Epetra_MpiComm = Epetra_MpiComm_Create(MPI_COMM_WORLD);
#else
  CT_Epetra_SerialComm_ID_Flex_t Comm;
  Comm.Epetra_SerialComm = Epetra_SerialComm_Create();
#endif

  NumMyElements = 1000;

  /* Construct a Map that puts same number of equations on each processor */
  Map.Epetra_Map = Epetra_Map_Create_Linear(-1, NumMyElements, 0, Comm.Epetra_Comm);

  NumGlobalElements = Epetra_BlockMap_NumGlobalElements(Map.Epetra_BlockMap);

  /* Create a Epetra_Matrix */
  A.Epetra_CrsMatrix = Epetra_CrsMatrix_Create(CT_Epetra_DataAccess_E_Copy, Map.Epetra_Map, 3, FALSE);
  
  /* Add rows one-at-a-time */
  negOne = -1.0;
  posTwo = 2.0;
  for (i=0; i<NumMyElements; i++) {
    GlobalRow = Epetra_CrsMatrix_GRID(A.Epetra_CrsMatrix, i);
    RowLess1 = GlobalRow - 1;
    RowPlus1 = GlobalRow + 1;

    if (RowLess1 != -1)
      Epetra_CrsMatrix_InsertGlobalValues(A.Epetra_CrsMatrix, GlobalRow, 1, &negOne, &RowLess1);
    if (RowPlus1 != NumGlobalElements)
      Epetra_CrsMatrix_InsertGlobalValues(A.Epetra_CrsMatrix, GlobalRow, 1, &negOne, &RowPlus1);
    Epetra_CrsMatrix_InsertGlobalValues(A.Epetra_CrsMatrix, GlobalRow, 1, &posTwo, &GlobalRow);
  }
  
  /* Finish up */
  Epetra_CrsMatrix_FillComplete(A.Epetra_CrsMatrix, TRUE);

  /* Create x and b vectors */
  x.Epetra_Vector = Epetra_Vector_Create(Map.Epetra_BlockMap, TRUE);
  b.Epetra_Vector = Epetra_Vector_Create(Map.Epetra_BlockMap, TRUE);
  Epetra_MultiVector_Random(b.Epetra_MultiVector);

  /* Create Linear Problem */
  problem = Epetra_LinearProblem_Create_FromMatrix(A.Epetra_RowMatrix, x.Epetra_MultiVector, b.Epetra_MultiVector);

  /* Create AztecOO instance */
  solver = AztecOO_Create_FromLinearProblem(problem);

  AztecOO_SetAztecOption(solver, AZ_precond, AZ_Jacobi);
  AztecOO_Iterate_Current(solver, 1000, 1.0E-8);

  printf("Solver performed %d iterations.\n", AztecOO_NumIters(solver));
  printf("Norm of true residual = %g\n", AztecOO_TrueResidual(solver));

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif
  return 0;
}

