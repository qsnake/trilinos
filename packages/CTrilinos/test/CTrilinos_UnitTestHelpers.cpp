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

#include <assert.h>

#include "CTrilinos_enums.h"
#include "CTrilinos_table_man.h"

#include "CEpetra_Comm_Cpp.hpp"
#include "CEpetra_Comm.h"
#include "CEpetra_SerialComm_Cpp.hpp"
#include "CEpetra_SerialComm.h"
#ifdef HAVE_MPI
#include "CEpetra_MpiComm_Cpp.hpp"
#include "CEpetra_MpiComm.h"
#endif /* HAVE_MPI */
#include "CEpetra_BlockMap_Cpp.hpp"
#include "CEpetra_BlockMap.h"
#include "CEpetra_Map_Cpp.hpp"
#include "CEpetra_Map.h"
#include "CEpetra_Object_Cpp.hpp"
#include "CEpetra_Distributor_Cpp.hpp"
#include "CEpetra_Directory_Cpp.hpp"
#include "CEpetra_Import_Cpp.hpp"
#include "CEpetra_Import.h"
#include "CEpetra_Export_Cpp.hpp"
#include "CEpetra_OffsetIndex_Cpp.hpp"
#include "CEpetra_DistObject_Cpp.hpp"
#include "CEpetra_SrcDistObject_Cpp.hpp"
#include "CEpetra_BLAS_Cpp.hpp"
#include "CEpetra_Flops_Cpp.hpp"
#include "CEpetra_CompObject_Cpp.hpp"
#include "CEpetra_MultiVector_Cpp.hpp"
#include "CEpetra_Vector_Cpp.hpp"
#include "CEpetra_CrsGraph_Cpp.hpp"
#include "CEpetra_CrsMatrix_Cpp.hpp"
#include "CEpetra_Operator_Cpp.hpp"
#include "CEpetra_RowMatrix_Cpp.hpp"
#include "CEpetra_Time_Cpp.hpp"
#include "CEpetra_JadMatrix_Cpp.hpp"
#include "CEpetra_LinearProblem_Cpp.hpp"
#include "CEpetra_LAPACK_Cpp.hpp"
#include "CAmesos_Cpp.hpp"
#include "CAmesos_BaseSolver_Cpp.hpp"
#include "CTeuchos_any_Cpp.hpp"
#include "CTeuchos_CommandLineProcessor_Cpp.hpp"
#include "CTeuchos_ParameterEntry_Cpp.hpp"
#include "CTeuchos_ParameterList_Cpp.hpp"

#include "CTrilinos_UnitTestHelpers.hpp"

void CEpetra_Test_CleanSlate()
{
  CEpetra::purgeComm();
  CEpetra::purgeSerialComm();
#ifdef HAVE_MPI
  CEpetra::purgeMpiComm();
#endif /* HAVE_MPI */
  CEpetra::purgeBlockMap();
  CEpetra::purgeMap();
  CEpetra::purgeObject();
  CEpetra::purgeDistributor();
  CEpetra::purgeDirectory();
  CEpetra::purgeImport();
  CEpetra::purgeExport();
  CEpetra::purgeOffsetIndex();
  CEpetra::purgeDistObject();
  CEpetra::purgeSrcDistObject();
  CEpetra::purgeBLAS();
  CEpetra::purgeFlops();
  CEpetra::purgeCompObject();
  CEpetra::purgeMultiVector();
  CEpetra::purgeVector();
  CEpetra::purgeCrsGraph();
  CEpetra::purgeCrsMatrix();
  CEpetra::purgeOperator();
  CEpetra::purgeRowMatrix();
  CEpetra::purgeTime();
  CEpetra::purgeJadMatrix();
  CEpetra::purgeLinearProblem();
  CEpetra::purgeLAPACK();
#ifdef HAVE_CTRILINOS_AMESOS
  CAmesos::purgeAmesos();
  CAmesos::purgeBaseSolver();
#endif
  CTeuchos::purgeany();
  CTeuchos::purgeCommandLineProcessor();
  CTeuchos::purgeParameterEntry();
  CTeuchos::purgeParameterList();
}

CT_Epetra_Comm_ID_t
UnitTest_Create_Comm()
{
#ifdef EPETRA_MPI
  CT_Epetra_MpiComm_ID_t comm = Epetra_MpiComm_Create(MPI_COMM_WORLD);
  CTrilinos_Universal_ID_t ucomm = Epetra_MpiComm_Generalize(comm);
#else /* EPETRA_MPI */
  CT_Epetra_SerialComm_ID_t comm = Epetra_SerialComm_Create();
  CTrilinos_Universal_ID_t ucomm = Epetra_SerialComm_Generalize(comm);
#endif /* EPETRA_MPI */
  CT_Migrate(&ucomm, CT_Epetra_Comm_ID);
  return Epetra_Comm_Degeneralize(ucomm);
}

CT_Epetra_Import_ID_t
initialize_doxygen_example(CT_Epetra_Comm_ID_t CommID)
{
  int MyPID = Epetra_Comm_MyPID(CommID);
  assert(MyPID >= 0);
  int NumProc = Epetra_Comm_NumProc(CommID);
  assert(NumProc > 0);

  CT_Epetra_Import_ID_t selfID;
  selfID.table = CT_Invalid_ID;
  selfID.index = -1;
  selfID.is_const = FALSE;

  /* Create the source map */
  int IndexBase = 0;
  const int NumMyElements = 3;
  int NumGlobalElements = NumMyElements * NumProc;
  int off = NumMyElements*MyPID;
  int MyGlobalElements[NumMyElements] = {0+off, 1+off, 2+off};
  CT_Epetra_Map_ID_t srcID = Epetra_Map_Create_Arbitrary(
       NumGlobalElements, NumMyElements, MyGlobalElements, IndexBase, CommID);

  /* Migrate it */
  CTrilinos_Universal_ID_t tmpID = Epetra_Map_Generalize(srcID);
  CT_Migrate(&tmpID, CT_Epetra_BlockMap_ID);
  CT_Epetra_BlockMap_ID_t bsrcID = Epetra_BlockMap_Degeneralize(tmpID);
  int els = Epetra_BlockMap_NumMyElements(bsrcID);
  assert(els == NumMyElements);

  /* Create the target map */
  const int NumMyElements2 = 5;
  int NumGlobalElements2 = NumMyElements2 * NumProc;
  int MyGlobalElements2a[NumMyElements2] = {0, 1, 2, 3, 8};
  int MyGlobalElements2b[NumMyElements2] = {2, 3, 4, 5, 6};
  int MyGlobalElements2c[NumMyElements2] = {0, 5, 6, 7, 8};
  int MyGlobalElements2[NumMyElements2];
  for (int i=0; i<NumMyElements2; i++) {
    switch (MyPID) {
    case 0:
      MyGlobalElements2[i] = MyGlobalElements2a[i];
      break;
    case 1:
      MyGlobalElements2[i] = MyGlobalElements2b[i];
      break;
    case 2:
      MyGlobalElements2[i] = MyGlobalElements2c[i];
      break;
    }
  }
  CT_Epetra_Map_ID_t tarID = Epetra_Map_Create_Arbitrary(
       NumGlobalElements2, NumMyElements2, MyGlobalElements2, IndexBase, CommID);

  /* Migrate it */
  tmpID = Epetra_Map_Generalize(tarID);
  CT_Migrate(&tmpID, CT_Epetra_BlockMap_ID);
  CT_Epetra_BlockMap_ID_t btarID = Epetra_BlockMap_Degeneralize(tmpID);
  int els2 = Epetra_BlockMap_NumMyElements(btarID);
  assert(els2 == NumMyElements2);

  /* Create an importer */
  selfID = Epetra_Import_Create(btarID, bsrcID);
  assert(selfID.table == CT_Epetra_Import_ID);
  assert(selfID.index >= 0);

  return selfID;
}

