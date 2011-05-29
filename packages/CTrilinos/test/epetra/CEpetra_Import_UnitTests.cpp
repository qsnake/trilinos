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
#include "CEpetra_Comm.h"
#include "CEpetra_BlockMap.h"
#include "CEpetra_Map.h"
#include "CEpetra_Distributor.h"
#include "CEpetra_Object.h"
#include "Epetra_Import.h"
#include "CEpetra_Import.h"
#include "CEpetra_Import_Cpp.hpp"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_flex_enums.h"
#include "CTrilinos_exceptions.hpp"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_test_utils.hpp"

#include "CTrilinos_UnitTestHelpers.hpp"
#include "Teuchos_UnitTestHarness.hpp"


namespace {


/**********************************************************************
CT_Epetra_Import_ID_t Epetra_Import_Create ( 
  CT_Epetra_BlockMap_ID_t TargetMapID, 
  CT_Epetra_BlockMap_ID_t SourceMapID );
 **********************************************************************/

TEUCHOS_UNIT_TEST ( Epetra_Import , Create )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Set up communication */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int NumProc = Epetra_Comm_NumProc(CommID));
  ECHO(int MyPID = Epetra_Comm_MyPID(CommID));

  /* Create the source map */
  ECHO(int IndexBase = 0);
  ECHO(const int NumMyElements = 4);
  ECHO(int NumGlobalElements = NumMyElements * NumProc);
  ECHO(int off = NumMyElements*MyPID);
  int MyGlobalElements[NumMyElements] = {0+off, 1+off, 2+off, 3+off};
  ECHO(CT_Epetra_Map_ID_Flex_t srcID);
  ECHO(srcID.Epetra_Map = Epetra_Map_Create_Arbitrary(
       NumGlobalElements, NumMyElements, MyGlobalElements, IndexBase, CommID));

  /* Create the target map */
  ECHO(int off2 = NumMyElements*(NumProc-MyPID));
  int GetMyGlobalElements[NumMyElements] = {off2-1, off2-2, off2-3, off2-4};
  ECHO(CT_Epetra_Map_ID_Flex_t tarID);
  ECHO(tarID.Epetra_Map = Epetra_Map_Create_Arbitrary(
       NumGlobalElements, NumMyElements, GetMyGlobalElements, IndexBase, CommID));

  /* Create an importer */
  ECHO(CT_Epetra_Import_ID_t selfID = Epetra_Import_Create(tarID.Epetra_BlockMap, srcID.Epetra_BlockMap));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(selfID.table, CT_Epetra_Import_ID);
  TEST_EQUALITY_CONST(selfID.index, 0);
  ECHO(Teuchos::RCP<Epetra_Import> rcp1 = CEpetra::getImport(selfID));
  TEST_EQUALITY_CONST(rcp1.is_null(), false);
}

/**********************************************************************
CT_Epetra_Import_ID_t Epetra_Import_Duplicate ( 
  CT_Epetra_Import_ID_t ImporterID );
 **********************************************************************/

TEUCHOS_UNIT_TEST ( Epetra_Import , Duplicate )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Set up communication */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int NumProc = Epetra_Comm_NumProc(CommID));
  ECHO(int MyPID = Epetra_Comm_MyPID(CommID));

  /* Create the source map */
  ECHO(int IndexBase = 0);
  ECHO(const int NumMyElements = 4);
  ECHO(int NumGlobalElements = NumMyElements * NumProc);
  ECHO(int off = NumMyElements*MyPID);
  int MyGlobalElements[NumMyElements] = {0+off, 1+off, 2+off, 3+off};
  ECHO(CT_Epetra_Map_ID_Flex_t srcID);
  ECHO(srcID.Epetra_Map = Epetra_Map_Create_Arbitrary(
       NumGlobalElements, NumMyElements, MyGlobalElements, IndexBase, CommID));

  /* Create the target map */
  ECHO(int off2 = NumMyElements*(NumProc-MyPID));
  int GetMyGlobalElements[NumMyElements] = {off2-1, off2-2, off2-3, off2-4};
  ECHO(CT_Epetra_Map_ID_Flex_t tarID);
  ECHO(tarID.Epetra_Map = Epetra_Map_Create_Arbitrary(
       NumGlobalElements, NumMyElements, GetMyGlobalElements, IndexBase, CommID));

  /* Create an importer */
  ECHO(CT_Epetra_Import_ID_t selfID = Epetra_Import_Create(tarID.Epetra_BlockMap, srcID.Epetra_BlockMap));

  /* Duplicate it */
  ECHO(CT_Epetra_Import_ID_t dupID = Epetra_Import_Duplicate(selfID));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(dupID.table, CT_Epetra_Import_ID);
  TEST_EQUALITY_CONST(dupID.index, 1);
  TEST_EQUALITY_CONST(CTrilinos::isSameObject(selfID, dupID), false);
}

/**********************************************************************
void Epetra_Import_Destroy ( CT_Epetra_Import_ID_t * selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST ( Epetra_Import , Destroy )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Set up communication */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int NumProc = Epetra_Comm_NumProc(CommID));
  ECHO(int MyPID = Epetra_Comm_MyPID(CommID));

  /* Create the source map */
  ECHO(int IndexBase = 0);
  ECHO(const int NumMyElements = 4);
  ECHO(int NumGlobalElements = NumMyElements * NumProc);
  ECHO(int off = NumMyElements*MyPID);
  int MyGlobalElements[NumMyElements] = {0+off, 1+off, 2+off, 3+off};
  ECHO(CT_Epetra_Map_ID_Flex_t srcID);
  ECHO(srcID.Epetra_Map = Epetra_Map_Create_Arbitrary(
       NumGlobalElements, NumMyElements, MyGlobalElements, IndexBase, CommID));

  /* Create the target map */
  ECHO(int off2 = NumMyElements*(NumProc-MyPID));
  int GetMyGlobalElements[NumMyElements] = {off2-1, off2-2, off2-3, off2-4};
  ECHO(CT_Epetra_Map_ID_Flex_t tarID);
  ECHO(tarID.Epetra_Map = Epetra_Map_Create_Arbitrary(
       NumGlobalElements, NumMyElements, GetMyGlobalElements, IndexBase, CommID));

  /* Create an importer */
  ECHO(CT_Epetra_Import_ID_t selfID = Epetra_Import_Create(tarID.Epetra_BlockMap, srcID.Epetra_BlockMap));

  ECHO(Epetra_Import_Destroy(&selfID));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(selfID.table, CT_Invalid_ID);
  TEST_EQUALITY_CONST(selfID.index, -1);
}

/**********************************************************************
CT_Epetra_BlockMap_ID_t Epetra_Import_SourceMap ( 
  CT_Epetra_Import_ID_t selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST ( Epetra_Import , SourceMap )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Set up communication */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int NumProc = Epetra_Comm_NumProc(CommID));
  ECHO(int MyPID = Epetra_Comm_MyPID(CommID));

  /* Create the source map */
  ECHO(int IndexBase = 0);
  ECHO(const int NumMyElements = 4);
  ECHO(int NumGlobalElements = NumMyElements * NumProc);
  ECHO(int off = NumMyElements*MyPID);
  int MyGlobalElements[NumMyElements] = {0+off, 1+off, 2+off, 3+off};
  ECHO(CT_Epetra_Map_ID_Flex_t srcID);
  ECHO(srcID.Epetra_Map = Epetra_Map_Create_Arbitrary(
       NumGlobalElements, NumMyElements, MyGlobalElements, IndexBase, CommID));

  /* Create the target map */
  ECHO(int off2 = NumMyElements*(NumProc-MyPID));
  int GetMyGlobalElements[NumMyElements] = {off2-1, off2-2, off2-3, off2-4};
  ECHO(CT_Epetra_Map_ID_Flex_t tarID);
  ECHO(tarID.Epetra_Map = Epetra_Map_Create_Arbitrary(
       NumGlobalElements, NumMyElements, GetMyGlobalElements, IndexBase, CommID));

  /* Create an importer */
  ECHO(CT_Epetra_Import_ID_t selfID = Epetra_Import_Create(tarID.Epetra_BlockMap, srcID.Epetra_BlockMap));

  /* Get the source map using wrapper */
  ECHO(CT_Epetra_BlockMap_ID_t bsrcID2 = Epetra_Import_SourceMap(selfID));
  TEST_EQUALITY(Epetra_BlockMap_NumGlobalElements(bsrcID2), NumGlobalElements);
}

/**********************************************************************
CT_Epetra_BlockMap_ID_t Epetra_Import_TargetMap ( 
  CT_Epetra_Import_ID_t selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST ( Epetra_Import , TargetMap )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Set up communication */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int NumProc = Epetra_Comm_NumProc(CommID));
  ECHO(int MyPID = Epetra_Comm_MyPID(CommID));

  /* Create the source map */
  ECHO(int IndexBase = 0);
  ECHO(const int NumMyElements = 4);
  ECHO(int NumGlobalElements = NumMyElements * NumProc);
  ECHO(int off = NumMyElements*MyPID);
  int MyGlobalElements[NumMyElements] = {0+off, 1+off, 2+off, 3+off};
  ECHO(CT_Epetra_Map_ID_Flex_t srcID);
  ECHO(srcID.Epetra_Map = Epetra_Map_Create_Arbitrary(
       NumGlobalElements, NumMyElements, MyGlobalElements, IndexBase, CommID));

  /* Create the target map */
  ECHO(int off2 = NumMyElements*(NumProc-MyPID));
  int GetMyGlobalElements[NumMyElements] = {off2-1, off2-2, off2-3, off2-4};
  ECHO(CT_Epetra_Map_ID_Flex_t tarID);
  ECHO(tarID.Epetra_Map = Epetra_Map_Create_Arbitrary(
       NumGlobalElements, NumMyElements, GetMyGlobalElements, IndexBase, CommID));

  /* Create an importer */
  ECHO(CT_Epetra_Import_ID_t selfID = Epetra_Import_Create(tarID.Epetra_BlockMap, srcID.Epetra_BlockMap));

  /* Get the source map using wrapper */
  ECHO(CT_Epetra_BlockMap_ID_t btarID2 = Epetra_Import_TargetMap(selfID));
  TEST_EQUALITY(Epetra_BlockMap_NumGlobalElements(btarID2), NumGlobalElements);
}

/**********************************************************************
CT_Epetra_Distributor_ID_t Epetra_Import_Distributor ( 
  CT_Epetra_Import_ID_t selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST ( Epetra_Import , Distributor )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Set up communication */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int NumProc = Epetra_Comm_NumProc(CommID));
  ECHO(int MyPID = Epetra_Comm_MyPID(CommID));

  /* Create the source map */
  ECHO(int IndexBase = 0);
  ECHO(const int NumMyElements = 4);
  ECHO(int NumGlobalElements = NumMyElements * NumProc);
  ECHO(int off = NumMyElements*MyPID);
  int MyGlobalElements[NumMyElements] = {0+off, 1+off, 2+off, 3+off};
  ECHO(CT_Epetra_Map_ID_Flex_t srcID);
  ECHO(srcID.Epetra_Map = Epetra_Map_Create_Arbitrary(
       NumGlobalElements, NumMyElements, MyGlobalElements, IndexBase, CommID));

  /* Create the target map */
  ECHO(int off2 = NumMyElements*(NumProc-MyPID));
  int GetMyGlobalElements[NumMyElements] = {off2-1, off2-2, off2-3, off2-4};
  ECHO(CT_Epetra_Map_ID_Flex_t tarID);
  ECHO(tarID.Epetra_Map = Epetra_Map_Create_Arbitrary(
       NumGlobalElements, NumMyElements, GetMyGlobalElements, IndexBase, CommID));

  /* Create an importer */
  ECHO(CT_Epetra_Import_ID_t selfID = Epetra_Import_Create(tarID.Epetra_BlockMap, srcID.Epetra_BlockMap));

  /* No distributor will exist if not distributedglobal */
  if (Epetra_BlockMap_DistributedGlobal(srcID.Epetra_BlockMap) == TRUE) {
    ECHO(CT_Epetra_Distributor_ID_t dID = Epetra_Import_Distributor(selfID));

    /* Now check the result of the call to the wrapper function */
    TEST_EQUALITY(dID.table, CT_Epetra_Distributor_ID);
    TEST_EQUALITY_CONST(dID.index, 0);
    TEST_EQUALITY_CONST(CEpetra::getDistributor(dID).is_null(), false);
  }
}

/**********************************************************************/

//
// Template Instantiations
//


#ifdef TEUCHOS_DEBUG

#  define DEBUG_UNIT_TEST_GROUP( T ) \

#else

#  define DEBUG_UNIT_TEST_GROUP( T )

#endif


#define UNIT_TEST_GROUP( T ) \
  DEBUG_UNIT_TEST_GROUP( T )


} // namespace

