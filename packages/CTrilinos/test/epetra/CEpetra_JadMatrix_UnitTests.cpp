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
#include "CEpetra_Map.h"
#include "CEpetra_MultiVector.h"
#include "CEpetra_RowMatrix.h"
#include "CEpetra_CrsMatrix.h"
#include "Epetra_JadMatrix.h"
#include "CEpetra_JadMatrix.h"
#include "CEpetra_JadMatrix_Cpp.hpp"
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
CT_Epetra_JadMatrix_ID_t Epetra_JadMatrix_Create ( 
  CT_Epetra_RowMatrix_ID_t MatrixID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_JadMatrix , Create )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(int NumGlobalElements = 4);
  ECHO(int IndexBase = 0);
  ECHO(CT_Epetra_Map_ID_t MapID = Epetra_Map_Create(NumGlobalElements, IndexBase, CommID));

  /* Create the source matrix */
  ECHO(int NumIndicesPerRow = 4);
  ECHO(CT_Epetra_DataAccess_E_t CV = CT_Epetra_DataAccess_E_Copy);
  ECHO(CT_Epetra_CrsMatrix_ID_Flex_t crsID);
  ECHO(crsID.Epetra_CrsMatrix = Epetra_CrsMatrix_Create(
       CV, MapID, NumIndicesPerRow, FALSE));

  /* Initialize the source matrix */
  ECHO(double val = 1.0);
  ECHO(int ret = Epetra_CrsMatrix_PutScalar(crsID.Epetra_CrsMatrix, val));
  TEST_EQUALITY(ret, 0);
  ECHO(ret = Epetra_CrsMatrix_FillComplete(crsID.Epetra_CrsMatrix, TRUE));
  TEST_EQUALITY(ret, 0);

  /* Create a JadMatrix from the RowMatrix */
  ECHO(CT_Epetra_JadMatrix_ID_t selfID = Epetra_JadMatrix_Create(crsID.Epetra_RowMatrix));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(selfID.table, CT_Epetra_JadMatrix_ID);
}

/**********************************************************************
void Epetra_JadMatrix_Destroy ( CT_Epetra_JadMatrix_ID_t * selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_JadMatrix_UpdateValues ( 
  CT_Epetra_JadMatrix_ID_t selfID, CT_Epetra_RowMatrix_ID_t MatrixID, 
  boolean CheckStructure );
 **********************************************************************/

/**********************************************************************
int Epetra_JadMatrix_ExtractMyRowCopy ( 
  CT_Epetra_JadMatrix_ID_t selfID, int MyRow, int Length, 
  int * NumEntries, double * Values, int * Indices );
 **********************************************************************/

/**********************************************************************
int Epetra_JadMatrix_ExtractMyEntryView ( 
  CT_Epetra_JadMatrix_ID_t selfID, int CurEntry, double * * Value, 
  int * RowIndex, int * ColIndex );
 **********************************************************************/

/**********************************************************************
int Epetra_JadMatrix_ExtractMyEntryView_Const ( 
  CT_Epetra_JadMatrix_ID_t selfID, int CurEntry, 
  double const ** Value, int * RowIndex, int * ColIndex );
 **********************************************************************/

/**********************************************************************
int Epetra_JadMatrix_NumMyRowEntries ( 
  CT_Epetra_JadMatrix_ID_t selfID, int MyRow, int * NumEntries );
 **********************************************************************/

/**********************************************************************
int Epetra_JadMatrix_Multiply ( 
  CT_Epetra_JadMatrix_ID_t selfID, boolean TransA, 
  CT_Epetra_MultiVector_ID_t XID, CT_Epetra_MultiVector_ID_t YID );
 **********************************************************************/

/**********************************************************************
int Epetra_JadMatrix_Solve ( 
  CT_Epetra_JadMatrix_ID_t selfID, boolean Upper, boolean Trans, 
  boolean UnitDiagonal, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t YID );
 **********************************************************************/

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

