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
#include "CEpetra_Import.h"
#include "CEpetra_Map.h"
#include "CEpetra_MultiVector.h"
#include "CEpetra_Vector.h"
#include "CEpetra_CrsMatrix.h"
#include "Epetra_RowMatrix.h"
#include "CEpetra_RowMatrix.h"
#include "CEpetra_RowMatrix_Cpp.hpp"
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
void Epetra_RowMatrix_Destroy ( CT_Epetra_RowMatrix_ID_t * selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_RowMatrix , Destroy )
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

  ECHO(Epetra_RowMatrix_Destroy(&crsID.Epetra_RowMatrix));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(crsID.Epetra_RowMatrix.table, CT_Invalid_ID);
  TEST_EQUALITY_CONST(crsID.Epetra_RowMatrix.index, -1);
}

/**********************************************************************
int Epetra_RowMatrix_NumMyRowEntries ( 
  CT_Epetra_RowMatrix_ID_t selfID, int MyRow, int * NumEntries );
 **********************************************************************/

/**********************************************************************
int Epetra_RowMatrix_MaxNumEntries ( 
  CT_Epetra_RowMatrix_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_RowMatrix_ExtractMyRowCopy ( 
  CT_Epetra_RowMatrix_ID_t selfID, int MyRow, int Length, 
  int * NumEntries, double * Values, int * Indices );
 **********************************************************************/

/**********************************************************************
int Epetra_RowMatrix_ExtractDiagonalCopy ( 
  CT_Epetra_RowMatrix_ID_t selfID, CT_Epetra_Vector_ID_t DiagonalID );
 **********************************************************************/

/**********************************************************************
int Epetra_RowMatrix_Multiply ( 
  CT_Epetra_RowMatrix_ID_t selfID, boolean TransA, 
  CT_Epetra_MultiVector_ID_t XID, CT_Epetra_MultiVector_ID_t YID );
 **********************************************************************/

/**********************************************************************
int Epetra_RowMatrix_Solve ( 
  CT_Epetra_RowMatrix_ID_t selfID, boolean Upper, boolean Trans, 
  boolean UnitDiagonal, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t YID );
 **********************************************************************/

/**********************************************************************
int Epetra_RowMatrix_InvRowSums ( 
  CT_Epetra_RowMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID );
 **********************************************************************/

/**********************************************************************
int Epetra_RowMatrix_LeftScale ( 
  CT_Epetra_RowMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID );
 **********************************************************************/

/**********************************************************************
int Epetra_RowMatrix_InvColSums ( 
  CT_Epetra_RowMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID );
 **********************************************************************/

/**********************************************************************
int Epetra_RowMatrix_RightScale ( 
  CT_Epetra_RowMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID );
 **********************************************************************/

/**********************************************************************
boolean Epetra_RowMatrix_Filled ( CT_Epetra_RowMatrix_ID_t selfID );
 **********************************************************************/

/**********************************************************************
double Epetra_RowMatrix_NormInf ( CT_Epetra_RowMatrix_ID_t selfID );
 **********************************************************************/

/**********************************************************************
double Epetra_RowMatrix_NormOne ( CT_Epetra_RowMatrix_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_RowMatrix_NumGlobalNonzeros ( 
  CT_Epetra_RowMatrix_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_RowMatrix_NumGlobalRows ( 
  CT_Epetra_RowMatrix_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_RowMatrix_NumGlobalCols ( 
  CT_Epetra_RowMatrix_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_RowMatrix_NumGlobalDiagonals ( 
  CT_Epetra_RowMatrix_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_RowMatrix_NumMyNonzeros ( 
  CT_Epetra_RowMatrix_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_RowMatrix_NumMyRows ( CT_Epetra_RowMatrix_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_RowMatrix_NumMyCols ( CT_Epetra_RowMatrix_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_RowMatrix_NumMyDiagonals ( 
  CT_Epetra_RowMatrix_ID_t selfID );
 **********************************************************************/

/**********************************************************************
boolean Epetra_RowMatrix_LowerTriangular ( 
  CT_Epetra_RowMatrix_ID_t selfID );
 **********************************************************************/

/**********************************************************************
boolean Epetra_RowMatrix_UpperTriangular ( 
  CT_Epetra_RowMatrix_ID_t selfID );
 **********************************************************************/

/**********************************************************************
CT_Epetra_Map_ID_t Epetra_RowMatrix_RowMatrixRowMap ( 
  CT_Epetra_RowMatrix_ID_t selfID );
 **********************************************************************/

/**********************************************************************
CT_Epetra_Map_ID_t Epetra_RowMatrix_RowMatrixColMap ( 
  CT_Epetra_RowMatrix_ID_t selfID );
 **********************************************************************/

/**********************************************************************
CT_Epetra_Import_ID_t Epetra_RowMatrix_RowMatrixImporter ( 
  CT_Epetra_RowMatrix_ID_t selfID );
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

