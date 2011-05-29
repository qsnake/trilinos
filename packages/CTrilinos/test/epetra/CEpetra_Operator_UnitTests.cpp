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
#include "CEpetra_Map.h"
#include "CEpetra_MultiVector.h"
#include "CEpetra_CrsMatrix.h"
#include "Epetra_Operator.h"
#include "CEpetra_Operator.h"
#include "CEpetra_Operator_Cpp.hpp"
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
void Epetra_Operator_Destroy ( CT_Epetra_Operator_ID_t * selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_Operator , Destroy )
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

  ECHO(Epetra_Operator_Destroy(&crsID.Epetra_Operator));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(crsID.Epetra_Operator.table, CT_Invalid_ID);
  TEST_EQUALITY_CONST(crsID.Epetra_Operator.index, -1);
}

/**********************************************************************
int Epetra_Operator_SetUseTranspose ( 
  CT_Epetra_Operator_ID_t selfID, boolean UseTranspose );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_Operator , SetUseTranspose )
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

  /* Test true */
  ECHO(boolean tr = TRUE);
  ECHO(ret = Epetra_Operator_SetUseTranspose(crsID.Epetra_Operator, tr));
  TEST_EQUALITY_CONST(ret, 0);
  ECHO(boolean tr2 = Epetra_Operator_UseTranspose(crsID.Epetra_Operator));
  TEST_EQUALITY(tr, tr2);

  /* Test false */
  ECHO(tr = FALSE);
  ECHO(ret = Epetra_Operator_SetUseTranspose(crsID.Epetra_Operator, tr));
  TEST_EQUALITY_CONST(ret, 0);
  ECHO(tr2 = Epetra_Operator_UseTranspose(crsID.Epetra_Operator));
  TEST_EQUALITY(tr, tr2);
}

/**********************************************************************
int Epetra_Operator_Apply ( 
  CT_Epetra_Operator_ID_t selfID, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t YID );
 **********************************************************************/

/**********************************************************************
int Epetra_Operator_ApplyInverse ( 
  CT_Epetra_Operator_ID_t selfID, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t YID );
 **********************************************************************/

/**********************************************************************
double Epetra_Operator_NormInf ( CT_Epetra_Operator_ID_t selfID );
 **********************************************************************/

/**********************************************************************
const char * Epetra_Operator_Label ( 
  CT_Epetra_Operator_ID_t selfID );
 **********************************************************************/

/**********************************************************************
boolean Epetra_Operator_UseTranspose ( 
  CT_Epetra_Operator_ID_t selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_Operator , UseTranspose )
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

  /* Test true */
  ECHO(boolean tr = TRUE);
  ECHO(ret = Epetra_Operator_SetUseTranspose(crsID.Epetra_Operator, tr));
  TEST_EQUALITY_CONST(ret, 0);
  ECHO(boolean tr2 = Epetra_Operator_UseTranspose(crsID.Epetra_Operator));
  TEST_EQUALITY(tr, tr2);

  /* Test false */
  ECHO(tr = FALSE);
  ECHO(ret = Epetra_Operator_SetUseTranspose(crsID.Epetra_Operator, tr));
  TEST_EQUALITY_CONST(ret, 0);
  ECHO(tr2 = Epetra_Operator_UseTranspose(crsID.Epetra_Operator));
  TEST_EQUALITY(tr, tr2);
}

/**********************************************************************
boolean Epetra_Operator_HasNormInf ( 
  CT_Epetra_Operator_ID_t selfID );
 **********************************************************************/

/**********************************************************************
CT_Epetra_Comm_ID_t Epetra_Operator_Comm ( 
  CT_Epetra_Operator_ID_t selfID );
 **********************************************************************/

/**********************************************************************
CT_Epetra_Map_ID_t Epetra_Operator_OperatorDomainMap ( 
  CT_Epetra_Operator_ID_t selfID );
 **********************************************************************/

/**********************************************************************
CT_Epetra_Map_ID_t Epetra_Operator_OperatorRangeMap ( 
  CT_Epetra_Operator_ID_t selfID );
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

