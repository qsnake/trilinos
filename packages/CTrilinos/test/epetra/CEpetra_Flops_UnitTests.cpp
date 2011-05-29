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
#include "CEpetra_CompObject.h"
#include "Epetra_Flops.h"
#include "CEpetra_Flops.h"
#include "CEpetra_Flops_Cpp.hpp"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_exceptions.hpp"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_test_utils.hpp"

#include "CTrilinos_UnitTestHelpers.hpp"
#include "Teuchos_UnitTestHarness.hpp"


namespace {


/**********************************************************************
CT_Epetra_Flops_ID_t Epetra_Flops_Create (  );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_Flops , Create )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(CT_Epetra_Flops_ID_t selfID = Epetra_Flops_Create());

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(selfID.table, CT_Epetra_Flops_ID);
  TEST_EQUALITY_CONST(selfID.index, 0);
}

/**********************************************************************
CT_Epetra_Flops_ID_t Epetra_Flops_Duplicate ( 
  CT_Epetra_Flops_ID_t Flops_inID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_Flops , Duplicate )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(CT_Epetra_Flops_ID_t selfID = Epetra_Flops_Create());
  ECHO(CT_Epetra_Flops_ID_t dupID = Epetra_Flops_Duplicate(selfID));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(dupID.table, CT_Epetra_Flops_ID);
  TEST_EQUALITY_CONST(dupID.index, 1);
  TEST_EQUALITY_CONST(CTrilinos::isSameObject(selfID, dupID), false);
}

/**********************************************************************
double Epetra_Flops_Flops ( CT_Epetra_Flops_ID_t selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_Flops , Flops )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(CT_Epetra_CompObject_ID_t compID = Epetra_CompObject_Create());
  ECHO(CT_Epetra_Flops_ID_t fID = Epetra_Flops_Create());
  ECHO(Epetra_CompObject_SetFlopCounter(compID, fID));
  ECHO(double f = 23456.125);
  ECHO(Epetra_CompObject_UpdateFlops_Double(compID, f));
  ECHO(CT_Epetra_Flops_ID_t fID2 = Epetra_CompObject_GetFlopCounter(compID));

  ECHO(double f2 = Epetra_Flops_Flops(fID2));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(f2, f);
}

/**********************************************************************
void Epetra_Flops_ResetFlops ( CT_Epetra_Flops_ID_t selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_Flops , ResetFlops )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(CT_Epetra_CompObject_ID_t compID = Epetra_CompObject_Create());
  ECHO(CT_Epetra_Flops_ID_t fID = Epetra_Flops_Create());
  ECHO(Epetra_CompObject_SetFlopCounter(compID, fID));
  ECHO(double f = 23456.125);
  ECHO(Epetra_CompObject_UpdateFlops_Double(compID, f));
  ECHO(CT_Epetra_Flops_ID_t fID2 = Epetra_CompObject_GetFlopCounter(compID));
  ECHO(double f2 = Epetra_Flops_Flops(fID2));
  TEST_EQUALITY(f2, f);

  ECHO(Epetra_Flops_ResetFlops(fID2));
  ECHO(double f3 = Epetra_Flops_Flops(fID2));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY_CONST(f3, 0.0);
}

/**********************************************************************
void Epetra_Flops_Destroy ( CT_Epetra_Flops_ID_t * selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_Flops , Destroy )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(CT_Epetra_Flops_ID_t selfID = Epetra_Flops_Create());

  ECHO(Epetra_Flops_Destroy(&selfID));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(selfID.table, CT_Invalid_ID);
  TEST_EQUALITY_CONST(selfID.index, -1);
}

/**********************************************************************
void Epetra_Flops_Assign ( 
  CT_Epetra_Flops_ID_t selfID, CT_Epetra_Flops_ID_t srcID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_Flops , Assign )
{
  ECHO(CEpetra_Test_CleanSlate());

  ECHO(CT_Epetra_CompObject_ID_t compID = Epetra_CompObject_Create());
  ECHO(CT_Epetra_Flops_ID_t fID = Epetra_Flops_Create());
  ECHO(Epetra_CompObject_SetFlopCounter(compID, fID));
  ECHO(double f = 8468.432);
  ECHO(Epetra_CompObject_UpdateFlops_Double(compID, f));
  ECHO(CT_Epetra_Flops_ID_t fID2 = Epetra_CompObject_GetFlopCounter(compID));

  ECHO(CT_Epetra_Flops_ID_t selfID = Epetra_Flops_Create());
  ECHO(Epetra_Flops_Assign(selfID, fID2));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY_CONST(CTrilinos::isSameObject(selfID, fID2), false);
  ECHO(double f2 = Epetra_Flops_Flops(selfID));
  TEST_EQUALITY(f2, f);
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

