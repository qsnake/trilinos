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
#include "Epetra_Distributor.h"
#include "CEpetra_Distributor.h"
#include "CEpetra_Distributor_Cpp.hpp"
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
CT_Epetra_Distributor_ID_t Epetra_Distributor_Clone ( 
  CT_Epetra_Distributor_ID_t selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_Distributor , Clone )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(CT_Epetra_Distributor_ID_t selfID = Epetra_Comm_CreateDistributor(CommID));

  ECHO(CT_Epetra_Distributor_ID_t dupID = Epetra_Distributor_Clone(selfID));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(dupID.table, CT_Epetra_Distributor_ID);
  TEST_EQUALITY_CONST(dupID.index, 1);
}

/**********************************************************************
void Epetra_Distributor_Destroy ( 
  CT_Epetra_Distributor_ID_t * selfID );
 **********************************************************************/

TEUCHOS_UNIT_TEST( Epetra_Distributor , Destroy )
{
  ECHO(CEpetra_Test_CleanSlate());

  /* Create everything we need to pass to the constructor */
  ECHO(CT_Epetra_Comm_ID_t CommID = UnitTest_Create_Comm());
  ECHO(CT_Epetra_Distributor_ID_t selfID = Epetra_Comm_CreateDistributor(CommID));

  ECHO(Epetra_Distributor_Destroy(&selfID));

  /* Now check the result of the call to the wrapper function */
  TEST_EQUALITY(selfID.table, CT_Invalid_ID);
  TEST_EQUALITY_CONST(selfID.index, -1);
}

/**********************************************************************
int Epetra_Distributor_CreateFromSends ( 
  CT_Epetra_Distributor_ID_t selfID, int NumExportIDs, 
  const int * ExportPIDs, boolean Deterministic, 
  int * NumRemoteIDs );
 **********************************************************************/

/**********************************************************************
int Epetra_Distributor_CreateFromRecvs ( 
  CT_Epetra_Distributor_ID_t selfID, int NumRemoteIDs, 
  const int * RemoteGIDs, const int * RemotePIDs, 
  boolean Deterministic, int * NumExportIDs, int ** ExportGIDs, 
  int ** ExportPIDs );
 **********************************************************************/

/**********************************************************************
int Epetra_Distributor_Do ( 
  CT_Epetra_Distributor_ID_t selfID, char * export_objs, int obj_size, 
  int * len_import_objs, char ** import_objs );
 **********************************************************************/

/**********************************************************************
int Epetra_Distributor_DoReverse ( 
  CT_Epetra_Distributor_ID_t selfID, char * export_objs, int obj_size, 
  int * len_import_objs, char ** import_objs );
 **********************************************************************/

/**********************************************************************
int Epetra_Distributor_DoPosts ( 
  CT_Epetra_Distributor_ID_t selfID, char * export_objs, int obj_size, 
  int * len_import_objs, char ** import_objs );
 **********************************************************************/

/**********************************************************************
int Epetra_Distributor_DoWaits ( CT_Epetra_Distributor_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_Distributor_DoReversePosts ( 
  CT_Epetra_Distributor_ID_t selfID, char * export_objs, int obj_size, 
  int * len_import_objs, char ** import_objs );
 **********************************************************************/

/**********************************************************************
int Epetra_Distributor_DoReverseWaits ( 
  CT_Epetra_Distributor_ID_t selfID );
 **********************************************************************/

/**********************************************************************
int Epetra_Distributor_Do_VarLen ( 
  CT_Epetra_Distributor_ID_t selfID, char * export_objs, int obj_size, 
  int ** sizes, int * len_import_objs, char ** import_objs );
 **********************************************************************/

/**********************************************************************
int Epetra_Distributor_DoReverse_VarLen ( 
  CT_Epetra_Distributor_ID_t selfID, char * export_objs, int obj_size, 
  int ** sizes, int * len_import_objs, char ** import_objs );
 **********************************************************************/

/**********************************************************************
int Epetra_Distributor_DoPosts_VarLen ( 
  CT_Epetra_Distributor_ID_t selfID, char * export_objs, int obj_size, 
  int ** sizes, int * len_import_objs, char ** import_objs );
 **********************************************************************/

/**********************************************************************
int Epetra_Distributor_DoReversePosts_VarLen ( 
  CT_Epetra_Distributor_ID_t selfID, char * export_objs, int obj_size, 
  int ** sizes, int * len_import_objs, char ** import_objs );
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

