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
#include "CTrilinos_enums.h"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_test_utils.hpp"
#include "CTrilinos_Table.hpp"
#include "CTrilinos_table_man.h"
#include "CEpetra_SerialComm.h"
#include "Epetra_SerialComm.h"
#include "CEpetra_Comm.h"
#include "Epetra_Comm.h"

#include "Teuchos_RCP.hpp"

#include "Teuchos_UnitTestHarness.hpp"


#define JOIN_SET_0(A, B, C) A ## B ## C
#define JOIN_SET(A, B, C)   JOIN_SET_0(A, B, C)

#define BUILD_CALL(A, F) JOIN_SET( A , _ , F )
#define CLASS_TYPE(A)    JOIN_SET( CT_ , A , _ID_t )
#define CLASS_ENUM(A)    JOIN_SET( CT_ , A , _ID )
#define CLASS_ESTR(A)    XSTRFY(CLASS_ENUM(A))
#define STRFY(A)         #A
#define XSTRFY(A)        STRFY(A)
#define CONSTRUCTOR(A)   A


#define T1 Epetra_SerialComm
#define T2 Epetra_Comm


namespace {


using Teuchos::null;
using Teuchos::RCP;
using Teuchos::rcp;
using CTrilinos::Table;

TEUCHOS_UNIT_TEST( Utils, isSameObjectRRTrue )
{
  ECHO(Table<T1> table1(CONSTRUCTOR(CLASS_ENUM(T1))));
  ECHO(Table<T2> table2(CONSTRUCTOR(CLASS_ENUM(T2))));

  ECHO(CTrilinos_Universal_ID_t id1 = table1.store<T1>(new T1, true));
  ECHO(CTrilinos_Universal_ID_t id2 = table2.alias(table1.get<T1>(id1)));

  TEST_EQUALITY_CONST(CTrilinos::isSameObject(table1.get<T1>(id1), table2.get<T1>(id2)), true);
  TEST_EQUALITY_CONST(CTrilinos::isSameObject(table2.get<T1>(id2), table1.get<T1>(id1)), true);
}

TEUCHOS_UNIT_TEST( Utils, isSameObjectRRFalse )
{
  ECHO(Table<T1> table1(CONSTRUCTOR(CLASS_ENUM(T1))));

  ECHO(CTrilinos_Universal_ID_t id1 = table1.store<T1>(new T1, true));
  ECHO(CTrilinos_Universal_ID_t id2 = table1.store<T1>(new T1, true));

  TEST_EQUALITY_CONST(CTrilinos::isSameObject(table1.get<T1>(id1), table1.get<T1>(id2)), false);
  TEST_EQUALITY_CONST(CTrilinos::isSameObject(table1.get<T1>(id2), table1.get<T1>(id1)), false);
}

TEUCHOS_UNIT_TEST( Utils, isSameObjectRATrue )
{
  ECHO(CT_Epetra_SerialComm_ID_t id1 = CEpetra::storeSerialComm(new Epetra_SerialComm));
  ECHO(CTrilinos_Universal_ID_t aid1 = CTrilinos::abstractType(id1));

  ECHO(Teuchos::RCP<Epetra_SerialComm> rcp1 = CEpetra::getSerialComm(id1));

  TEST_EQUALITY_CONST(CTrilinos::isSameObject(rcp1, aid1), true);
  TEST_EQUALITY_CONST(CTrilinos::isSameObject(aid1, rcp1), true);

  ECHO(CTrilinos_Universal_ID_t aid2 = CT_Alias(aid1, CT_Epetra_Comm_ID));

  TEST_EQUALITY_CONST(CTrilinos::isSameObject(rcp1, aid2), true);
  TEST_EQUALITY_CONST(CTrilinos::isSameObject(aid2, rcp1), true);

  ECHO(Teuchos::RCP<Epetra_Comm> rcp2 = CEpetra::getComm(aid2));

  TEST_EQUALITY_CONST(CTrilinos::isSameObject(rcp2, aid1), true);
  TEST_EQUALITY_CONST(CTrilinos::isSameObject(aid1, rcp2), true);

  TEST_EQUALITY_CONST(CTrilinos::isSameObject(rcp2, aid2), true);
  TEST_EQUALITY_CONST(CTrilinos::isSameObject(aid2, rcp2), true);
}

TEUCHOS_UNIT_TEST( Utils, isSameObjectRITrue )
{
  ECHO(CT_Epetra_SerialComm_ID_t id1 = CEpetra::storeSerialComm(new Epetra_SerialComm));
  ECHO(CTrilinos_Universal_ID_t aid1 = CTrilinos::abstractType(id1));

  ECHO(Teuchos::RCP<Epetra_SerialComm> rcp1 = CEpetra::getSerialComm(id1));

  TEST_EQUALITY_CONST(CTrilinos::isSameObject(rcp1, id1), true);
  TEST_EQUALITY_CONST(CTrilinos::isSameObject(id1, rcp1), true);

  ECHO(CTrilinos_Universal_ID_t aid2 = CT_Alias(aid1, CT_Epetra_Comm_ID));
  ECHO(CT_Epetra_Comm_ID_t id2 = CTrilinos::concreteType<CT_Epetra_Comm_ID_t>(aid2));

  TEST_EQUALITY_CONST(CTrilinos::isSameObject(rcp1, id2), true);
  TEST_EQUALITY_CONST(CTrilinos::isSameObject(id2, rcp1), true);

  ECHO(Teuchos::RCP<Epetra_Comm> rcp2 = CEpetra::getComm(aid2));

  TEST_EQUALITY_CONST(CTrilinos::isSameObject(rcp2, id1), true);
  TEST_EQUALITY_CONST(CTrilinos::isSameObject(id1, rcp2), true);

  TEST_EQUALITY_CONST(CTrilinos::isSameObject(rcp2, id2), true);
  TEST_EQUALITY_CONST(CTrilinos::isSameObject(id2, rcp2), true);
}

TEUCHOS_UNIT_TEST( Utils, isSameObjectRIFalse )
{
  ECHO(CT_Epetra_SerialComm_ID_t id1 = CEpetra::storeSerialComm(new Epetra_SerialComm));
  ECHO(CT_Epetra_SerialComm_ID_t id2 = CEpetra::storeSerialComm(new Epetra_SerialComm));

  ECHO(Teuchos::RCP<Epetra_SerialComm> rcp1 = CEpetra::getSerialComm(id1));
  ECHO(Teuchos::RCP<Epetra_SerialComm> rcp2 = CEpetra::getSerialComm(id2));

  TEST_EQUALITY_CONST(CTrilinos::isSameObject(rcp1, id2), false);
  TEST_EQUALITY_CONST(CTrilinos::isSameObject(id2, rcp1), false);

  TEST_EQUALITY_CONST(CTrilinos::isSameObject(rcp2, id1), false);
  TEST_EQUALITY_CONST(CTrilinos::isSameObject(id1, rcp2), false);
}

TEUCHOS_UNIT_TEST( Utils, isSameObjectAATrue )
{
  ECHO(CT_Epetra_SerialComm_ID_t id1 = CEpetra::storeSerialComm(new Epetra_SerialComm));
  ECHO(CTrilinos_Universal_ID_t aid1 = CTrilinos::abstractType(id1));

  ECHO(CTrilinos_Universal_ID_t aid2 = CT_Alias(aid1, CT_Epetra_Comm_ID));

  TEST_EQUALITY_CONST(CTrilinos::isSameObject(aid1, aid2), true);
  TEST_EQUALITY_CONST(CTrilinos::isSameObject(aid2, aid1), true);
}

TEUCHOS_UNIT_TEST( Utils, isSameObjectAITrue )
{
  ECHO(CT_Epetra_SerialComm_ID_t id1 = CEpetra::storeSerialComm(new Epetra_SerialComm));
  ECHO(CTrilinos_Universal_ID_t aid1 = CTrilinos::abstractType(id1));

  ECHO(CTrilinos_Universal_ID_t aid2 = CT_Alias(aid1, CT_Epetra_Comm_ID));
  ECHO(CT_Epetra_Comm_ID_t id2 = CTrilinos::concreteType<CT_Epetra_Comm_ID_t>(aid2));

  TEST_EQUALITY_CONST(CTrilinos::isSameObject(aid1, id2), true);
  TEST_EQUALITY_CONST(CTrilinos::isSameObject(id1, aid2), true);

  TEST_EQUALITY_CONST(CTrilinos::isSameObject(aid2, id1), true);
  TEST_EQUALITY_CONST(CTrilinos::isSameObject(id2, aid1), true);
}

TEUCHOS_UNIT_TEST( Utils, isSameObjectIITrue )
{
  ECHO(CT_Epetra_SerialComm_ID_t id1 = CEpetra::storeSerialComm(new Epetra_SerialComm));
  ECHO(CTrilinos_Universal_ID_t aid1 = CTrilinos::abstractType(id1));

  ECHO(CTrilinos_Universal_ID_t aid2 = CT_Alias(aid1, CT_Epetra_Comm_ID));
  ECHO(CT_Epetra_Comm_ID_t id2 = CTrilinos::concreteType<CT_Epetra_Comm_ID_t>(aid2));

  TEST_EQUALITY_CONST(CTrilinos::isSameObject(id1, id2), true);
  TEST_EQUALITY_CONST(CTrilinos::isSameObject(id2, id1), true);
}

TEUCHOS_UNIT_TEST( Utils, isSameObjectIIFalse )
{
  ECHO(CT_Epetra_SerialComm_ID_t id1 = CEpetra::storeSerialComm(new Epetra_SerialComm));
  ECHO(CT_Epetra_SerialComm_ID_t id2 = CEpetra::storeSerialComm(new Epetra_SerialComm));

  TEST_EQUALITY_CONST(CTrilinos::isSameObject(id1, id2), false);
  TEST_EQUALITY_CONST(CTrilinos::isSameObject(id2, id1), false);
}


//
// Template Instantiations
//


#ifdef TEUCHOS_DEBUG

#  define DEBUG_UNIT_TEST_GROUP( TT ) \

#else

#  define DEBUG_UNIT_TEST_GROUP( TT )

#endif


#define UNIT_TEST_GROUP( TT ) \
  DEBUG_UNIT_TEST_GROUP( TT )


} // namespace
