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
#include "CTrilinos_TableRepos.hpp"
#include "Teuchos_RCP.hpp"

#include "Teuchos_UnitTestHarness.hpp"

#include "CEpetra_SerialComm.h"
#include "CEpetra_Comm.h"
#include "CEpetra_Vector.h"
#include "Epetra_SerialComm.h"
#include "Epetra_Comm.h"
#include "Epetra_Vector.h"


#define JOIN_SET_0(A, B, C) A ## B ## C
#define JOIN_SET(A, B, C)   JOIN_SET_0(A, B, C)

#define BUILD_CALL(A, F) JOIN_SET( A , _ , F )
#define CLASS_TYPE(A)    JOIN_SET( CT_ , A , _ID_t )
#define CLASS_ENUM(A)    JOIN_SET( CT_ , A , _ID )
#define CLASS_ESTR(A)    XSTRFY(CLASS_ENUM(A))
#define STRFY(A)         #A
#define XSTRFY(A)        STRFY(A)
#define CONSTRUCTOR(A)   A


#define T Epetra_SerialComm
#define T1 Epetra_SerialComm
#define T2 Epetra_Comm
#define T3 Epetra_SerialComm
#define T4 Epetra_Vector


namespace {


using Teuchos::null;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::RangeError;
using Teuchos::NullReferenceError;
using Teuchos::m_bad_cast;
using CTrilinos::CTrilinosTypeMismatchError;
using CTrilinos::CTrilinosConstCastError;
using CTrilinos::TableRepos;


/* Table::store() owned */

TEUCHOS_UNIT_TEST( TableRepos, store )
{
  ECHO(TableRepos repos);
  ECHO(CTrilinos_Universal_ID_t id = repos.store<T>(new T, true));
  TEST_EQUALITY_CONST(id.index, 0);
  TEST_EQUALITY_CONST(id.is_const, FALSE);
  TEST_EQUALITY(id.table, CLASS_ENUM(T));
  TEST_EQUALITY_CONST(nonnull(repos.get<T>(id)), true);
  TEST_EQUALITY_CONST(is_null(repos.get<T>(id)), false);
  TEST_EQUALITY_CONST(nonnull(repos.getConst<T>(id)), true);
  TEST_EQUALITY_CONST(is_null(repos.getConst<T>(id)), false);
}

TEUCHOS_UNIT_TEST( TableRepos, storeBase )
{
  ECHO(TableRepos repos);
  ECHO(CTrilinos_Universal_ID_t id = repos.store<T2>(new T1, true));
  TEST_EQUALITY_CONST(id.index, 0);
  TEST_EQUALITY_CONST(id.is_const, FALSE);
  TEST_EQUALITY(id.table, CLASS_ENUM(T2));
  TEST_EQUALITY_CONST(nonnull(repos.get<T2>(id)), true);
  TEST_EQUALITY_CONST(is_null(repos.get<T2>(id)), false);
  TEST_EQUALITY_CONST(nonnull(repos.getConst<T2>(id)), true);
  TEST_EQUALITY_CONST(is_null(repos.getConst<T2>(id)), false);
  TEST_EQUALITY_CONST(nonnull(repos.get<T1>(id)), true);
  TEST_EQUALITY_CONST(is_null(repos.get<T1>(id)), false);
  TEST_EQUALITY_CONST(nonnull(repos.getConst<T1>(id)), true);
  TEST_EQUALITY_CONST(is_null(repos.getConst<T1>(id)), false);
}

TEUCHOS_UNIT_TEST( TableRepos, storeNull )
{
  ECHO(TableRepos repos);
  ECHO(T* pobj = NULL);
  TEST_THROW(repos.store<T>(pobj, false), NullReferenceError); 
}


/* Table::store() non-owned */

TEUCHOS_UNIT_TEST( TableRepos, storeShared )
{
  ECHO(TableRepos repos);
  ECHO(T *pobj = new T);
  ECHO(CTrilinos_Universal_ID_t id = repos.store<T>(pobj, false));
  TEST_EQUALITY_CONST(id.index, 0);
  TEST_EQUALITY(id.table, CLASS_ENUM(T));
  TEST_EQUALITY(id.is_const, FALSE);
  TEST_EQUALITY_CONST(nonnull(repos.get<T>(id)), true);
  TEST_EQUALITY_CONST(is_null(repos.get<T>(id)), false);
  TEST_EQUALITY_CONST(nonnull(repos.getConst<T>(id)), true);
  TEST_EQUALITY_CONST(is_null(repos.getConst<T>(id)), false);
  ECHO(repos.remove(&id));
  ECHO(delete pobj);
}

TEUCHOS_UNIT_TEST( TableRepos, storeConstShared )
{
  ECHO(TableRepos repos);
  ECHO(const T *pobj = new T);
  ECHO(CTrilinos_Universal_ID_t id = repos.store<T>(pobj, false));
  TEST_EQUALITY_CONST(id.index, 0);
  TEST_EQUALITY(id.table, CLASS_ENUM(T));
  TEST_EQUALITY(id.is_const, TRUE);
  TEST_EQUALITY_CONST(nonnull(repos.getConst<T>(id)), true);
  TEST_EQUALITY_CONST(is_null(repos.getConst<T>(id)), false);
  TEST_THROW(nonnull(repos.get<T>(id)), CTrilinosConstCastError);
  TEST_THROW(is_null(repos.get<T>(id)), CTrilinosConstCastError);
  ECHO(repos.remove(&id));
  ECHO(delete pobj);
}

TEUCHOS_UNIT_TEST( TableRepos, storeSharedBase )
{
  ECHO(TableRepos repos);
  ECHO(T1 *pobj = new T1);
  ECHO(CTrilinos_Universal_ID_t id = repos.store<T2>(pobj, false));
  TEST_EQUALITY_CONST(id.index, 0);
  TEST_EQUALITY(id.table, CLASS_ENUM(T2));
  TEST_EQUALITY(id.is_const, FALSE);
  TEST_EQUALITY_CONST(nonnull(repos.get<T2>(id)), true);
  TEST_EQUALITY_CONST(is_null(repos.get<T2>(id)), false);
  TEST_EQUALITY_CONST(nonnull(repos.getConst<T2>(id)), true);
  TEST_EQUALITY_CONST(is_null(repos.getConst<T2>(id)), false);
  TEST_EQUALITY_CONST(nonnull(repos.get<T1>(id)), true);
  TEST_EQUALITY_CONST(is_null(repos.get<T1>(id)), false);
  TEST_EQUALITY_CONST(nonnull(repos.getConst<T1>(id)), true);
  TEST_EQUALITY_CONST(is_null(repos.getConst<T1>(id)), false);
  ECHO(repos.remove(&id));
  ECHO(delete pobj);
}

TEUCHOS_UNIT_TEST( TableRepos, storeConstSharedBase )
{
  ECHO(TableRepos repos);
  ECHO(const T1 *pobj = new T1);
  ECHO(CTrilinos_Universal_ID_t id = repos.store<T2>(pobj, false));
  TEST_EQUALITY_CONST(id.index, 0);
  TEST_EQUALITY(id.table, CLASS_ENUM(T2));
  TEST_EQUALITY(id.is_const, TRUE);
  TEST_EQUALITY_CONST(nonnull(repos.getConst<T2>(id)), true);
  TEST_EQUALITY_CONST(is_null(repos.getConst<T2>(id)), false);
  TEST_THROW(nonnull(repos.get<T2>(id)), CTrilinosConstCastError);
  TEST_THROW(is_null(repos.get<T2>(id)), CTrilinosConstCastError);
  TEST_EQUALITY_CONST(nonnull(repos.getConst<T1>(id)), true);
  TEST_EQUALITY_CONST(is_null(repos.getConst<T1>(id)), false);
  TEST_THROW(nonnull(repos.get<T1>(id)), CTrilinosConstCastError);
  TEST_THROW(is_null(repos.get<T1>(id)), CTrilinosConstCastError);
  ECHO(repos.remove(&id));
  ECHO(delete pobj);
}

TEUCHOS_UNIT_TEST( TableRepos, storeSharedNull )
{
  ECHO(TableRepos repos);
  ECHO(T* pobj = NULL);
  TEST_THROW(repos.store<T>(pobj, false), NullReferenceError); 
}

TEUCHOS_UNIT_TEST( TableRepos, storeConstSharedNull )
{
  ECHO(TableRepos repos);
  ECHO(const T* pobj = NULL);
  TEST_THROW(repos.store<T>(pobj, false), NullReferenceError); 
}


/* Table::remove() */

TEUCHOS_UNIT_TEST( TableRepos, remove )
{
  ECHO(TableRepos repos);
  ECHO(CTrilinos_Universal_ID_t id = repos.store<T>(new T, true));
  TEST_EQUALITY_CONST(id.index, 0);
  TEST_EQUALITY(id.table, CLASS_ENUM(T));
  TEST_EQUALITY_CONST(id.is_const, FALSE);
  ECHO(repos.remove(&id));
  TEST_EQUALITY_CONST(id.index, -1);
  TEST_EQUALITY(id.table, CLASS_ENUM(Invalid));
}

TEUCHOS_UNIT_TEST( TableRepos, removeConst )
{
  ECHO(TableRepos repos);
  ECHO(const T* pobj = new T);
  ECHO(CTrilinos_Universal_ID_t id = repos.store<T>(pobj, false));
  TEST_EQUALITY_CONST(id.index, 0);
  TEST_EQUALITY(id.table, CLASS_ENUM(T));
  TEST_EQUALITY_CONST(id.is_const, TRUE);
  ECHO(repos.remove(&id));
  TEST_EQUALITY_CONST(id.index, -1);
  TEST_EQUALITY(id.table, CLASS_ENUM(Invalid));
  ECHO(delete pobj);
}

#ifdef TEUCHOS_DEBUG

TEUCHOS_UNIT_TEST( TableRepos, removeInvalid )
{
  ECHO(TableRepos repos);
  ECHO(CTrilinos_Universal_ID_t id);
  ECHO(id.index = -1);
  ECHO(id.table = CLASS_ENUM(T));
  ECHO(id.is_const = FALSE);
  TEST_THROW(repos.remove(&id), RangeError);
}

#endif /* TEUCHOS_DEBUG */


/* Table::get() */

TEUCHOS_UNIT_TEST( TableRepos, get )
{
  ECHO(TableRepos repos);
  ECHO(CTrilinos_Universal_ID_t id = repos.store<T>(new T, true));
  ECHO(RCP<T> rcpT = repos.get<T>(id));
  TEST_EQUALITY_CONST(nonnull(rcpT), true);
  TEST_EQUALITY_CONST(is_null(rcpT), false);
  ECHO(RCP<const T> rcpCT = repos.getConst<T>(id));
  TEST_EQUALITY_CONST(nonnull(rcpCT), true);
  TEST_EQUALITY_CONST(is_null(rcpCT), false);
}

#ifdef TEUCHOS_DEBUG

TEUCHOS_UNIT_TEST( TableRepos, getInvalid )
{
  ECHO(TableRepos repos);
  ECHO(CTrilinos_Universal_ID_t id);
  ECHO(id.index = 0);
  ECHO(id.table = CLASS_ENUM(T));
  ECHO(id.is_const = FALSE);
  TEST_THROW(repos.get<T>(id), RangeError);
}

#endif /* TEUCHOS_DEBUG */


/* Table::alias() */

TEUCHOS_UNIT_TEST( TableRepos, alias )
{
  ECHO(TableRepos repos);
  ECHO(CTrilinos_Universal_ID_t id1 = repos.store<T1>(new T1, true));
  ECHO(CTrilinos_Universal_ID_t id2 = repos.alias(id1, CLASS_ENUM(T2), true));
  TEST_EQUALITY(id2.table, CLASS_ENUM(T2));
  TEST_EQUALITY_CONST(id2.index, 0);
  TEST_EQUALITY_CONST(id2.is_const, FALSE);
}

TEUCHOS_UNIT_TEST( TableRepos, aliasConst )
{
  ECHO(TableRepos repos);
  ECHO(CTrilinos_Universal_ID_t id1 = repos.store<T1>(new const T1, true));
  ECHO(CTrilinos_Universal_ID_t id2 = repos.alias(id1, CLASS_ENUM(T2), true));
  TEST_EQUALITY(id2.table, CLASS_ENUM(T2));
  TEST_EQUALITY_CONST(id2.index, 0);
  TEST_EQUALITY_CONST(id2.is_const, TRUE);
}

TEUCHOS_UNIT_TEST( TableRepos, aliasBad )
{
  ECHO(TableRepos repos);
  ECHO(CTrilinos_Universal_ID_t id3 = repos.store<T3>(new T3, true));
  TEST_THROW(repos.alias(id3, CLASS_ENUM(T4), true), m_bad_cast);
}


/* Table::purge() */

#ifdef TEUCHOS_DEBUG

TEUCHOS_UNIT_TEST( TableRepos, purge )
{
  ECHO(TableRepos repos);
  ECHO(CTrilinos_Universal_ID_t id = repos.store<T>(new T, true));
  TEST_EQUALITY_CONST(nonnull(repos.get<T>(id)), true);
  ECHO(repos.purge<T>());
  TEST_THROW(repos.get<T>(id), RangeError);
}

TEUCHOS_UNIT_TEST( TableRepos, purgeAll )
{
  ECHO(TableRepos repos);
  ECHO(CTrilinos_Universal_ID_t id = repos.store<T>(new T, true));
  TEST_EQUALITY_CONST(nonnull(repos.get<T>(id)), true);
  ECHO(repos.purgeAll());
  TEST_THROW(repos.get<T>(id), RangeError);
}

#endif /* TEUCHOS_DEBUG */


/* Table::typeCheck() */

TEUCHOS_UNIT_TEST( TableRepos, typeCheck )
{
  ECHO(TableRepos repos);
  ECHO(CTrilinos_Universal_ID_t id = repos.store<T1>(new T1, true));
  TEST_EQUALITY_CONST(repos.typeCheck(id, CLASS_ENUM(T1)), true);
  TEST_EQUALITY_CONST(repos.typeCheck(id, CLASS_ENUM(T2)), true);
  TEST_EQUALITY_CONST(repos.typeCheck(id, CLASS_ENUM(T4)), false);
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
