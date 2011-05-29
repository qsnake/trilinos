
/*! @HEADER */
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
Questions? Contact M. Nicole Lemaster (mnlemas@sandia.gov)

************************************************************************
*/
/*! @HEADER */


/*! @file CTrilinos_table_calls.cpp
 * @brief Calls to pull class instances out of the tables. */


#include "CTrilinos_TableRepos.hpp"
#include "CTrilinos_utils_templ.hpp"

#include "CEpetra_Distributor_Cpp.hpp"
#include "CEpetra_SerialComm_Cpp.hpp"
#include "CEpetra_BLAS_Cpp.hpp"
#include "CEpetra_Comm_Cpp.hpp"
#include "CEpetra_Operator_Cpp.hpp"
#include "CEpetra_MultiVector_Cpp.hpp"
#include "CEpetra_OffsetIndex_Cpp.hpp"
#include "CEpetra_Object_Cpp.hpp"
#include "CEpetra_RowMatrix_Cpp.hpp"
#include "CEpetra_CompObject_Cpp.hpp"
#include "CEpetra_Directory_Cpp.hpp"
#include "CEpetra_Flops_Cpp.hpp"
#include "CEpetra_SrcDistObject_Cpp.hpp"
#ifdef HAVE_MPI
#include "CEpetra_MpiComm_Cpp.hpp"
#endif /* HAVE_MPI */
#include "CEpetra_CrsMatrix_Cpp.hpp"
#include "CEpetra_CrsGraph_Cpp.hpp"
#include "CEpetra_DistObject_Cpp.hpp"
#include "CEpetra_Vector_Cpp.hpp"
#include "CEpetra_Export_Cpp.hpp"
#include "CEpetra_Map_Cpp.hpp"
#include "CEpetra_BlockMap_Cpp.hpp"
#include "CEpetra_Import_Cpp.hpp"
#include "CEpetra_Time_Cpp.hpp"
#include "CEpetra_JadMatrix_Cpp.hpp"
#include "CEpetra_LinearProblem_Cpp.hpp"
#include "CEpetra_LAPACK_Cpp.hpp"
#include "CTeuchos_ParameterList_Cpp.hpp"
#ifdef HAVE_CTRILINOS_AMESOS
#include "CAmesos_BaseSolver_Cpp.hpp"
#endif /* HAVE_CTRILINOS_AMESOS */
#include "CEpetra_FECrsMatrix_Cpp.hpp"
#include "CEpetra_IntSerialDenseVector_Cpp.hpp"
#include "CEpetra_SerialDenseMatrix_Cpp.hpp"
#ifdef HAVE_CTRILINOS_AZTECOO
#include "CAztecOO_StatusTest_Cpp.hpp"
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
#include "CAztecOO_StatusTestCombo_Cpp.hpp"
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
#include "CAztecOO_StatusTestMaxIters_Cpp.hpp"
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
#include "CAztecOO_StatusTestResNorm_Cpp.hpp"
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_IFPACK
#include "CIfpack_Preconditioner_Cpp.hpp"
#endif /* HAVE_CTRILINOS_IFPACK */


//
// Definitions for Epetra_Distributor
//

/* get Epetra_Distributor from non-const table using CT_Epetra_Distributor_ID */
const Teuchos::RCP<Epetra_Distributor>
CEpetra::getDistributor( CT_Epetra_Distributor_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_Distributor>(
        CTrilinos::abstractType<CT_Epetra_Distributor_ID_t>(id));
}

/* get Epetra_Distributor from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_Distributor>
CEpetra::getDistributor( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_Distributor>(id);
}

/* get const Epetra_Distributor from either the const or non-const table
 * using CT_Epetra_Distributor_ID */
const Teuchos::RCP<const Epetra_Distributor>
CEpetra::getConstDistributor( CT_Epetra_Distributor_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_Distributor>(
        CTrilinos::abstractType<CT_Epetra_Distributor_ID_t>(id));
}

/* get const Epetra_Distributor from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_Distributor>
CEpetra::getConstDistributor( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_Distributor>(id);
}

/* store Epetra_Distributor (owned) in non-const table */
CT_Epetra_Distributor_ID_t
CEpetra::storeNewDistributor( Epetra_Distributor *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Distributor_ID_t>(
        CTrilinos::tableRepos().store<Epetra_Distributor>(pobj, true));
}

/* store Epetra_Distributor in non-const table */
CT_Epetra_Distributor_ID_t
CEpetra::storeDistributor( Epetra_Distributor *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Distributor_ID_t>(
        CTrilinos::tableRepos().store<Epetra_Distributor>(pobj, false));
}

/* store const Epetra_Distributor in const table */
CT_Epetra_Distributor_ID_t
CEpetra::storeConstDistributor( const Epetra_Distributor *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Distributor_ID_t>(
        CTrilinos::tableRepos().store<Epetra_Distributor>(pobj, false));
}

/* remove Epetra_Distributor from table using CT_Epetra_Distributor_ID */
void
CEpetra::removeDistributor( CT_Epetra_Distributor_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_Distributor_ID_t>(*id);
    CTrilinos::tableRepos().remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_Distributor_ID_t>(aid);
}

/* purge Epetra_Distributor table */
void
CEpetra::purgeDistributor(  )
{
    CTrilinos::tableRepos().purge<Epetra_Distributor>();
}


//
// Definitions for Epetra_SerialComm
//

/* get Epetra_SerialComm from non-const table using CT_Epetra_SerialComm_ID */
const Teuchos::RCP<Epetra_SerialComm>
CEpetra::getSerialComm( CT_Epetra_SerialComm_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_SerialComm>(
        CTrilinos::abstractType<CT_Epetra_SerialComm_ID_t>(id));
}

/* get Epetra_SerialComm from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_SerialComm>
CEpetra::getSerialComm( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_SerialComm>(id);
}

/* get const Epetra_SerialComm from either the const or non-const table
 * using CT_Epetra_SerialComm_ID */
const Teuchos::RCP<const Epetra_SerialComm>
CEpetra::getConstSerialComm( CT_Epetra_SerialComm_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_SerialComm>(
        CTrilinos::abstractType<CT_Epetra_SerialComm_ID_t>(id));
}

/* get const Epetra_SerialComm from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_SerialComm>
CEpetra::getConstSerialComm( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_SerialComm>(id);
}

/* store Epetra_SerialComm (owned) in non-const table */
CT_Epetra_SerialComm_ID_t
CEpetra::storeNewSerialComm( Epetra_SerialComm *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_SerialComm_ID_t>(
        CTrilinos::tableRepos().store<Epetra_SerialComm>(pobj, true));
}

/* store Epetra_SerialComm in non-const table */
CT_Epetra_SerialComm_ID_t
CEpetra::storeSerialComm( Epetra_SerialComm *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_SerialComm_ID_t>(
        CTrilinos::tableRepos().store<Epetra_SerialComm>(pobj, false));
}

/* store const Epetra_SerialComm in const table */
CT_Epetra_SerialComm_ID_t
CEpetra::storeConstSerialComm( const Epetra_SerialComm *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_SerialComm_ID_t>(
        CTrilinos::tableRepos().store<Epetra_SerialComm>(pobj, false));
}

/* remove Epetra_SerialComm from table using CT_Epetra_SerialComm_ID */
void
CEpetra::removeSerialComm( CT_Epetra_SerialComm_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_SerialComm_ID_t>(*id);
    CTrilinos::tableRepos().remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_SerialComm_ID_t>(aid);
}

/* purge Epetra_SerialComm table */
void
CEpetra::purgeSerialComm(  )
{
    CTrilinos::tableRepos().purge<Epetra_SerialComm>();
}


//
// Definitions for Epetra_BLAS
//

/* get Epetra_BLAS from non-const table using CT_Epetra_BLAS_ID */
const Teuchos::RCP<Epetra_BLAS>
CEpetra::getBLAS( CT_Epetra_BLAS_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_BLAS>(
        CTrilinos::abstractType<CT_Epetra_BLAS_ID_t>(id));
}

/* get Epetra_BLAS from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_BLAS>
CEpetra::getBLAS( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_BLAS>(id);
}

/* get const Epetra_BLAS from either the const or non-const table
 * using CT_Epetra_BLAS_ID */
const Teuchos::RCP<const Epetra_BLAS>
CEpetra::getConstBLAS( CT_Epetra_BLAS_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_BLAS>(
        CTrilinos::abstractType<CT_Epetra_BLAS_ID_t>(id));
}

/* get const Epetra_BLAS from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_BLAS>
CEpetra::getConstBLAS( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_BLAS>(id);
}

/* store Epetra_BLAS (owned) in non-const table */
CT_Epetra_BLAS_ID_t
CEpetra::storeNewBLAS( Epetra_BLAS *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_BLAS_ID_t>(
        CTrilinos::tableRepos().store<Epetra_BLAS>(pobj, true));
}

/* store Epetra_BLAS in non-const table */
CT_Epetra_BLAS_ID_t
CEpetra::storeBLAS( Epetra_BLAS *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_BLAS_ID_t>(
        CTrilinos::tableRepos().store<Epetra_BLAS>(pobj, false));
}

/* store const Epetra_BLAS in const table */
CT_Epetra_BLAS_ID_t
CEpetra::storeConstBLAS( const Epetra_BLAS *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_BLAS_ID_t>(
        CTrilinos::tableRepos().store<Epetra_BLAS>(pobj, false));
}

/* remove Epetra_BLAS from table using CT_Epetra_BLAS_ID */
void
CEpetra::removeBLAS( CT_Epetra_BLAS_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_BLAS_ID_t>(*id);
    CTrilinos::tableRepos().remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_BLAS_ID_t>(aid);
}

/* purge Epetra_BLAS table */
void
CEpetra::purgeBLAS(  )
{
    CTrilinos::tableRepos().purge<Epetra_BLAS>();
}


//
// Definitions for Epetra_Comm
//

/* get Epetra_Comm from non-const table using CT_Epetra_Comm_ID */
const Teuchos::RCP<Epetra_Comm>
CEpetra::getComm( CT_Epetra_Comm_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_Comm>(
        CTrilinos::abstractType<CT_Epetra_Comm_ID_t>(id));
}

/* get Epetra_Comm from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_Comm>
CEpetra::getComm( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_Comm>(id);
}

/* get const Epetra_Comm from either the const or non-const table
 * using CT_Epetra_Comm_ID */
const Teuchos::RCP<const Epetra_Comm>
CEpetra::getConstComm( CT_Epetra_Comm_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_Comm>(
        CTrilinos::abstractType<CT_Epetra_Comm_ID_t>(id));
}

/* get const Epetra_Comm from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_Comm>
CEpetra::getConstComm( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_Comm>(id);
}

/* store Epetra_Comm (owned) in non-const table */
CT_Epetra_Comm_ID_t
CEpetra::storeNewComm( Epetra_Comm *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Comm_ID_t>(
        CTrilinos::tableRepos().store<Epetra_Comm>(pobj, true));
}

/* store Epetra_Comm in non-const table */
CT_Epetra_Comm_ID_t
CEpetra::storeComm( Epetra_Comm *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Comm_ID_t>(
        CTrilinos::tableRepos().store<Epetra_Comm>(pobj, false));
}

/* store const Epetra_Comm in const table */
CT_Epetra_Comm_ID_t
CEpetra::storeConstComm( const Epetra_Comm *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Comm_ID_t>(
        CTrilinos::tableRepos().store<Epetra_Comm>(pobj, false));
}

/* remove Epetra_Comm from table using CT_Epetra_Comm_ID */
void
CEpetra::removeComm( CT_Epetra_Comm_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_Comm_ID_t>(*id);
    CTrilinos::tableRepos().remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_Comm_ID_t>(aid);
}

/* purge Epetra_Comm table */
void
CEpetra::purgeComm(  )
{
    CTrilinos::tableRepos().purge<Epetra_Comm>();
}


//
// Definitions for Epetra_Operator
//

/* get Epetra_Operator from non-const table using CT_Epetra_Operator_ID */
const Teuchos::RCP<Epetra_Operator>
CEpetra::getOperator( CT_Epetra_Operator_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_Operator>(
        CTrilinos::abstractType<CT_Epetra_Operator_ID_t>(id));
}

/* get Epetra_Operator from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_Operator>
CEpetra::getOperator( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_Operator>(id);
}

/* get const Epetra_Operator from either the const or non-const table
 * using CT_Epetra_Operator_ID */
const Teuchos::RCP<const Epetra_Operator>
CEpetra::getConstOperator( CT_Epetra_Operator_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_Operator>(
        CTrilinos::abstractType<CT_Epetra_Operator_ID_t>(id));
}

/* get const Epetra_Operator from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_Operator>
CEpetra::getConstOperator( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_Operator>(id);
}

/* store Epetra_Operator (owned) in non-const table */
CT_Epetra_Operator_ID_t
CEpetra::storeNewOperator( Epetra_Operator *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Operator_ID_t>(
        CTrilinos::tableRepos().store<Epetra_Operator>(pobj, true));
}

/* store Epetra_Operator in non-const table */
CT_Epetra_Operator_ID_t
CEpetra::storeOperator( Epetra_Operator *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Operator_ID_t>(
        CTrilinos::tableRepos().store<Epetra_Operator>(pobj, false));
}

/* store const Epetra_Operator in const table */
CT_Epetra_Operator_ID_t
CEpetra::storeConstOperator( const Epetra_Operator *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Operator_ID_t>(
        CTrilinos::tableRepos().store<Epetra_Operator>(pobj, false));
}

/* remove Epetra_Operator from table using CT_Epetra_Operator_ID */
void
CEpetra::removeOperator( CT_Epetra_Operator_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_Operator_ID_t>(*id);
    CTrilinos::tableRepos().remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_Operator_ID_t>(aid);
}

/* purge Epetra_Operator table */
void
CEpetra::purgeOperator(  )
{
    CTrilinos::tableRepos().purge<Epetra_Operator>();
}


//
// Definitions for Epetra_MultiVector
//

/* get Epetra_MultiVector from non-const table using CT_Epetra_MultiVector_ID */
const Teuchos::RCP<Epetra_MultiVector>
CEpetra::getMultiVector( CT_Epetra_MultiVector_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_MultiVector>(
        CTrilinos::abstractType<CT_Epetra_MultiVector_ID_t>(id));
}

/* get Epetra_MultiVector from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_MultiVector>
CEpetra::getMultiVector( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_MultiVector>(id);
}

/* get const Epetra_MultiVector from either the const or non-const table
 * using CT_Epetra_MultiVector_ID */
const Teuchos::RCP<const Epetra_MultiVector>
CEpetra::getConstMultiVector( CT_Epetra_MultiVector_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_MultiVector>(
        CTrilinos::abstractType<CT_Epetra_MultiVector_ID_t>(id));
}

/* get const Epetra_MultiVector from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_MultiVector>
CEpetra::getConstMultiVector( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_MultiVector>(id);
}

/* store Epetra_MultiVector (owned) in non-const table */
CT_Epetra_MultiVector_ID_t
CEpetra::storeNewMultiVector( Epetra_MultiVector *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_MultiVector_ID_t>(
        CTrilinos::tableRepos().store<Epetra_MultiVector>(pobj, true));
}

/* store Epetra_MultiVector in non-const table */
CT_Epetra_MultiVector_ID_t
CEpetra::storeMultiVector( Epetra_MultiVector *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_MultiVector_ID_t>(
        CTrilinos::tableRepos().store<Epetra_MultiVector>(pobj, false));
}

/* store const Epetra_MultiVector in const table */
CT_Epetra_MultiVector_ID_t
CEpetra::storeConstMultiVector( const Epetra_MultiVector *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_MultiVector_ID_t>(
        CTrilinos::tableRepos().store<Epetra_MultiVector>(pobj, false));
}

/* remove Epetra_MultiVector from table using CT_Epetra_MultiVector_ID */
void
CEpetra::removeMultiVector( CT_Epetra_MultiVector_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_MultiVector_ID_t>(*id);
    CTrilinos::tableRepos().remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_MultiVector_ID_t>(aid);
}

/* purge Epetra_MultiVector table */
void
CEpetra::purgeMultiVector(  )
{
    CTrilinos::tableRepos().purge<Epetra_MultiVector>();
}


//
// Definitions for Epetra_OffsetIndex
//

/* get Epetra_OffsetIndex from non-const table using CT_Epetra_OffsetIndex_ID */
const Teuchos::RCP<Epetra_OffsetIndex>
CEpetra::getOffsetIndex( CT_Epetra_OffsetIndex_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_OffsetIndex>(
        CTrilinos::abstractType<CT_Epetra_OffsetIndex_ID_t>(id));
}

/* get Epetra_OffsetIndex from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_OffsetIndex>
CEpetra::getOffsetIndex( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_OffsetIndex>(id);
}

/* get const Epetra_OffsetIndex from either the const or non-const table
 * using CT_Epetra_OffsetIndex_ID */
const Teuchos::RCP<const Epetra_OffsetIndex>
CEpetra::getConstOffsetIndex( CT_Epetra_OffsetIndex_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_OffsetIndex>(
        CTrilinos::abstractType<CT_Epetra_OffsetIndex_ID_t>(id));
}

/* get const Epetra_OffsetIndex from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_OffsetIndex>
CEpetra::getConstOffsetIndex( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_OffsetIndex>(id);
}

/* store Epetra_OffsetIndex (owned) in non-const table */
CT_Epetra_OffsetIndex_ID_t
CEpetra::storeNewOffsetIndex( Epetra_OffsetIndex *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_OffsetIndex_ID_t>(
        CTrilinos::tableRepos().store<Epetra_OffsetIndex>(pobj, true));
}

/* store Epetra_OffsetIndex in non-const table */
CT_Epetra_OffsetIndex_ID_t
CEpetra::storeOffsetIndex( Epetra_OffsetIndex *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_OffsetIndex_ID_t>(
        CTrilinos::tableRepos().store<Epetra_OffsetIndex>(pobj, false));
}

/* store const Epetra_OffsetIndex in const table */
CT_Epetra_OffsetIndex_ID_t
CEpetra::storeConstOffsetIndex( const Epetra_OffsetIndex *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_OffsetIndex_ID_t>(
        CTrilinos::tableRepos().store<Epetra_OffsetIndex>(pobj, false));
}

/* remove Epetra_OffsetIndex from table using CT_Epetra_OffsetIndex_ID */
void
CEpetra::removeOffsetIndex( CT_Epetra_OffsetIndex_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_OffsetIndex_ID_t>(*id);
    CTrilinos::tableRepos().remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_OffsetIndex_ID_t>(aid);
}

/* purge Epetra_OffsetIndex table */
void
CEpetra::purgeOffsetIndex(  )
{
    CTrilinos::tableRepos().purge<Epetra_OffsetIndex>();
}


//
// Definitions for Epetra_Object
//

/* get Epetra_Object from non-const table using CT_Epetra_Object_ID */
const Teuchos::RCP<Epetra_Object>
CEpetra::getObject( CT_Epetra_Object_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_Object>(
        CTrilinos::abstractType<CT_Epetra_Object_ID_t>(id));
}

/* get Epetra_Object from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_Object>
CEpetra::getObject( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_Object>(id);
}

/* get const Epetra_Object from either the const or non-const table
 * using CT_Epetra_Object_ID */
const Teuchos::RCP<const Epetra_Object>
CEpetra::getConstObject( CT_Epetra_Object_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_Object>(
        CTrilinos::abstractType<CT_Epetra_Object_ID_t>(id));
}

/* get const Epetra_Object from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_Object>
CEpetra::getConstObject( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_Object>(id);
}

/* store Epetra_Object (owned) in non-const table */
CT_Epetra_Object_ID_t
CEpetra::storeNewObject( Epetra_Object *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Object_ID_t>(
        CTrilinos::tableRepos().store<Epetra_Object>(pobj, true));
}

/* store Epetra_Object in non-const table */
CT_Epetra_Object_ID_t
CEpetra::storeObject( Epetra_Object *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Object_ID_t>(
        CTrilinos::tableRepos().store<Epetra_Object>(pobj, false));
}

/* store const Epetra_Object in const table */
CT_Epetra_Object_ID_t
CEpetra::storeConstObject( const Epetra_Object *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Object_ID_t>(
        CTrilinos::tableRepos().store<Epetra_Object>(pobj, false));
}

/* remove Epetra_Object from table using CT_Epetra_Object_ID */
void
CEpetra::removeObject( CT_Epetra_Object_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_Object_ID_t>(*id);
    CTrilinos::tableRepos().remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_Object_ID_t>(aid);
}

/* purge Epetra_Object table */
void
CEpetra::purgeObject(  )
{
    CTrilinos::tableRepos().purge<Epetra_Object>();
}


//
// Definitions for Epetra_RowMatrix
//

/* get Epetra_RowMatrix from non-const table using CT_Epetra_RowMatrix_ID */
const Teuchos::RCP<Epetra_RowMatrix>
CEpetra::getRowMatrix( CT_Epetra_RowMatrix_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_RowMatrix>(
        CTrilinos::abstractType<CT_Epetra_RowMatrix_ID_t>(id));
}

/* get Epetra_RowMatrix from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_RowMatrix>
CEpetra::getRowMatrix( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_RowMatrix>(id);
}

/* get const Epetra_RowMatrix from either the const or non-const table
 * using CT_Epetra_RowMatrix_ID */
const Teuchos::RCP<const Epetra_RowMatrix>
CEpetra::getConstRowMatrix( CT_Epetra_RowMatrix_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_RowMatrix>(
        CTrilinos::abstractType<CT_Epetra_RowMatrix_ID_t>(id));
}

/* get const Epetra_RowMatrix from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_RowMatrix>
CEpetra::getConstRowMatrix( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_RowMatrix>(id);
}

/* store Epetra_RowMatrix (owned) in non-const table */
CT_Epetra_RowMatrix_ID_t
CEpetra::storeNewRowMatrix( Epetra_RowMatrix *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_RowMatrix_ID_t>(
        CTrilinos::tableRepos().store<Epetra_RowMatrix>(pobj, true));
}

/* store Epetra_RowMatrix in non-const table */
CT_Epetra_RowMatrix_ID_t
CEpetra::storeRowMatrix( Epetra_RowMatrix *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_RowMatrix_ID_t>(
        CTrilinos::tableRepos().store<Epetra_RowMatrix>(pobj, false));
}

/* store const Epetra_RowMatrix in const table */
CT_Epetra_RowMatrix_ID_t
CEpetra::storeConstRowMatrix( const Epetra_RowMatrix *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_RowMatrix_ID_t>(
        CTrilinos::tableRepos().store<Epetra_RowMatrix>(pobj, false));
}

/* remove Epetra_RowMatrix from table using CT_Epetra_RowMatrix_ID */
void
CEpetra::removeRowMatrix( CT_Epetra_RowMatrix_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_RowMatrix_ID_t>(*id);
    CTrilinos::tableRepos().remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_RowMatrix_ID_t>(aid);
}

/* purge Epetra_RowMatrix table */
void
CEpetra::purgeRowMatrix(  )
{
    CTrilinos::tableRepos().purge<Epetra_RowMatrix>();
}


//
// Definitions for Epetra_CompObject
//

/* get Epetra_CompObject from non-const table using CT_Epetra_CompObject_ID */
const Teuchos::RCP<Epetra_CompObject>
CEpetra::getCompObject( CT_Epetra_CompObject_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_CompObject>(
        CTrilinos::abstractType<CT_Epetra_CompObject_ID_t>(id));
}

/* get Epetra_CompObject from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_CompObject>
CEpetra::getCompObject( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_CompObject>(id);
}

/* get const Epetra_CompObject from either the const or non-const table
 * using CT_Epetra_CompObject_ID */
const Teuchos::RCP<const Epetra_CompObject>
CEpetra::getConstCompObject( CT_Epetra_CompObject_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_CompObject>(
        CTrilinos::abstractType<CT_Epetra_CompObject_ID_t>(id));
}

/* get const Epetra_CompObject from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_CompObject>
CEpetra::getConstCompObject( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_CompObject>(id);
}

/* store Epetra_CompObject (owned) in non-const table */
CT_Epetra_CompObject_ID_t
CEpetra::storeNewCompObject( Epetra_CompObject *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_CompObject_ID_t>(
        CTrilinos::tableRepos().store<Epetra_CompObject>(pobj, true));
}

/* store Epetra_CompObject in non-const table */
CT_Epetra_CompObject_ID_t
CEpetra::storeCompObject( Epetra_CompObject *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_CompObject_ID_t>(
        CTrilinos::tableRepos().store<Epetra_CompObject>(pobj, false));
}

/* store const Epetra_CompObject in const table */
CT_Epetra_CompObject_ID_t
CEpetra::storeConstCompObject( const Epetra_CompObject *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_CompObject_ID_t>(
        CTrilinos::tableRepos().store<Epetra_CompObject>(pobj, false));
}

/* remove Epetra_CompObject from table using CT_Epetra_CompObject_ID */
void
CEpetra::removeCompObject( CT_Epetra_CompObject_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_CompObject_ID_t>(*id);
    CTrilinos::tableRepos().remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_CompObject_ID_t>(aid);
}

/* purge Epetra_CompObject table */
void
CEpetra::purgeCompObject(  )
{
    CTrilinos::tableRepos().purge<Epetra_CompObject>();
}


//
// Definitions for Epetra_Directory
//

/* get Epetra_Directory from non-const table using CT_Epetra_Directory_ID */
const Teuchos::RCP<Epetra_Directory>
CEpetra::getDirectory( CT_Epetra_Directory_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_Directory>(
        CTrilinos::abstractType<CT_Epetra_Directory_ID_t>(id));
}

/* get Epetra_Directory from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_Directory>
CEpetra::getDirectory( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_Directory>(id);
}

/* get const Epetra_Directory from either the const or non-const table
 * using CT_Epetra_Directory_ID */
const Teuchos::RCP<const Epetra_Directory>
CEpetra::getConstDirectory( CT_Epetra_Directory_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_Directory>(
        CTrilinos::abstractType<CT_Epetra_Directory_ID_t>(id));
}

/* get const Epetra_Directory from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_Directory>
CEpetra::getConstDirectory( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_Directory>(id);
}

/* store Epetra_Directory (owned) in non-const table */
CT_Epetra_Directory_ID_t
CEpetra::storeNewDirectory( Epetra_Directory *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Directory_ID_t>(
        CTrilinos::tableRepos().store<Epetra_Directory>(pobj, true));
}

/* store Epetra_Directory in non-const table */
CT_Epetra_Directory_ID_t
CEpetra::storeDirectory( Epetra_Directory *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Directory_ID_t>(
        CTrilinos::tableRepos().store<Epetra_Directory>(pobj, false));
}

/* store const Epetra_Directory in const table */
CT_Epetra_Directory_ID_t
CEpetra::storeConstDirectory( const Epetra_Directory *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Directory_ID_t>(
        CTrilinos::tableRepos().store<Epetra_Directory>(pobj, false));
}

/* remove Epetra_Directory from table using CT_Epetra_Directory_ID */
void
CEpetra::removeDirectory( CT_Epetra_Directory_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_Directory_ID_t>(*id);
    CTrilinos::tableRepos().remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_Directory_ID_t>(aid);
}

/* purge Epetra_Directory table */
void
CEpetra::purgeDirectory(  )
{
    CTrilinos::tableRepos().purge<Epetra_Directory>();
}


//
// Definitions for Epetra_Flops
//

/* get Epetra_Flops from non-const table using CT_Epetra_Flops_ID */
const Teuchos::RCP<Epetra_Flops>
CEpetra::getFlops( CT_Epetra_Flops_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_Flops>(
        CTrilinos::abstractType<CT_Epetra_Flops_ID_t>(id));
}

/* get Epetra_Flops from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_Flops>
CEpetra::getFlops( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_Flops>(id);
}

/* get const Epetra_Flops from either the const or non-const table
 * using CT_Epetra_Flops_ID */
const Teuchos::RCP<const Epetra_Flops>
CEpetra::getConstFlops( CT_Epetra_Flops_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_Flops>(
        CTrilinos::abstractType<CT_Epetra_Flops_ID_t>(id));
}

/* get const Epetra_Flops from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_Flops>
CEpetra::getConstFlops( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_Flops>(id);
}

/* store Epetra_Flops (owned) in non-const table */
CT_Epetra_Flops_ID_t
CEpetra::storeNewFlops( Epetra_Flops *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Flops_ID_t>(
        CTrilinos::tableRepos().store<Epetra_Flops>(pobj, true));
}

/* store Epetra_Flops in non-const table */
CT_Epetra_Flops_ID_t
CEpetra::storeFlops( Epetra_Flops *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Flops_ID_t>(
        CTrilinos::tableRepos().store<Epetra_Flops>(pobj, false));
}

/* store const Epetra_Flops in const table */
CT_Epetra_Flops_ID_t
CEpetra::storeConstFlops( const Epetra_Flops *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Flops_ID_t>(
        CTrilinos::tableRepos().store<Epetra_Flops>(pobj, false));
}

/* remove Epetra_Flops from table using CT_Epetra_Flops_ID */
void
CEpetra::removeFlops( CT_Epetra_Flops_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_Flops_ID_t>(*id);
    CTrilinos::tableRepos().remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_Flops_ID_t>(aid);
}

/* purge Epetra_Flops table */
void
CEpetra::purgeFlops(  )
{
    CTrilinos::tableRepos().purge<Epetra_Flops>();
}


//
// Definitions for Epetra_SrcDistObject
//

/* get Epetra_SrcDistObject from non-const table using CT_Epetra_SrcDistObject_ID */
const Teuchos::RCP<Epetra_SrcDistObject>
CEpetra::getSrcDistObject( CT_Epetra_SrcDistObject_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_SrcDistObject>(
        CTrilinos::abstractType<CT_Epetra_SrcDistObject_ID_t>(id));
}

/* get Epetra_SrcDistObject from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_SrcDistObject>
CEpetra::getSrcDistObject( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_SrcDistObject>(id);
}

/* get const Epetra_SrcDistObject from either the const or non-const table
 * using CT_Epetra_SrcDistObject_ID */
const Teuchos::RCP<const Epetra_SrcDistObject>
CEpetra::getConstSrcDistObject( CT_Epetra_SrcDistObject_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_SrcDistObject>(
        CTrilinos::abstractType<CT_Epetra_SrcDistObject_ID_t>(id));
}

/* get const Epetra_SrcDistObject from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_SrcDistObject>
CEpetra::getConstSrcDistObject( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_SrcDistObject>(id);
}

/* store Epetra_SrcDistObject (owned) in non-const table */
CT_Epetra_SrcDistObject_ID_t
CEpetra::storeNewSrcDistObject( Epetra_SrcDistObject *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_SrcDistObject_ID_t>(
        CTrilinos::tableRepos().store<Epetra_SrcDistObject>(pobj, true));
}

/* store Epetra_SrcDistObject in non-const table */
CT_Epetra_SrcDistObject_ID_t
CEpetra::storeSrcDistObject( Epetra_SrcDistObject *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_SrcDistObject_ID_t>(
        CTrilinos::tableRepos().store<Epetra_SrcDistObject>(pobj, false));
}

/* store const Epetra_SrcDistObject in const table */
CT_Epetra_SrcDistObject_ID_t
CEpetra::storeConstSrcDistObject( const Epetra_SrcDistObject *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_SrcDistObject_ID_t>(
        CTrilinos::tableRepos().store<Epetra_SrcDistObject>(pobj, false));
}

/* remove Epetra_SrcDistObject from table using CT_Epetra_SrcDistObject_ID */
void
CEpetra::removeSrcDistObject( CT_Epetra_SrcDistObject_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_SrcDistObject_ID_t>(*id);
    CTrilinos::tableRepos().remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_SrcDistObject_ID_t>(aid);
}

/* purge Epetra_SrcDistObject table */
void
CEpetra::purgeSrcDistObject(  )
{
    CTrilinos::tableRepos().purge<Epetra_SrcDistObject>();
}


//
// Definitions for Epetra_MpiComm
//

#ifdef HAVE_MPI

/* get Epetra_MpiComm from non-const table using CT_Epetra_MpiComm_ID */
const Teuchos::RCP<Epetra_MpiComm>
CEpetra::getMpiComm( CT_Epetra_MpiComm_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_MpiComm>(
        CTrilinos::abstractType<CT_Epetra_MpiComm_ID_t>(id));
}

/* get Epetra_MpiComm from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_MpiComm>
CEpetra::getMpiComm( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_MpiComm>(id);
}

/* get const Epetra_MpiComm from either the const or non-const table
 * using CT_Epetra_MpiComm_ID */
const Teuchos::RCP<const Epetra_MpiComm>
CEpetra::getConstMpiComm( CT_Epetra_MpiComm_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_MpiComm>(
        CTrilinos::abstractType<CT_Epetra_MpiComm_ID_t>(id));
}

/* get const Epetra_MpiComm from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_MpiComm>
CEpetra::getConstMpiComm( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_MpiComm>(id);
}

/* store Epetra_MpiComm (owned) in non-const table */
CT_Epetra_MpiComm_ID_t
CEpetra::storeNewMpiComm( Epetra_MpiComm *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_MpiComm_ID_t>(
        CTrilinos::tableRepos().store<Epetra_MpiComm>(pobj, true));
}

/* store Epetra_MpiComm in non-const table */
CT_Epetra_MpiComm_ID_t
CEpetra::storeMpiComm( Epetra_MpiComm *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_MpiComm_ID_t>(
        CTrilinos::tableRepos().store<Epetra_MpiComm>(pobj, false));
}

/* store const Epetra_MpiComm in const table */
CT_Epetra_MpiComm_ID_t
CEpetra::storeConstMpiComm( const Epetra_MpiComm *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_MpiComm_ID_t>(
        CTrilinos::tableRepos().store<Epetra_MpiComm>(pobj, false));
}

/* remove Epetra_MpiComm from table using CT_Epetra_MpiComm_ID */
void
CEpetra::removeMpiComm( CT_Epetra_MpiComm_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_MpiComm_ID_t>(*id);
    CTrilinos::tableRepos().remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_MpiComm_ID_t>(aid);
}

/* purge Epetra_MpiComm table */
void
CEpetra::purgeMpiComm(  )
{
    CTrilinos::tableRepos().purge<Epetra_MpiComm>();
}

#endif /* HAVE_MPI */


//
// Definitions for Epetra_CrsMatrix
//

/* get Epetra_CrsMatrix from non-const table using CT_Epetra_CrsMatrix_ID */
const Teuchos::RCP<Epetra_CrsMatrix>
CEpetra::getCrsMatrix( CT_Epetra_CrsMatrix_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_CrsMatrix>(
        CTrilinos::abstractType<CT_Epetra_CrsMatrix_ID_t>(id));
}

/* get Epetra_CrsMatrix from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_CrsMatrix>
CEpetra::getCrsMatrix( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_CrsMatrix>(id);
}

/* get const Epetra_CrsMatrix from either the const or non-const table
 * using CT_Epetra_CrsMatrix_ID */
const Teuchos::RCP<const Epetra_CrsMatrix>
CEpetra::getConstCrsMatrix( CT_Epetra_CrsMatrix_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_CrsMatrix>(
        CTrilinos::abstractType<CT_Epetra_CrsMatrix_ID_t>(id));
}

/* get const Epetra_CrsMatrix from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_CrsMatrix>
CEpetra::getConstCrsMatrix( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_CrsMatrix>(id);
}

/* store Epetra_CrsMatrix (owned) in non-const table */
CT_Epetra_CrsMatrix_ID_t
CEpetra::storeNewCrsMatrix( Epetra_CrsMatrix *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_CrsMatrix_ID_t>(
        CTrilinos::tableRepos().store<Epetra_CrsMatrix>(pobj, true));
}

/* store Epetra_CrsMatrix in non-const table */
CT_Epetra_CrsMatrix_ID_t
CEpetra::storeCrsMatrix( Epetra_CrsMatrix *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_CrsMatrix_ID_t>(
        CTrilinos::tableRepos().store<Epetra_CrsMatrix>(pobj, false));
}

/* store const Epetra_CrsMatrix in const table */
CT_Epetra_CrsMatrix_ID_t
CEpetra::storeConstCrsMatrix( const Epetra_CrsMatrix *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_CrsMatrix_ID_t>(
        CTrilinos::tableRepos().store<Epetra_CrsMatrix>(pobj, false));
}

/* remove Epetra_CrsMatrix from table using CT_Epetra_CrsMatrix_ID */
void
CEpetra::removeCrsMatrix( CT_Epetra_CrsMatrix_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_CrsMatrix_ID_t>(*id);
    CTrilinos::tableRepos().remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_CrsMatrix_ID_t>(aid);
}

/* purge Epetra_CrsMatrix table */
void
CEpetra::purgeCrsMatrix(  )
{
    CTrilinos::tableRepos().purge<Epetra_CrsMatrix>();
}


//
// Definitions for Epetra_CrsGraph
//

/* get Epetra_CrsGraph from non-const table using CT_Epetra_CrsGraph_ID */
const Teuchos::RCP<Epetra_CrsGraph>
CEpetra::getCrsGraph( CT_Epetra_CrsGraph_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_CrsGraph>(
        CTrilinos::abstractType<CT_Epetra_CrsGraph_ID_t>(id));
}

/* get Epetra_CrsGraph from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_CrsGraph>
CEpetra::getCrsGraph( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_CrsGraph>(id);
}

/* get const Epetra_CrsGraph from either the const or non-const table
 * using CT_Epetra_CrsGraph_ID */
const Teuchos::RCP<const Epetra_CrsGraph>
CEpetra::getConstCrsGraph( CT_Epetra_CrsGraph_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_CrsGraph>(
        CTrilinos::abstractType<CT_Epetra_CrsGraph_ID_t>(id));
}

/* get const Epetra_CrsGraph from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_CrsGraph>
CEpetra::getConstCrsGraph( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_CrsGraph>(id);
}

/* store Epetra_CrsGraph (owned) in non-const table */
CT_Epetra_CrsGraph_ID_t
CEpetra::storeNewCrsGraph( Epetra_CrsGraph *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_CrsGraph_ID_t>(
        CTrilinos::tableRepos().store<Epetra_CrsGraph>(pobj, true));
}

/* store Epetra_CrsGraph in non-const table */
CT_Epetra_CrsGraph_ID_t
CEpetra::storeCrsGraph( Epetra_CrsGraph *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_CrsGraph_ID_t>(
        CTrilinos::tableRepos().store<Epetra_CrsGraph>(pobj, false));
}

/* store const Epetra_CrsGraph in const table */
CT_Epetra_CrsGraph_ID_t
CEpetra::storeConstCrsGraph( const Epetra_CrsGraph *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_CrsGraph_ID_t>(
        CTrilinos::tableRepos().store<Epetra_CrsGraph>(pobj, false));
}

/* remove Epetra_CrsGraph from table using CT_Epetra_CrsGraph_ID */
void
CEpetra::removeCrsGraph( CT_Epetra_CrsGraph_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_CrsGraph_ID_t>(*id);
    CTrilinos::tableRepos().remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_CrsGraph_ID_t>(aid);
}

/* purge Epetra_CrsGraph table */
void
CEpetra::purgeCrsGraph(  )
{
    CTrilinos::tableRepos().purge<Epetra_CrsGraph>();
}


//
// Definitions for Epetra_DistObject
//

/* get Epetra_DistObject from non-const table using CT_Epetra_DistObject_ID */
const Teuchos::RCP<Epetra_DistObject>
CEpetra::getDistObject( CT_Epetra_DistObject_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_DistObject>(
        CTrilinos::abstractType<CT_Epetra_DistObject_ID_t>(id));
}

/* get Epetra_DistObject from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_DistObject>
CEpetra::getDistObject( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_DistObject>(id);
}

/* get const Epetra_DistObject from either the const or non-const table
 * using CT_Epetra_DistObject_ID */
const Teuchos::RCP<const Epetra_DistObject>
CEpetra::getConstDistObject( CT_Epetra_DistObject_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_DistObject>(
        CTrilinos::abstractType<CT_Epetra_DistObject_ID_t>(id));
}

/* get const Epetra_DistObject from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_DistObject>
CEpetra::getConstDistObject( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_DistObject>(id);
}

/* store Epetra_DistObject (owned) in non-const table */
CT_Epetra_DistObject_ID_t
CEpetra::storeNewDistObject( Epetra_DistObject *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_DistObject_ID_t>(
        CTrilinos::tableRepos().store<Epetra_DistObject>(pobj, true));
}

/* store Epetra_DistObject in non-const table */
CT_Epetra_DistObject_ID_t
CEpetra::storeDistObject( Epetra_DistObject *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_DistObject_ID_t>(
        CTrilinos::tableRepos().store<Epetra_DistObject>(pobj, false));
}

/* store const Epetra_DistObject in const table */
CT_Epetra_DistObject_ID_t
CEpetra::storeConstDistObject( const Epetra_DistObject *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_DistObject_ID_t>(
        CTrilinos::tableRepos().store<Epetra_DistObject>(pobj, false));
}

/* remove Epetra_DistObject from table using CT_Epetra_DistObject_ID */
void
CEpetra::removeDistObject( CT_Epetra_DistObject_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_DistObject_ID_t>(*id);
    CTrilinos::tableRepos().remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_DistObject_ID_t>(aid);
}

/* purge Epetra_DistObject table */
void
CEpetra::purgeDistObject(  )
{
    CTrilinos::tableRepos().purge<Epetra_DistObject>();
}


//
// Definitions for Epetra_Vector
//

/* get Epetra_Vector from non-const table using CT_Epetra_Vector_ID */
const Teuchos::RCP<Epetra_Vector>
CEpetra::getVector( CT_Epetra_Vector_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_Vector>(
        CTrilinos::abstractType<CT_Epetra_Vector_ID_t>(id));
}

/* get Epetra_Vector from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_Vector>
CEpetra::getVector( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_Vector>(id);
}

/* get const Epetra_Vector from either the const or non-const table
 * using CT_Epetra_Vector_ID */
const Teuchos::RCP<const Epetra_Vector>
CEpetra::getConstVector( CT_Epetra_Vector_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_Vector>(
        CTrilinos::abstractType<CT_Epetra_Vector_ID_t>(id));
}

/* get const Epetra_Vector from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_Vector>
CEpetra::getConstVector( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_Vector>(id);
}

/* store Epetra_Vector (owned) in non-const table */
CT_Epetra_Vector_ID_t
CEpetra::storeNewVector( Epetra_Vector *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Vector_ID_t>(
        CTrilinos::tableRepos().store<Epetra_Vector>(pobj, true));
}

/* store Epetra_Vector in non-const table */
CT_Epetra_Vector_ID_t
CEpetra::storeVector( Epetra_Vector *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Vector_ID_t>(
        CTrilinos::tableRepos().store<Epetra_Vector>(pobj, false));
}

/* store const Epetra_Vector in const table */
CT_Epetra_Vector_ID_t
CEpetra::storeConstVector( const Epetra_Vector *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Vector_ID_t>(
        CTrilinos::tableRepos().store<Epetra_Vector>(pobj, false));
}

/* remove Epetra_Vector from table using CT_Epetra_Vector_ID */
void
CEpetra::removeVector( CT_Epetra_Vector_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_Vector_ID_t>(*id);
    CTrilinos::tableRepos().remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_Vector_ID_t>(aid);
}

/* purge Epetra_Vector table */
void
CEpetra::purgeVector(  )
{
    CTrilinos::tableRepos().purge<Epetra_Vector>();
}


//
// Definitions for Epetra_Export
//

/* get Epetra_Export from non-const table using CT_Epetra_Export_ID */
const Teuchos::RCP<Epetra_Export>
CEpetra::getExport( CT_Epetra_Export_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_Export>(
        CTrilinos::abstractType<CT_Epetra_Export_ID_t>(id));
}

/* get Epetra_Export from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_Export>
CEpetra::getExport( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_Export>(id);
}

/* get const Epetra_Export from either the const or non-const table
 * using CT_Epetra_Export_ID */
const Teuchos::RCP<const Epetra_Export>
CEpetra::getConstExport( CT_Epetra_Export_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_Export>(
        CTrilinos::abstractType<CT_Epetra_Export_ID_t>(id));
}

/* get const Epetra_Export from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_Export>
CEpetra::getConstExport( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_Export>(id);
}

/* store Epetra_Export (owned) in non-const table */
CT_Epetra_Export_ID_t
CEpetra::storeNewExport( Epetra_Export *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Export_ID_t>(
        CTrilinos::tableRepos().store<Epetra_Export>(pobj, true));
}

/* store Epetra_Export in non-const table */
CT_Epetra_Export_ID_t
CEpetra::storeExport( Epetra_Export *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Export_ID_t>(
        CTrilinos::tableRepos().store<Epetra_Export>(pobj, false));
}

/* store const Epetra_Export in const table */
CT_Epetra_Export_ID_t
CEpetra::storeConstExport( const Epetra_Export *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Export_ID_t>(
        CTrilinos::tableRepos().store<Epetra_Export>(pobj, false));
}

/* remove Epetra_Export from table using CT_Epetra_Export_ID */
void
CEpetra::removeExport( CT_Epetra_Export_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_Export_ID_t>(*id);
    CTrilinos::tableRepos().remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_Export_ID_t>(aid);
}

/* purge Epetra_Export table */
void
CEpetra::purgeExport(  )
{
    CTrilinos::tableRepos().purge<Epetra_Export>();
}


//
// Definitions for Epetra_Map
//

/* get Epetra_Map from non-const table using CT_Epetra_Map_ID */
const Teuchos::RCP<Epetra_Map>
CEpetra::getMap( CT_Epetra_Map_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_Map>(
        CTrilinos::abstractType<CT_Epetra_Map_ID_t>(id));
}

/* get Epetra_Map from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_Map>
CEpetra::getMap( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_Map>(id);
}

/* get const Epetra_Map from either the const or non-const table
 * using CT_Epetra_Map_ID */
const Teuchos::RCP<const Epetra_Map>
CEpetra::getConstMap( CT_Epetra_Map_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_Map>(
        CTrilinos::abstractType<CT_Epetra_Map_ID_t>(id));
}

/* get const Epetra_Map from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_Map>
CEpetra::getConstMap( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_Map>(id);
}

/* store Epetra_Map (owned) in non-const table */
CT_Epetra_Map_ID_t
CEpetra::storeNewMap( Epetra_Map *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Map_ID_t>(
        CTrilinos::tableRepos().store<Epetra_Map>(pobj, true));
}

/* store Epetra_Map in non-const table */
CT_Epetra_Map_ID_t
CEpetra::storeMap( Epetra_Map *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Map_ID_t>(
        CTrilinos::tableRepos().store<Epetra_Map>(pobj, false));
}

/* store const Epetra_Map in const table */
CT_Epetra_Map_ID_t
CEpetra::storeConstMap( const Epetra_Map *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Map_ID_t>(
        CTrilinos::tableRepos().store<Epetra_Map>(pobj, false));
}

/* remove Epetra_Map from table using CT_Epetra_Map_ID */
void
CEpetra::removeMap( CT_Epetra_Map_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_Map_ID_t>(*id);
    CTrilinos::tableRepos().remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_Map_ID_t>(aid);
}

/* purge Epetra_Map table */
void
CEpetra::purgeMap(  )
{
    CTrilinos::tableRepos().purge<Epetra_Map>();
}


//
// Definitions for Epetra_BlockMap
//

/* get Epetra_BlockMap from non-const table using CT_Epetra_BlockMap_ID */
const Teuchos::RCP<Epetra_BlockMap>
CEpetra::getBlockMap( CT_Epetra_BlockMap_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_BlockMap>(
        CTrilinos::abstractType<CT_Epetra_BlockMap_ID_t>(id));
}

/* get Epetra_BlockMap from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_BlockMap>
CEpetra::getBlockMap( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_BlockMap>(id);
}

/* get const Epetra_BlockMap from either the const or non-const table
 * using CT_Epetra_BlockMap_ID */
const Teuchos::RCP<const Epetra_BlockMap>
CEpetra::getConstBlockMap( CT_Epetra_BlockMap_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_BlockMap>(
        CTrilinos::abstractType<CT_Epetra_BlockMap_ID_t>(id));
}

/* get const Epetra_BlockMap from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_BlockMap>
CEpetra::getConstBlockMap( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_BlockMap>(id);
}

/* store Epetra_BlockMap (owned) in non-const table */
CT_Epetra_BlockMap_ID_t
CEpetra::storeNewBlockMap( Epetra_BlockMap *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_BlockMap_ID_t>(
        CTrilinos::tableRepos().store<Epetra_BlockMap>(pobj, true));
}

/* store Epetra_BlockMap in non-const table */
CT_Epetra_BlockMap_ID_t
CEpetra::storeBlockMap( Epetra_BlockMap *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_BlockMap_ID_t>(
        CTrilinos::tableRepos().store<Epetra_BlockMap>(pobj, false));
}

/* store const Epetra_BlockMap in const table */
CT_Epetra_BlockMap_ID_t
CEpetra::storeConstBlockMap( const Epetra_BlockMap *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_BlockMap_ID_t>(
        CTrilinos::tableRepos().store<Epetra_BlockMap>(pobj, false));
}

/* remove Epetra_BlockMap from table using CT_Epetra_BlockMap_ID */
void
CEpetra::removeBlockMap( CT_Epetra_BlockMap_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_BlockMap_ID_t>(*id);
    CTrilinos::tableRepos().remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_BlockMap_ID_t>(aid);
}

/* purge Epetra_BlockMap table */
void
CEpetra::purgeBlockMap(  )
{
    CTrilinos::tableRepos().purge<Epetra_BlockMap>();
}


//
// Definitions for Epetra_Import
//

/* get Epetra_Import from non-const table using CT_Epetra_Import_ID */
const Teuchos::RCP<Epetra_Import>
CEpetra::getImport( CT_Epetra_Import_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_Import>(
        CTrilinos::abstractType<CT_Epetra_Import_ID_t>(id));
}

/* get Epetra_Import from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_Import>
CEpetra::getImport( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_Import>(id);
}

/* get const Epetra_Import from either the const or non-const table
 * using CT_Epetra_Import_ID */
const Teuchos::RCP<const Epetra_Import>
CEpetra::getConstImport( CT_Epetra_Import_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_Import>(
        CTrilinos::abstractType<CT_Epetra_Import_ID_t>(id));
}

/* get const Epetra_Import from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_Import>
CEpetra::getConstImport( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_Import>(id);
}

/* store Epetra_Import (owned) in non-const table */
CT_Epetra_Import_ID_t
CEpetra::storeNewImport( Epetra_Import *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Import_ID_t>(
        CTrilinos::tableRepos().store<Epetra_Import>(pobj, true));
}

/* store Epetra_Import in non-const table */
CT_Epetra_Import_ID_t
CEpetra::storeImport( Epetra_Import *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Import_ID_t>(
        CTrilinos::tableRepos().store<Epetra_Import>(pobj, false));
}

/* store const Epetra_Import in const table */
CT_Epetra_Import_ID_t
CEpetra::storeConstImport( const Epetra_Import *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Import_ID_t>(
        CTrilinos::tableRepos().store<Epetra_Import>(pobj, false));
}

/* remove Epetra_Import from table using CT_Epetra_Import_ID */
void
CEpetra::removeImport( CT_Epetra_Import_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_Import_ID_t>(*id);
    CTrilinos::tableRepos().remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_Import_ID_t>(aid);
}

/* purge Epetra_Import table */
void
CEpetra::purgeImport(  )
{
    CTrilinos::tableRepos().purge<Epetra_Import>();
}


//
// Definitions for Epetra_Time
//

/* get Epetra_Time from non-const table using CT_Epetra_Time_ID */
const Teuchos::RCP<Epetra_Time>
CEpetra::getTime( CT_Epetra_Time_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_Time>(
        CTrilinos::abstractType<CT_Epetra_Time_ID_t>(id));
}

/* get Epetra_Time from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_Time>
CEpetra::getTime( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_Time>(id);
}

/* get const Epetra_Time from either the const or non-const table
 * using CT_Epetra_Time_ID */
const Teuchos::RCP<const Epetra_Time>
CEpetra::getConstTime( CT_Epetra_Time_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_Time>(
        CTrilinos::abstractType<CT_Epetra_Time_ID_t>(id));
}

/* get const Epetra_Time from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_Time>
CEpetra::getConstTime( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_Time>(id);
}

/* store Epetra_Time (owned) in non-const table */
CT_Epetra_Time_ID_t
CEpetra::storeNewTime( Epetra_Time *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Time_ID_t>(
        CTrilinos::tableRepos().store<Epetra_Time>(pobj, true));
}

/* store Epetra_Time in non-const table */
CT_Epetra_Time_ID_t
CEpetra::storeTime( Epetra_Time *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Time_ID_t>(
        CTrilinos::tableRepos().store<Epetra_Time>(pobj, false));
}

/* store const Epetra_Time in const table */
CT_Epetra_Time_ID_t
CEpetra::storeConstTime( const Epetra_Time *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Time_ID_t>(
        CTrilinos::tableRepos().store<Epetra_Time>(pobj, false));
}

/* remove Epetra_Time from table using CT_Epetra_Time_ID */
void
CEpetra::removeTime( CT_Epetra_Time_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_Time_ID_t>(*id);
    CTrilinos::tableRepos().remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_Time_ID_t>(aid);
}

/* purge Epetra_Time table */
void
CEpetra::purgeTime(  )
{
    CTrilinos::tableRepos().purge<Epetra_Time>();
}


//
// Definitions for Epetra_JadMatrix
//

/* get Epetra_JadMatrix from non-const table using CT_Epetra_JadMatrix_ID */
const Teuchos::RCP<Epetra_JadMatrix>
CEpetra::getJadMatrix( CT_Epetra_JadMatrix_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_JadMatrix>(
        CTrilinos::abstractType<CT_Epetra_JadMatrix_ID_t>(id));
}

/* get Epetra_JadMatrix from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_JadMatrix>
CEpetra::getJadMatrix( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_JadMatrix>(id);
}

/* get const Epetra_JadMatrix from either the const or non-const table
 * using CT_Epetra_JadMatrix_ID */
const Teuchos::RCP<const Epetra_JadMatrix>
CEpetra::getConstJadMatrix( CT_Epetra_JadMatrix_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_JadMatrix>(
        CTrilinos::abstractType<CT_Epetra_JadMatrix_ID_t>(id));
}

/* get const Epetra_JadMatrix from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_JadMatrix>
CEpetra::getConstJadMatrix( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_JadMatrix>(id);
}

/* store Epetra_JadMatrix (owned) in non-const table */
CT_Epetra_JadMatrix_ID_t
CEpetra::storeNewJadMatrix( Epetra_JadMatrix *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_JadMatrix_ID_t>(
        CTrilinos::tableRepos().store<Epetra_JadMatrix>(pobj, true));
}

/* store Epetra_JadMatrix in non-const table */
CT_Epetra_JadMatrix_ID_t
CEpetra::storeJadMatrix( Epetra_JadMatrix *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_JadMatrix_ID_t>(
        CTrilinos::tableRepos().store<Epetra_JadMatrix>(pobj, false));
}

/* store const Epetra_JadMatrix in const table */
CT_Epetra_JadMatrix_ID_t
CEpetra::storeConstJadMatrix( const Epetra_JadMatrix *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_JadMatrix_ID_t>(
        CTrilinos::tableRepos().store<Epetra_JadMatrix>(pobj, false));
}

/* remove Epetra_JadMatrix from table using CT_Epetra_JadMatrix_ID */
void
CEpetra::removeJadMatrix( CT_Epetra_JadMatrix_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_JadMatrix_ID_t>(*id);
    CTrilinos::tableRepos().remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_JadMatrix_ID_t>(aid);
}

/* purge Epetra_JadMatrix table */
void
CEpetra::purgeJadMatrix(  )
{
    CTrilinos::tableRepos().purge<Epetra_JadMatrix>();
}


//
// Definitions for Epetra_LinearProblem
//

/* get Epetra_LinearProblem from non-const table using CT_Epetra_LinearProblem_ID */
const Teuchos::RCP<Epetra_LinearProblem>
CEpetra::getLinearProblem( CT_Epetra_LinearProblem_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_LinearProblem>(
        CTrilinos::abstractType<CT_Epetra_LinearProblem_ID_t>(id));
}

/* get Epetra_LinearProblem from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_LinearProblem>
CEpetra::getLinearProblem( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_LinearProblem>(id);
}

/* get const Epetra_LinearProblem from either the const or non-const table
 * using CT_Epetra_LinearProblem_ID */
const Teuchos::RCP<const Epetra_LinearProblem>
CEpetra::getConstLinearProblem( CT_Epetra_LinearProblem_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_LinearProblem>(
        CTrilinos::abstractType<CT_Epetra_LinearProblem_ID_t>(id));
}

/* get const Epetra_LinearProblem from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_LinearProblem>
CEpetra::getConstLinearProblem( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_LinearProblem>(id);
}

/* store Epetra_LinearProblem (owned) in non-const table */
CT_Epetra_LinearProblem_ID_t
CEpetra::storeNewLinearProblem( Epetra_LinearProblem *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_LinearProblem_ID_t>(
        CTrilinos::tableRepos().store<Epetra_LinearProblem>(pobj, true));
}

/* store Epetra_LinearProblem in non-const table */
CT_Epetra_LinearProblem_ID_t
CEpetra::storeLinearProblem( Epetra_LinearProblem *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_LinearProblem_ID_t>(
        CTrilinos::tableRepos().store<Epetra_LinearProblem>(pobj, false));
}

/* store const Epetra_LinearProblem in const table */
CT_Epetra_LinearProblem_ID_t
CEpetra::storeConstLinearProblem( const Epetra_LinearProblem *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_LinearProblem_ID_t>(
        CTrilinos::tableRepos().store<Epetra_LinearProblem>(pobj, false));
}

/* remove Epetra_LinearProblem from table using CT_Epetra_LinearProblem_ID */
void
CEpetra::removeLinearProblem( CT_Epetra_LinearProblem_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_LinearProblem_ID_t>(*id);
    CTrilinos::tableRepos().remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_LinearProblem_ID_t>(aid);
}

/* purge Epetra_LinearProblem table */
void
CEpetra::purgeLinearProblem(  )
{
    CTrilinos::tableRepos().purge<Epetra_LinearProblem>();
}


//
// Definitions for Epetra_LAPACK
//

/* get Epetra_LAPACK from non-const table using CT_Epetra_LAPACK_ID */
const Teuchos::RCP<Epetra_LAPACK>
CEpetra::getLAPACK( CT_Epetra_LAPACK_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_LAPACK>(
        CTrilinos::abstractType<CT_Epetra_LAPACK_ID_t>(id));
}

/* get Epetra_LAPACK from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_LAPACK>
CEpetra::getLAPACK( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_LAPACK>(id);
}

/* get const Epetra_LAPACK from either the const or non-const table
 * using CT_Epetra_LAPACK_ID */
const Teuchos::RCP<const Epetra_LAPACK>
CEpetra::getConstLAPACK( CT_Epetra_LAPACK_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_LAPACK>(
        CTrilinos::abstractType<CT_Epetra_LAPACK_ID_t>(id));
}

/* get const Epetra_LAPACK from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_LAPACK>
CEpetra::getConstLAPACK( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_LAPACK>(id);
}

/* store Epetra_LAPACK (owned) in non-const table */
CT_Epetra_LAPACK_ID_t
CEpetra::storeNewLAPACK( Epetra_LAPACK *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_LAPACK_ID_t>(
        CTrilinos::tableRepos().store<Epetra_LAPACK>(pobj, true));
}

/* store Epetra_LAPACK in non-const table */
CT_Epetra_LAPACK_ID_t
CEpetra::storeLAPACK( Epetra_LAPACK *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_LAPACK_ID_t>(
        CTrilinos::tableRepos().store<Epetra_LAPACK>(pobj, false));
}

/* store const Epetra_LAPACK in const table */
CT_Epetra_LAPACK_ID_t
CEpetra::storeConstLAPACK( const Epetra_LAPACK *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_LAPACK_ID_t>(
        CTrilinos::tableRepos().store<Epetra_LAPACK>(pobj, false));
}

/* remove Epetra_LAPACK from table using CT_Epetra_LAPACK_ID */
void
CEpetra::removeLAPACK( CT_Epetra_LAPACK_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_LAPACK_ID_t>(*id);
    CTrilinos::tableRepos().remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_LAPACK_ID_t>(aid);
}

/* purge Epetra_LAPACK table */
void
CEpetra::purgeLAPACK(  )
{
    CTrilinos::tableRepos().purge<Epetra_LAPACK>();
}


//
// Definitions for ParameterList
//

/* get Teuchos::ParameterList from non-const table using CT_Teuchos_ParameterList_ID */
const Teuchos::RCP<Teuchos::ParameterList>
CTeuchos::getParameterList( CT_Teuchos_ParameterList_ID_t id )
{
    return CTrilinos::tableRepos().get<Teuchos::ParameterList>(
        CTrilinos::abstractType<CT_Teuchos_ParameterList_ID_t>(id));
}

/* get Teuchos::ParameterList from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Teuchos::ParameterList>
CTeuchos::getParameterList( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().get<Teuchos::ParameterList>(id);
}

/* get const Teuchos::ParameterList from either the const or non-const table
 * using CT_Teuchos_ParameterList_ID */
const Teuchos::RCP<const Teuchos::ParameterList>
CTeuchos::getConstParameterList( CT_Teuchos_ParameterList_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Teuchos::ParameterList>(
        CTrilinos::abstractType<CT_Teuchos_ParameterList_ID_t>(id));
}

/* get const Teuchos::ParameterList from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Teuchos::ParameterList>
CTeuchos::getConstParameterList( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Teuchos::ParameterList>(id);
}

/* store Teuchos::ParameterList (owned) in non-const table */
CT_Teuchos_ParameterList_ID_t
CTeuchos::storeNewParameterList( Teuchos::ParameterList *pobj )
{
    return CTrilinos::concreteType<CT_Teuchos_ParameterList_ID_t>(
        CTrilinos::tableRepos().store<Teuchos::ParameterList>(pobj, true));
}

/* store Teuchos::ParameterList in non-const table */
CT_Teuchos_ParameterList_ID_t
CTeuchos::storeParameterList( Teuchos::ParameterList *pobj )
{
    return CTrilinos::concreteType<CT_Teuchos_ParameterList_ID_t>(
        CTrilinos::tableRepos().store<Teuchos::ParameterList>(pobj, false));
}

/* store const Teuchos::ParameterList in const table */
CT_Teuchos_ParameterList_ID_t
CTeuchos::storeConstParameterList( const Teuchos::ParameterList *pobj )
{
    return CTrilinos::concreteType<CT_Teuchos_ParameterList_ID_t>(
        CTrilinos::tableRepos().store<Teuchos::ParameterList>(pobj, false));
}

/* remove Teuchos::ParameterList from table using CT_Teuchos_ParameterList_ID */
void
CTeuchos::removeParameterList( CT_Teuchos_ParameterList_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Teuchos_ParameterList_ID_t>(*id);
    CTrilinos::tableRepos().remove(&aid);
    *id = CTrilinos::concreteType<CT_Teuchos_ParameterList_ID_t>(aid);
}

/* purge Teuchos::ParameterList table */
void
CTeuchos::purgeParameterList(  )
{
    CTrilinos::tableRepos().purge<Teuchos::ParameterList>();
}


//
// Definitions for Amesos_BaseSolver
//

#ifdef HAVE_CTRILINOS_AMESOS

/* get Amesos_BaseSolver from non-const table using CT_Amesos_BaseSolver_ID */
const Teuchos::RCP<Amesos_BaseSolver>
CAmesos::getBaseSolver( CT_Amesos_BaseSolver_ID_t id )
{
    return CTrilinos::tableRepos().get<Amesos_BaseSolver>(
        CTrilinos::abstractType<CT_Amesos_BaseSolver_ID_t>(id));
}

/* get Amesos_BaseSolver from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Amesos_BaseSolver>
CAmesos::getBaseSolver( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().get<Amesos_BaseSolver>(id);
}

/* get const Amesos_BaseSolver from either the const or non-const table
 * using CT_Amesos_BaseSolver_ID */
const Teuchos::RCP<const Amesos_BaseSolver>
CAmesos::getConstBaseSolver( CT_Amesos_BaseSolver_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Amesos_BaseSolver>(
        CTrilinos::abstractType<CT_Amesos_BaseSolver_ID_t>(id));
}

/* get const Amesos_BaseSolver from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Amesos_BaseSolver>
CAmesos::getConstBaseSolver( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Amesos_BaseSolver>(id);
}

/* store Amesos_BaseSolver (owned) in non-const table */
CT_Amesos_BaseSolver_ID_t
CAmesos::storeNewBaseSolver( Amesos_BaseSolver *pobj )
{
    return CTrilinos::concreteType<CT_Amesos_BaseSolver_ID_t>(
        CTrilinos::tableRepos().store<Amesos_BaseSolver>(pobj, true));
}

/* store Amesos_BaseSolver in non-const table */
CT_Amesos_BaseSolver_ID_t
CAmesos::storeBaseSolver( Amesos_BaseSolver *pobj )
{
    return CTrilinos::concreteType<CT_Amesos_BaseSolver_ID_t>(
        CTrilinos::tableRepos().store<Amesos_BaseSolver>(pobj, false));
}

/* store const Amesos_BaseSolver in const table */
CT_Amesos_BaseSolver_ID_t
CAmesos::storeConstBaseSolver( const Amesos_BaseSolver *pobj )
{
    return CTrilinos::concreteType<CT_Amesos_BaseSolver_ID_t>(
        CTrilinos::tableRepos().store<Amesos_BaseSolver>(pobj, false));
}

/* remove Amesos_BaseSolver from table using CT_Amesos_BaseSolver_ID */
void
CAmesos::removeBaseSolver( CT_Amesos_BaseSolver_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Amesos_BaseSolver_ID_t>(*id);
    CTrilinos::tableRepos().remove(&aid);
    *id = CTrilinos::concreteType<CT_Amesos_BaseSolver_ID_t>(aid);
}

/* purge Amesos_BaseSolver table */
void
CAmesos::purgeBaseSolver(  )
{
    CTrilinos::tableRepos().purge<Amesos_BaseSolver>();
}

#endif /* HAVE_CTRILINOS_AMESOS */


//
// Definitions for Epetra_FECrsMatrix
//

/* get Epetra_FECrsMatrix from non-const table using CT_Epetra_FECrsMatrix_ID */
const Teuchos::RCP<Epetra_FECrsMatrix>
CEpetra::getFECrsMatrix( CT_Epetra_FECrsMatrix_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_FECrsMatrix>(
        CTrilinos::abstractType<CT_Epetra_FECrsMatrix_ID_t>(id));
}

/* get Epetra_FECrsMatrix from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_FECrsMatrix>
CEpetra::getFECrsMatrix( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_FECrsMatrix>(id);
}

/* get const Epetra_FECrsMatrix from either the const or non-const table
 * using CT_Epetra_FECrsMatrix_ID */
const Teuchos::RCP<const Epetra_FECrsMatrix>
CEpetra::getConstFECrsMatrix( CT_Epetra_FECrsMatrix_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_FECrsMatrix>(
        CTrilinos::abstractType<CT_Epetra_FECrsMatrix_ID_t>(id));
}

/* get const Epetra_FECrsMatrix from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_FECrsMatrix>
CEpetra::getConstFECrsMatrix( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_FECrsMatrix>(id);
}

/* store Epetra_FECrsMatrix (owned) in non-const table */
CT_Epetra_FECrsMatrix_ID_t
CEpetra::storeNewFECrsMatrix( Epetra_FECrsMatrix *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_FECrsMatrix_ID_t>(
        CTrilinos::tableRepos().store<Epetra_FECrsMatrix>(pobj, true));
}

/* store Epetra_FECrsMatrix in non-const table */
CT_Epetra_FECrsMatrix_ID_t
CEpetra::storeFECrsMatrix( Epetra_FECrsMatrix *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_FECrsMatrix_ID_t>(
        CTrilinos::tableRepos().store<Epetra_FECrsMatrix>(pobj, false));
}

/* store const Epetra_FECrsMatrix in const table */
CT_Epetra_FECrsMatrix_ID_t
CEpetra::storeConstFECrsMatrix( const Epetra_FECrsMatrix *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_FECrsMatrix_ID_t>(
        CTrilinos::tableRepos().store<Epetra_FECrsMatrix>(pobj, false));
}

/* remove Epetra_FECrsMatrix from table using CT_Epetra_FECrsMatrix_ID */
void
CEpetra::removeFECrsMatrix( CT_Epetra_FECrsMatrix_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_FECrsMatrix_ID_t>(*id);
    CTrilinos::tableRepos().remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_FECrsMatrix_ID_t>(aid);
}

/* purge Epetra_FECrsMatrix table */
void
CEpetra::purgeFECrsMatrix(  )
{
    CTrilinos::tableRepos().purge<Epetra_FECrsMatrix>();
}


//
// Definitions for Epetra_IntSerialDenseVector
//

/* get Epetra_IntSerialDenseVector from non-const table using CT_Epetra_IntSerialDenseVector_ID */
const Teuchos::RCP<Epetra_IntSerialDenseVector>
CEpetra::getIntSerialDenseVector( CT_Epetra_IntSerialDenseVector_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_IntSerialDenseVector>(
        CTrilinos::abstractType<CT_Epetra_IntSerialDenseVector_ID_t>(id));
}

/* get Epetra_IntSerialDenseVector from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_IntSerialDenseVector>
CEpetra::getIntSerialDenseVector( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_IntSerialDenseVector>(id);
}

/* get const Epetra_IntSerialDenseVector from either the const or non-const table
 * using CT_Epetra_IntSerialDenseVector_ID */
const Teuchos::RCP<const Epetra_IntSerialDenseVector>
CEpetra::getConstIntSerialDenseVector( CT_Epetra_IntSerialDenseVector_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_IntSerialDenseVector>(
        CTrilinos::abstractType<CT_Epetra_IntSerialDenseVector_ID_t>(id));
}

/* get const Epetra_IntSerialDenseVector from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_IntSerialDenseVector>
CEpetra::getConstIntSerialDenseVector( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_IntSerialDenseVector>(id);
}

/* store Epetra_IntSerialDenseVector (owned) in non-const table */
CT_Epetra_IntSerialDenseVector_ID_t
CEpetra::storeNewIntSerialDenseVector( Epetra_IntSerialDenseVector *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_IntSerialDenseVector_ID_t>(
        CTrilinos::tableRepos().store<Epetra_IntSerialDenseVector>(pobj, true));
}

/* store Epetra_IntSerialDenseVector in non-const table */
CT_Epetra_IntSerialDenseVector_ID_t
CEpetra::storeIntSerialDenseVector( Epetra_IntSerialDenseVector *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_IntSerialDenseVector_ID_t>(
        CTrilinos::tableRepos().store<Epetra_IntSerialDenseVector>(pobj, false));
}

/* store const Epetra_IntSerialDenseVector in const table */
CT_Epetra_IntSerialDenseVector_ID_t
CEpetra::storeConstIntSerialDenseVector( const Epetra_IntSerialDenseVector *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_IntSerialDenseVector_ID_t>(
        CTrilinos::tableRepos().store<Epetra_IntSerialDenseVector>(pobj, false));
}

/* remove Epetra_IntSerialDenseVector from table using CT_Epetra_IntSerialDenseVector_ID */
void
CEpetra::removeIntSerialDenseVector( CT_Epetra_IntSerialDenseVector_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_IntSerialDenseVector_ID_t>(*id);
    CTrilinos::tableRepos().remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_IntSerialDenseVector_ID_t>(aid);
}

/* purge Epetra_IntSerialDenseVector table */
void
CEpetra::purgeIntSerialDenseVector(  )
{
    CTrilinos::tableRepos().purge<Epetra_IntSerialDenseVector>();
}


//
// Definitions for Epetra_SerialDenseMatrix
//

/* get Epetra_SerialDenseMatrix from non-const table using CT_Epetra_SerialDenseMatrix_ID */
const Teuchos::RCP<Epetra_SerialDenseMatrix>
CEpetra::getSerialDenseMatrix( CT_Epetra_SerialDenseMatrix_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_SerialDenseMatrix>(
        CTrilinos::abstractType<CT_Epetra_SerialDenseMatrix_ID_t>(id));
}

/* get Epetra_SerialDenseMatrix from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_SerialDenseMatrix>
CEpetra::getSerialDenseMatrix( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().get<Epetra_SerialDenseMatrix>(id);
}

/* get const Epetra_SerialDenseMatrix from either the const or non-const table
 * using CT_Epetra_SerialDenseMatrix_ID */
const Teuchos::RCP<const Epetra_SerialDenseMatrix>
CEpetra::getConstSerialDenseMatrix( CT_Epetra_SerialDenseMatrix_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_SerialDenseMatrix>(
        CTrilinos::abstractType<CT_Epetra_SerialDenseMatrix_ID_t>(id));
}

/* get const Epetra_SerialDenseMatrix from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_SerialDenseMatrix>
CEpetra::getConstSerialDenseMatrix( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Epetra_SerialDenseMatrix>(id);
}

/* store Epetra_SerialDenseMatrix (owned) in non-const table */
CT_Epetra_SerialDenseMatrix_ID_t
CEpetra::storeNewSerialDenseMatrix( Epetra_SerialDenseMatrix *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_SerialDenseMatrix_ID_t>(
        CTrilinos::tableRepos().store<Epetra_SerialDenseMatrix>(pobj, true));
}

/* store Epetra_SerialDenseMatrix in non-const table */
CT_Epetra_SerialDenseMatrix_ID_t
CEpetra::storeSerialDenseMatrix( Epetra_SerialDenseMatrix *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_SerialDenseMatrix_ID_t>(
        CTrilinos::tableRepos().store<Epetra_SerialDenseMatrix>(pobj, false));
}

/* store const Epetra_SerialDenseMatrix in const table */
CT_Epetra_SerialDenseMatrix_ID_t
CEpetra::storeConstSerialDenseMatrix( const Epetra_SerialDenseMatrix *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_SerialDenseMatrix_ID_t>(
        CTrilinos::tableRepos().store<Epetra_SerialDenseMatrix>(pobj, false));
}

/* remove Epetra_SerialDenseMatrix from table using CT_Epetra_SerialDenseMatrix_ID */
void
CEpetra::removeSerialDenseMatrix( CT_Epetra_SerialDenseMatrix_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_SerialDenseMatrix_ID_t>(*id);
    CTrilinos::tableRepos().remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_SerialDenseMatrix_ID_t>(aid);
}

/* purge Epetra_SerialDenseMatrix table */
void
CEpetra::purgeSerialDenseMatrix(  )
{
    CTrilinos::tableRepos().purge<Epetra_SerialDenseMatrix>();
}


//
// Definitions for AztecOO_StatusTest
//

#ifdef HAVE_CTRILINOS_AZTECOO

/* get AztecOO_StatusTest from non-const table using CT_AztecOO_StatusTest_ID */
const Teuchos::RCP<AztecOO_StatusTest>
CAztecOO::getStatusTest( CT_AztecOO_StatusTest_ID_t id )
{
    return CTrilinos::tableRepos().get<AztecOO_StatusTest>(
        CTrilinos::abstractType<CT_AztecOO_StatusTest_ID_t>(id));
}

/* get AztecOO_StatusTest from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<AztecOO_StatusTest>
CAztecOO::getStatusTest( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().get<AztecOO_StatusTest>(id);
}

/* get const AztecOO_StatusTest from either the const or non-const table
 * using CT_AztecOO_StatusTest_ID */
const Teuchos::RCP<const AztecOO_StatusTest>
CAztecOO::getConstStatusTest( CT_AztecOO_StatusTest_ID_t id )
{
    return CTrilinos::tableRepos().getConst<AztecOO_StatusTest>(
        CTrilinos::abstractType<CT_AztecOO_StatusTest_ID_t>(id));
}

/* get const AztecOO_StatusTest from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const AztecOO_StatusTest>
CAztecOO::getConstStatusTest( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().getConst<AztecOO_StatusTest>(id);
}

/* store AztecOO_StatusTest (owned) in non-const table */
CT_AztecOO_StatusTest_ID_t
CAztecOO::storeNewStatusTest( AztecOO_StatusTest *pobj )
{
    return CTrilinos::concreteType<CT_AztecOO_StatusTest_ID_t>(
        CTrilinos::tableRepos().store<AztecOO_StatusTest>(pobj, true));
}

/* store AztecOO_StatusTest in non-const table */
CT_AztecOO_StatusTest_ID_t
CAztecOO::storeStatusTest( AztecOO_StatusTest *pobj )
{
    return CTrilinos::concreteType<CT_AztecOO_StatusTest_ID_t>(
        CTrilinos::tableRepos().store<AztecOO_StatusTest>(pobj, false));
}

/* store const AztecOO_StatusTest in const table */
CT_AztecOO_StatusTest_ID_t
CAztecOO::storeConstStatusTest( const AztecOO_StatusTest *pobj )
{
    return CTrilinos::concreteType<CT_AztecOO_StatusTest_ID_t>(
        CTrilinos::tableRepos().store<AztecOO_StatusTest>(pobj, false));
}

/* remove AztecOO_StatusTest from table using CT_AztecOO_StatusTest_ID */
void
CAztecOO::removeStatusTest( CT_AztecOO_StatusTest_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_AztecOO_StatusTest_ID_t>(*id);
    CTrilinos::tableRepos().remove(&aid);
    *id = CTrilinos::concreteType<CT_AztecOO_StatusTest_ID_t>(aid);
}

/* purge AztecOO_StatusTest table */
void
CAztecOO::purgeStatusTest(  )
{
    CTrilinos::tableRepos().purge<AztecOO_StatusTest>();
}

#endif /* HAVE_CTRILINOS_AZTECOO */


//
// Definitions for AztecOO_StatusTestCombo
//

#ifdef HAVE_CTRILINOS_AZTECOO

/* get AztecOO_StatusTestCombo from non-const table using CT_AztecOO_StatusTestCombo_ID */
const Teuchos::RCP<AztecOO_StatusTestCombo>
CAztecOO::getStatusTestCombo( CT_AztecOO_StatusTestCombo_ID_t id )
{
    return CTrilinos::tableRepos().get<AztecOO_StatusTestCombo>(
        CTrilinos::abstractType<CT_AztecOO_StatusTestCombo_ID_t>(id));
}

/* get AztecOO_StatusTestCombo from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<AztecOO_StatusTestCombo>
CAztecOO::getStatusTestCombo( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().get<AztecOO_StatusTestCombo>(id);
}

/* get const AztecOO_StatusTestCombo from either the const or non-const table
 * using CT_AztecOO_StatusTestCombo_ID */
const Teuchos::RCP<const AztecOO_StatusTestCombo>
CAztecOO::getConstStatusTestCombo( CT_AztecOO_StatusTestCombo_ID_t id )
{
    return CTrilinos::tableRepos().getConst<AztecOO_StatusTestCombo>(
        CTrilinos::abstractType<CT_AztecOO_StatusTestCombo_ID_t>(id));
}

/* get const AztecOO_StatusTestCombo from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const AztecOO_StatusTestCombo>
CAztecOO::getConstStatusTestCombo( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().getConst<AztecOO_StatusTestCombo>(id);
}

/* store AztecOO_StatusTestCombo (owned) in non-const table */
CT_AztecOO_StatusTestCombo_ID_t
CAztecOO::storeNewStatusTestCombo( AztecOO_StatusTestCombo *pobj )
{
    return CTrilinos::concreteType<CT_AztecOO_StatusTestCombo_ID_t>(
        CTrilinos::tableRepos().store<AztecOO_StatusTestCombo>(pobj, true));
}

/* store AztecOO_StatusTestCombo in non-const table */
CT_AztecOO_StatusTestCombo_ID_t
CAztecOO::storeStatusTestCombo( AztecOO_StatusTestCombo *pobj )
{
    return CTrilinos::concreteType<CT_AztecOO_StatusTestCombo_ID_t>(
        CTrilinos::tableRepos().store<AztecOO_StatusTestCombo>(pobj, false));
}

/* store const AztecOO_StatusTestCombo in const table */
CT_AztecOO_StatusTestCombo_ID_t
CAztecOO::storeConstStatusTestCombo( const AztecOO_StatusTestCombo *pobj )
{
    return CTrilinos::concreteType<CT_AztecOO_StatusTestCombo_ID_t>(
        CTrilinos::tableRepos().store<AztecOO_StatusTestCombo>(pobj, false));
}

/* remove AztecOO_StatusTestCombo from table using CT_AztecOO_StatusTestCombo_ID */
void
CAztecOO::removeStatusTestCombo( CT_AztecOO_StatusTestCombo_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_AztecOO_StatusTestCombo_ID_t>(*id);
    CTrilinos::tableRepos().remove(&aid);
    *id = CTrilinos::concreteType<CT_AztecOO_StatusTestCombo_ID_t>(aid);
}

/* purge AztecOO_StatusTestCombo table */
void
CAztecOO::purgeStatusTestCombo(  )
{
    CTrilinos::tableRepos().purge<AztecOO_StatusTestCombo>();
}

#endif /* HAVE_CTRILINOS_AZTECOO */


//
// Definitions for AztecOO_StatusTestMaxIters
//

#ifdef HAVE_CTRILINOS_AZTECOO

/* get AztecOO_StatusTestMaxIters from non-const table using CT_AztecOO_StatusTestMaxIters_ID */
const Teuchos::RCP<AztecOO_StatusTestMaxIters>
CAztecOO::getStatusTestMaxIters( CT_AztecOO_StatusTestMaxIters_ID_t id )
{
    return CTrilinos::tableRepos().get<AztecOO_StatusTestMaxIters>(
        CTrilinos::abstractType<CT_AztecOO_StatusTestMaxIters_ID_t>(id));
}

/* get AztecOO_StatusTestMaxIters from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<AztecOO_StatusTestMaxIters>
CAztecOO::getStatusTestMaxIters( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().get<AztecOO_StatusTestMaxIters>(id);
}

/* get const AztecOO_StatusTestMaxIters from either the const or non-const table
 * using CT_AztecOO_StatusTestMaxIters_ID */
const Teuchos::RCP<const AztecOO_StatusTestMaxIters>
CAztecOO::getConstStatusTestMaxIters( CT_AztecOO_StatusTestMaxIters_ID_t id )
{
    return CTrilinos::tableRepos().getConst<AztecOO_StatusTestMaxIters>(
        CTrilinos::abstractType<CT_AztecOO_StatusTestMaxIters_ID_t>(id));
}

/* get const AztecOO_StatusTestMaxIters from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const AztecOO_StatusTestMaxIters>
CAztecOO::getConstStatusTestMaxIters( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().getConst<AztecOO_StatusTestMaxIters>(id);
}

/* store AztecOO_StatusTestMaxIters (owned) in non-const table */
CT_AztecOO_StatusTestMaxIters_ID_t
CAztecOO::storeNewStatusTestMaxIters( AztecOO_StatusTestMaxIters *pobj )
{
    return CTrilinos::concreteType<CT_AztecOO_StatusTestMaxIters_ID_t>(
        CTrilinos::tableRepos().store<AztecOO_StatusTestMaxIters>(pobj, true));
}

/* store AztecOO_StatusTestMaxIters in non-const table */
CT_AztecOO_StatusTestMaxIters_ID_t
CAztecOO::storeStatusTestMaxIters( AztecOO_StatusTestMaxIters *pobj )
{
    return CTrilinos::concreteType<CT_AztecOO_StatusTestMaxIters_ID_t>(
        CTrilinos::tableRepos().store<AztecOO_StatusTestMaxIters>(pobj, false));
}

/* store const AztecOO_StatusTestMaxIters in const table */
CT_AztecOO_StatusTestMaxIters_ID_t
CAztecOO::storeConstStatusTestMaxIters( const AztecOO_StatusTestMaxIters *pobj )
{
    return CTrilinos::concreteType<CT_AztecOO_StatusTestMaxIters_ID_t>(
        CTrilinos::tableRepos().store<AztecOO_StatusTestMaxIters>(pobj, false));
}

/* remove AztecOO_StatusTestMaxIters from table using CT_AztecOO_StatusTestMaxIters_ID */
void
CAztecOO::removeStatusTestMaxIters( CT_AztecOO_StatusTestMaxIters_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_AztecOO_StatusTestMaxIters_ID_t>(*id);
    CTrilinos::tableRepos().remove(&aid);
    *id = CTrilinos::concreteType<CT_AztecOO_StatusTestMaxIters_ID_t>(aid);
}

/* purge AztecOO_StatusTestMaxIters table */
void
CAztecOO::purgeStatusTestMaxIters(  )
{
    CTrilinos::tableRepos().purge<AztecOO_StatusTestMaxIters>();
}

#endif /* HAVE_CTRILINOS_AZTECOO */


//
// Definitions for AztecOO_StatusTestResNorm
//

#ifdef HAVE_CTRILINOS_AZTECOO

/* get AztecOO_StatusTestResNorm from non-const table using CT_AztecOO_StatusTestResNorm_ID */
const Teuchos::RCP<AztecOO_StatusTestResNorm>
CAztecOO::getStatusTestResNorm( CT_AztecOO_StatusTestResNorm_ID_t id )
{
    return CTrilinos::tableRepos().get<AztecOO_StatusTestResNorm>(
        CTrilinos::abstractType<CT_AztecOO_StatusTestResNorm_ID_t>(id));
}

/* get AztecOO_StatusTestResNorm from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<AztecOO_StatusTestResNorm>
CAztecOO::getStatusTestResNorm( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().get<AztecOO_StatusTestResNorm>(id);
}

/* get const AztecOO_StatusTestResNorm from either the const or non-const table
 * using CT_AztecOO_StatusTestResNorm_ID */
const Teuchos::RCP<const AztecOO_StatusTestResNorm>
CAztecOO::getConstStatusTestResNorm( CT_AztecOO_StatusTestResNorm_ID_t id )
{
    return CTrilinos::tableRepos().getConst<AztecOO_StatusTestResNorm>(
        CTrilinos::abstractType<CT_AztecOO_StatusTestResNorm_ID_t>(id));
}

/* get const AztecOO_StatusTestResNorm from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const AztecOO_StatusTestResNorm>
CAztecOO::getConstStatusTestResNorm( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().getConst<AztecOO_StatusTestResNorm>(id);
}

/* store AztecOO_StatusTestResNorm (owned) in non-const table */
CT_AztecOO_StatusTestResNorm_ID_t
CAztecOO::storeNewStatusTestResNorm( AztecOO_StatusTestResNorm *pobj )
{
    return CTrilinos::concreteType<CT_AztecOO_StatusTestResNorm_ID_t>(
        CTrilinos::tableRepos().store<AztecOO_StatusTestResNorm>(pobj, true));
}

/* store AztecOO_StatusTestResNorm in non-const table */
CT_AztecOO_StatusTestResNorm_ID_t
CAztecOO::storeStatusTestResNorm( AztecOO_StatusTestResNorm *pobj )
{
    return CTrilinos::concreteType<CT_AztecOO_StatusTestResNorm_ID_t>(
        CTrilinos::tableRepos().store<AztecOO_StatusTestResNorm>(pobj, false));
}

/* store const AztecOO_StatusTestResNorm in const table */
CT_AztecOO_StatusTestResNorm_ID_t
CAztecOO::storeConstStatusTestResNorm( const AztecOO_StatusTestResNorm *pobj )
{
    return CTrilinos::concreteType<CT_AztecOO_StatusTestResNorm_ID_t>(
        CTrilinos::tableRepos().store<AztecOO_StatusTestResNorm>(pobj, false));
}

/* remove AztecOO_StatusTestResNorm from table using CT_AztecOO_StatusTestResNorm_ID */
void
CAztecOO::removeStatusTestResNorm( CT_AztecOO_StatusTestResNorm_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_AztecOO_StatusTestResNorm_ID_t>(*id);
    CTrilinos::tableRepos().remove(&aid);
    *id = CTrilinos::concreteType<CT_AztecOO_StatusTestResNorm_ID_t>(aid);
}

/* purge AztecOO_StatusTestResNorm table */
void
CAztecOO::purgeStatusTestResNorm(  )
{
    CTrilinos::tableRepos().purge<AztecOO_StatusTestResNorm>();
}

#endif /* HAVE_CTRILINOS_AZTECOO */


//
// Definitions for Ifpack_Preconditioner
//

#ifdef HAVE_CTRILINOS_IFPACK

/* get Ifpack_Preconditioner from non-const table using CT_Ifpack_Preconditioner_ID */
const Teuchos::RCP<Ifpack_Preconditioner>
CIfpack::getPreconditioner( CT_Ifpack_Preconditioner_ID_t id )
{
    return CTrilinos::tableRepos().get<Ifpack_Preconditioner>(
        CTrilinos::abstractType<CT_Ifpack_Preconditioner_ID_t>(id));
}

/* get Ifpack_Preconditioner from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Ifpack_Preconditioner>
CIfpack::getPreconditioner( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().get<Ifpack_Preconditioner>(id);
}

/* get const Ifpack_Preconditioner from either the const or non-const table
 * using CT_Ifpack_Preconditioner_ID */
const Teuchos::RCP<const Ifpack_Preconditioner>
CIfpack::getConstPreconditioner( CT_Ifpack_Preconditioner_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Ifpack_Preconditioner>(
        CTrilinos::abstractType<CT_Ifpack_Preconditioner_ID_t>(id));
}

/* get const Ifpack_Preconditioner from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Ifpack_Preconditioner>
CIfpack::getConstPreconditioner( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::tableRepos().getConst<Ifpack_Preconditioner>(id);
}

/* store Ifpack_Preconditioner (owned) in non-const table */
CT_Ifpack_Preconditioner_ID_t
CIfpack::storeNewPreconditioner( Ifpack_Preconditioner *pobj )
{
    return CTrilinos::concreteType<CT_Ifpack_Preconditioner_ID_t>(
        CTrilinos::tableRepos().store<Ifpack_Preconditioner>(pobj, true));
}

/* store Ifpack_Preconditioner in non-const table */
CT_Ifpack_Preconditioner_ID_t
CIfpack::storePreconditioner( Ifpack_Preconditioner *pobj )
{
    return CTrilinos::concreteType<CT_Ifpack_Preconditioner_ID_t>(
        CTrilinos::tableRepos().store<Ifpack_Preconditioner>(pobj, false));
}

/* store const Ifpack_Preconditioner in const table */
CT_Ifpack_Preconditioner_ID_t
CIfpack::storeConstPreconditioner( const Ifpack_Preconditioner *pobj )
{
    return CTrilinos::concreteType<CT_Ifpack_Preconditioner_ID_t>(
        CTrilinos::tableRepos().store<Ifpack_Preconditioner>(pobj, false));
}

/* remove Ifpack_Preconditioner from table using CT_Ifpack_Preconditioner_ID */
void
CIfpack::removePreconditioner( CT_Ifpack_Preconditioner_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Ifpack_Preconditioner_ID_t>(*id);
    CTrilinos::tableRepos().remove(&aid);
    *id = CTrilinos::concreteType<CT_Ifpack_Preconditioner_ID_t>(aid);
}

/* purge Ifpack_Preconditioner table */
void
CIfpack::purgePreconditioner(  )
{
    CTrilinos::tableRepos().purge<Ifpack_Preconditioner>();
}

#endif /* HAVE_CTRILINOS_IFPACK */


