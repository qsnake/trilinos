
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


/*! @file CTrilinos_TableRepos.hpp
 * @brief Central table repository for CTrilinos. */


#ifndef CTRILINOS_TABLEREPOS_HPP
#define CTRILINOS_TABLEREPOS_HPP


#include "CTrilinos_config.h"


#include <string>

#include "Epetra_Distributor.h"
#include "Epetra_SerialComm.h"
#include "Epetra_BLAS.h"
#include "Epetra_Comm.h"
#include "Epetra_Operator.h"
#include "Epetra_MultiVector.h"
#include "Epetra_OffsetIndex.h"
#include "Epetra_Object.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CompObject.h"
#include "Epetra_Directory.h"
#include "Epetra_Flops.h"
#include "Epetra_SrcDistObject.h"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#endif /* HAVE_MPI */
#include "Epetra_CrsMatrix.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_DistObject.h"
#include "Epetra_Vector.h"
#include "Epetra_Export.h"
#include "Epetra_Map.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Import.h"
#include "Epetra_Time.h"
#include "Epetra_JadMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_LAPACK.h"
#include "Teuchos_ParameterList.hpp"
#ifdef HAVE_CTRILINOS_AMESOS
#include "Amesos_BaseSolver.h"
#endif /* HAVE_CTRILINOS_AMESOS */
#include "Epetra_FECrsMatrix.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_SerialDenseMatrix.h"
#ifdef HAVE_CTRILINOS_AZTECOO
#include "AztecOO_StatusTest.h"
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
#include "AztecOO_StatusTestCombo.h"
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
#include "AztecOO_StatusTestMaxIters.h"
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
#include "AztecOO_StatusTestResNorm.h"
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_IFPACK
#include "Ifpack_Preconditioner.h"
#endif /* HAVE_CTRILINOS_IFPACK */

#include "CTrilinos_enums.h"
#include "CTrilinos_exceptions.hpp"
#include "CTrilinos_Table.hpp"
#include "Teuchos_RCP.hpp"


namespace CTrilinos {


class TableRepos
{
  public:

    /*! constructor */
    TableRepos();

    /*! destructor */
    ~TableRepos();

    /*! retrieve the object */
    template <class T>
    const Teuchos::RCP<T> get(CTrilinos_Universal_ID_t id);

    /*! retrieve the object */
    template <class T>
    const Teuchos::RCP<const T> getConst(CTrilinos_Universal_ID_t id);

    /*! store a non-const RCP to object of type T */
    template <class T>
    CTrilinos_Universal_ID_t store(T* pobj, bool owned);

    /*! store a const RCP to object of type T */
    template <class T>
    CTrilinos_Universal_ID_t store(const T* pobj, bool owned);

    /*! remove an object from the table and invalidate the id struct */
    void remove(CTrilinos_Universal_ID_t * id);

    /*! dump a specific table's content but keep its properties */
    template <class T>
    void purge();

    /*! dump the tables' content but keep their properties */
    void purgeAll();

    /*! get the CTrilinos::Table for the given type */
    void getTable(Table<Epetra_Distributor> *& tab);
    void getTable(Table<Epetra_SerialComm> *& tab);
    void getTable(Table<Epetra_BLAS> *& tab);
    void getTable(Table<Epetra_Comm> *& tab);
    void getTable(Table<Epetra_Operator> *& tab);
    void getTable(Table<Epetra_MultiVector> *& tab);
    void getTable(Table<Epetra_OffsetIndex> *& tab);
    void getTable(Table<Epetra_Object> *& tab);
    void getTable(Table<Epetra_RowMatrix> *& tab);
    void getTable(Table<Epetra_CompObject> *& tab);
    void getTable(Table<Epetra_Directory> *& tab);
    void getTable(Table<Epetra_Flops> *& tab);
    void getTable(Table<Epetra_SrcDistObject> *& tab);
#ifdef HAVE_MPI
    void getTable(Table<Epetra_MpiComm> *& tab);
#endif /* HAVE_MPI */
    void getTable(Table<Epetra_CrsMatrix> *& tab);
    void getTable(Table<Epetra_CrsGraph> *& tab);
    void getTable(Table<Epetra_DistObject> *& tab);
    void getTable(Table<Epetra_Vector> *& tab);
    void getTable(Table<Epetra_Export> *& tab);
    void getTable(Table<Epetra_Map> *& tab);
    void getTable(Table<Epetra_BlockMap> *& tab);
    void getTable(Table<Epetra_Import> *& tab);
    void getTable(Table<Epetra_Time> *& tab);
    void getTable(Table<Epetra_JadMatrix> *& tab);
    void getTable(Table<Epetra_LinearProblem> *& tab);
    void getTable(Table<Epetra_LAPACK> *& tab);
    void getTable(Table<Teuchos::ParameterList> *& tab);
#ifdef HAVE_CTRILINOS_AMESOS
    void getTable(Table<Amesos_BaseSolver> *& tab);
#endif /* HAVE_CTRILINOS_AMESOS */
    void getTable(Table<Epetra_FECrsMatrix> *& tab);
    void getTable(Table<Epetra_IntSerialDenseVector> *& tab);
    void getTable(Table<Epetra_SerialDenseMatrix> *& tab);
#ifdef HAVE_CTRILINOS_AZTECOO
    void getTable(Table<AztecOO_StatusTest> *& tab);
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
    void getTable(Table<AztecOO_StatusTestCombo> *& tab);
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
    void getTable(Table<AztecOO_StatusTestMaxIters> *& tab);
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
    void getTable(Table<AztecOO_StatusTestResNorm> *& tab);
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_IFPACK
    void getTable(Table<Ifpack_Preconditioner> *& tab);
#endif /* HAVE_CTRILINOS_IFPACK */

    /*! create an alias for the object in another table */
    CTrilinos_Universal_ID_t alias(CTrilinos_Universal_ID_t id, CTrilinos_Table_ID_t tab, bool keepold = true);

    /*! create an alias for the object in another table */
    template <class T>
    CTrilinos_Universal_ID_t do_alias(Table<T> &tab, CTrilinos_Universal_ID_t &aid, bool keepold);

    /*! create an alias for the object in another table */
    template <class T>
    CTrilinos_Universal_ID_t do_alias_const(Table<T> &tab, CTrilinos_Universal_ID_t &aid, bool keepold);

    /*! retrieve the object */
    template <class T>
    const Teuchos::RCP<T> getPoly(CTrilinos_Universal_ID_t aid);

    /*! retrieve the object */
    template <class T>
    const Teuchos::RCP<const T> getConstPoly(CTrilinos_Universal_ID_t aid);

    /*! see if the object is dynamic_cast'able */
    bool typeCheck(CTrilinos_Universal_ID_t aid, CTrilinos_Table_ID_t type);

    Table<Epetra_Distributor> tab_Epetra_Distributor;
    Table<Epetra_SerialComm> tab_Epetra_SerialComm;
    Table<Epetra_BLAS> tab_Epetra_BLAS;
    Table<Epetra_Comm> tab_Epetra_Comm;
    Table<Epetra_Operator> tab_Epetra_Operator;
    Table<Epetra_MultiVector> tab_Epetra_MultiVector;
    Table<Epetra_OffsetIndex> tab_Epetra_OffsetIndex;
    Table<Epetra_Object> tab_Epetra_Object;
    Table<Epetra_RowMatrix> tab_Epetra_RowMatrix;
    Table<Epetra_CompObject> tab_Epetra_CompObject;
    Table<Epetra_Directory> tab_Epetra_Directory;
    Table<Epetra_Flops> tab_Epetra_Flops;
    Table<Epetra_SrcDistObject> tab_Epetra_SrcDistObject;
#ifdef HAVE_MPI
    Table<Epetra_MpiComm> tab_Epetra_MpiComm;
#endif /* HAVE_MPI */
    Table<Epetra_CrsMatrix> tab_Epetra_CrsMatrix;
    Table<Epetra_CrsGraph> tab_Epetra_CrsGraph;
    Table<Epetra_DistObject> tab_Epetra_DistObject;
    Table<Epetra_Vector> tab_Epetra_Vector;
    Table<Epetra_Export> tab_Epetra_Export;
    Table<Epetra_Map> tab_Epetra_Map;
    Table<Epetra_BlockMap> tab_Epetra_BlockMap;
    Table<Epetra_Import> tab_Epetra_Import;
    Table<Epetra_Time> tab_Epetra_Time;
    Table<Epetra_JadMatrix> tab_Epetra_JadMatrix;
    Table<Epetra_LinearProblem> tab_Epetra_LinearProblem;
    Table<Epetra_LAPACK> tab_Epetra_LAPACK;
    Table<Teuchos::ParameterList> tab_Teuchos_ParameterList;
#ifdef HAVE_CTRILINOS_AMESOS
    Table<Amesos_BaseSolver> tab_Amesos_BaseSolver;
#endif /* HAVE_CTRILINOS_AMESOS */
    Table<Epetra_FECrsMatrix> tab_Epetra_FECrsMatrix;
    Table<Epetra_IntSerialDenseVector> tab_Epetra_IntSerialDenseVector;
    Table<Epetra_SerialDenseMatrix> tab_Epetra_SerialDenseMatrix;
#ifdef HAVE_CTRILINOS_AZTECOO
    Table<AztecOO_StatusTest> tab_AztecOO_StatusTest;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
    Table<AztecOO_StatusTestCombo> tab_AztecOO_StatusTestCombo;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
    Table<AztecOO_StatusTestMaxIters> tab_AztecOO_StatusTestMaxIters;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
    Table<AztecOO_StatusTestResNorm> tab_AztecOO_StatusTestResNorm;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_IFPACK
    Table<Ifpack_Preconditioner> tab_Ifpack_Preconditioner;
#endif /* HAVE_CTRILINOS_IFPACK */

    bool call_me_lazy;  /* I was too lazy to deal with the commas in the init list, so... */
};

TableRepos & tableRepos(  );


template <class T>
const Teuchos::RCP<T> TableRepos::get(CTrilinos_Universal_ID_t aid)
{
    /* Shortcut if stored in the most obvious table */
    Table<T> *tab_maybe = 0;
    getTable(tab_maybe);
    if (tab_maybe->isType(aid.table))
        return tab_maybe->get<T>(aid);
    else
        return getPoly<T>(aid);
}

template <class T>
const Teuchos::RCP<T> TableRepos::getPoly(CTrilinos_Universal_ID_t aid)
{
    switch (aid.table) {
    case CT_Epetra_Distributor_ID:
        return tab_Epetra_Distributor.get<T>(aid);
        break;
    case CT_Epetra_SerialComm_ID:
        return tab_Epetra_SerialComm.get<T>(aid);
        break;
    case CT_Epetra_BLAS_ID:
        return tab_Epetra_BLAS.get<T>(aid);
        break;
    case CT_Epetra_Comm_ID:
        return tab_Epetra_Comm.get<T>(aid);
        break;
    case CT_Epetra_Operator_ID:
        return tab_Epetra_Operator.get<T>(aid);
        break;
    case CT_Epetra_MultiVector_ID:
        return tab_Epetra_MultiVector.get<T>(aid);
        break;
    case CT_Epetra_OffsetIndex_ID:
        return tab_Epetra_OffsetIndex.get<T>(aid);
        break;
    case CT_Epetra_Object_ID:
        return tab_Epetra_Object.get<T>(aid);
        break;
    case CT_Epetra_RowMatrix_ID:
        return tab_Epetra_RowMatrix.get<T>(aid);
        break;
    case CT_Epetra_CompObject_ID:
        return tab_Epetra_CompObject.get<T>(aid);
        break;
    case CT_Epetra_Directory_ID:
        return tab_Epetra_Directory.get<T>(aid);
        break;
    case CT_Epetra_Flops_ID:
        return tab_Epetra_Flops.get<T>(aid);
        break;
    case CT_Epetra_SrcDistObject_ID:
        return tab_Epetra_SrcDistObject.get<T>(aid);
        break;
#ifdef HAVE_MPI
    case CT_Epetra_MpiComm_ID:
        return tab_Epetra_MpiComm.get<T>(aid);
        break;
#endif /* HAVE_MPI */
    case CT_Epetra_CrsMatrix_ID:
        return tab_Epetra_CrsMatrix.get<T>(aid);
        break;
    case CT_Epetra_CrsGraph_ID:
        return tab_Epetra_CrsGraph.get<T>(aid);
        break;
    case CT_Epetra_DistObject_ID:
        return tab_Epetra_DistObject.get<T>(aid);
        break;
    case CT_Epetra_Vector_ID:
        return tab_Epetra_Vector.get<T>(aid);
        break;
    case CT_Epetra_Export_ID:
        return tab_Epetra_Export.get<T>(aid);
        break;
    case CT_Epetra_Map_ID:
        return tab_Epetra_Map.get<T>(aid);
        break;
    case CT_Epetra_BlockMap_ID:
        return tab_Epetra_BlockMap.get<T>(aid);
        break;
    case CT_Epetra_Import_ID:
        return tab_Epetra_Import.get<T>(aid);
        break;
    case CT_Epetra_Time_ID:
        return tab_Epetra_Time.get<T>(aid);
        break;
    case CT_Epetra_JadMatrix_ID:
        return tab_Epetra_JadMatrix.get<T>(aid);
        break;
    case CT_Epetra_LinearProblem_ID:
        return tab_Epetra_LinearProblem.get<T>(aid);
        break;
    case CT_Epetra_LAPACK_ID:
        return tab_Epetra_LAPACK.get<T>(aid);
        break;
    case CT_Teuchos_ParameterList_ID:
        return tab_Teuchos_ParameterList.get<T>(aid);
        break;
#ifdef HAVE_CTRILINOS_AMESOS
    case CT_Amesos_BaseSolver_ID:
        return tab_Amesos_BaseSolver.get<T>(aid);
        break;
#endif /* HAVE_CTRILINOS_AMESOS */
    case CT_Epetra_FECrsMatrix_ID:
        return tab_Epetra_FECrsMatrix.get<T>(aid);
        break;
    case CT_Epetra_IntSerialDenseVector_ID:
        return tab_Epetra_IntSerialDenseVector.get<T>(aid);
        break;
    case CT_Epetra_SerialDenseMatrix_ID:
        return tab_Epetra_SerialDenseMatrix.get<T>(aid);
        break;
#ifdef HAVE_CTRILINOS_AZTECOO
    case CT_AztecOO_StatusTest_ID:
        return tab_AztecOO_StatusTest.get<T>(aid);
        break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
    case CT_AztecOO_StatusTestCombo_ID:
        return tab_AztecOO_StatusTestCombo.get<T>(aid);
        break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
    case CT_AztecOO_StatusTestMaxIters_ID:
        return tab_AztecOO_StatusTestMaxIters.get<T>(aid);
        break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
    case CT_AztecOO_StatusTestResNorm_ID:
        return tab_AztecOO_StatusTestResNorm.get<T>(aid);
        break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_IFPACK
    case CT_Ifpack_Preconditioner_ID:
        return tab_Ifpack_Preconditioner.get<T>(aid);
        break;
#endif /* HAVE_CTRILINOS_IFPACK */
    default:
        throw CTrilinosInvalidTypeError("invalid table id");
    }

    return Teuchos::null;
}

template <class T>
const Teuchos::RCP<const T> TableRepos::getConst(CTrilinos_Universal_ID_t aid)
{
    /* Shortcut if stored in the most obvious table */
    Table<T> *tab_maybe = 0;
    getTable(tab_maybe);
    if (tab_maybe->isType(aid.table))
        return tab_maybe->getConst<T>(aid);
    else
        return getConstPoly<T>(aid);
}

template <class T>
const Teuchos::RCP<const T> TableRepos::getConstPoly(CTrilinos_Universal_ID_t aid)
{
    switch (aid.table) {
    case CT_Epetra_Distributor_ID:
        return tab_Epetra_Distributor.getConst<T>(aid);
        break;
    case CT_Epetra_SerialComm_ID:
        return tab_Epetra_SerialComm.getConst<T>(aid);
        break;
    case CT_Epetra_BLAS_ID:
        return tab_Epetra_BLAS.getConst<T>(aid);
        break;
    case CT_Epetra_Comm_ID:
        return tab_Epetra_Comm.getConst<T>(aid);
        break;
    case CT_Epetra_Operator_ID:
        return tab_Epetra_Operator.getConst<T>(aid);
        break;
    case CT_Epetra_MultiVector_ID:
        return tab_Epetra_MultiVector.getConst<T>(aid);
        break;
    case CT_Epetra_OffsetIndex_ID:
        return tab_Epetra_OffsetIndex.getConst<T>(aid);
        break;
    case CT_Epetra_Object_ID:
        return tab_Epetra_Object.getConst<T>(aid);
        break;
    case CT_Epetra_RowMatrix_ID:
        return tab_Epetra_RowMatrix.getConst<T>(aid);
        break;
    case CT_Epetra_CompObject_ID:
        return tab_Epetra_CompObject.getConst<T>(aid);
        break;
    case CT_Epetra_Directory_ID:
        return tab_Epetra_Directory.getConst<T>(aid);
        break;
    case CT_Epetra_Flops_ID:
        return tab_Epetra_Flops.getConst<T>(aid);
        break;
    case CT_Epetra_SrcDistObject_ID:
        return tab_Epetra_SrcDistObject.getConst<T>(aid);
        break;
#ifdef HAVE_MPI
    case CT_Epetra_MpiComm_ID:
        return tab_Epetra_MpiComm.getConst<T>(aid);
        break;
#endif /* HAVE_MPI */
    case CT_Epetra_CrsMatrix_ID:
        return tab_Epetra_CrsMatrix.getConst<T>(aid);
        break;
    case CT_Epetra_CrsGraph_ID:
        return tab_Epetra_CrsGraph.getConst<T>(aid);
        break;
    case CT_Epetra_DistObject_ID:
        return tab_Epetra_DistObject.getConst<T>(aid);
        break;
    case CT_Epetra_Vector_ID:
        return tab_Epetra_Vector.getConst<T>(aid);
        break;
    case CT_Epetra_Export_ID:
        return tab_Epetra_Export.getConst<T>(aid);
        break;
    case CT_Epetra_Map_ID:
        return tab_Epetra_Map.getConst<T>(aid);
        break;
    case CT_Epetra_BlockMap_ID:
        return tab_Epetra_BlockMap.getConst<T>(aid);
        break;
    case CT_Epetra_Import_ID:
        return tab_Epetra_Import.getConst<T>(aid);
        break;
    case CT_Epetra_Time_ID:
        return tab_Epetra_Time.getConst<T>(aid);
        break;
    case CT_Epetra_JadMatrix_ID:
        return tab_Epetra_JadMatrix.getConst<T>(aid);
        break;
    case CT_Epetra_LinearProblem_ID:
        return tab_Epetra_LinearProblem.getConst<T>(aid);
        break;
    case CT_Epetra_LAPACK_ID:
        return tab_Epetra_LAPACK.getConst<T>(aid);
        break;
    case CT_Teuchos_ParameterList_ID:
        return tab_Teuchos_ParameterList.getConst<T>(aid);
        break;
#ifdef HAVE_CTRILINOS_AMESOS
    case CT_Amesos_BaseSolver_ID:
        return tab_Amesos_BaseSolver.getConst<T>(aid);
        break;
#endif /* HAVE_CTRILINOS_AMESOS */
    case CT_Epetra_FECrsMatrix_ID:
        return tab_Epetra_FECrsMatrix.getConst<T>(aid);
        break;
    case CT_Epetra_IntSerialDenseVector_ID:
        return tab_Epetra_IntSerialDenseVector.getConst<T>(aid);
        break;
    case CT_Epetra_SerialDenseMatrix_ID:
        return tab_Epetra_SerialDenseMatrix.getConst<T>(aid);
        break;
#ifdef HAVE_CTRILINOS_AZTECOO
    case CT_AztecOO_StatusTest_ID:
        return tab_AztecOO_StatusTest.getConst<T>(aid);
        break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
    case CT_AztecOO_StatusTestCombo_ID:
        return tab_AztecOO_StatusTestCombo.getConst<T>(aid);
        break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
    case CT_AztecOO_StatusTestMaxIters_ID:
        return tab_AztecOO_StatusTestMaxIters.getConst<T>(aid);
        break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
    case CT_AztecOO_StatusTestResNorm_ID:
        return tab_AztecOO_StatusTestResNorm.getConst<T>(aid);
        break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_IFPACK
    case CT_Ifpack_Preconditioner_ID:
        return tab_Ifpack_Preconditioner.getConst<T>(aid);
        break;
#endif /* HAVE_CTRILINOS_IFPACK */
    default:
        throw CTrilinosInvalidTypeError("invalid table id");
    }

    return Teuchos::null;
}

template <class T>
CTrilinos_Universal_ID_t TableRepos::do_alias(
    Table<T> &tab, CTrilinos_Universal_ID_t &aid, bool keepold)
{
    CTrilinos_Universal_ID_t newid;

    switch (aid.table) {
    case CT_Epetra_Distributor_ID:
        newid = tab.alias(tab_Epetra_Distributor.get<T>(aid));
        if (!keepold) tab_Epetra_Distributor.remove(&aid);
        break;
    case CT_Epetra_SerialComm_ID:
        newid = tab.alias(tab_Epetra_SerialComm.get<T>(aid));
        if (!keepold) tab_Epetra_SerialComm.remove(&aid);
        break;
    case CT_Epetra_BLAS_ID:
        newid = tab.alias(tab_Epetra_BLAS.get<T>(aid));
        if (!keepold) tab_Epetra_BLAS.remove(&aid);
        break;
    case CT_Epetra_Comm_ID:
        newid = tab.alias(tab_Epetra_Comm.get<T>(aid));
        if (!keepold) tab_Epetra_Comm.remove(&aid);
        break;
    case CT_Epetra_Operator_ID:
        newid = tab.alias(tab_Epetra_Operator.get<T>(aid));
        if (!keepold) tab_Epetra_Operator.remove(&aid);
        break;
    case CT_Epetra_MultiVector_ID:
        newid = tab.alias(tab_Epetra_MultiVector.get<T>(aid));
        if (!keepold) tab_Epetra_MultiVector.remove(&aid);
        break;
    case CT_Epetra_OffsetIndex_ID:
        newid = tab.alias(tab_Epetra_OffsetIndex.get<T>(aid));
        if (!keepold) tab_Epetra_OffsetIndex.remove(&aid);
        break;
    case CT_Epetra_Object_ID:
        newid = tab.alias(tab_Epetra_Object.get<T>(aid));
        if (!keepold) tab_Epetra_Object.remove(&aid);
        break;
    case CT_Epetra_RowMatrix_ID:
        newid = tab.alias(tab_Epetra_RowMatrix.get<T>(aid));
        if (!keepold) tab_Epetra_RowMatrix.remove(&aid);
        break;
    case CT_Epetra_CompObject_ID:
        newid = tab.alias(tab_Epetra_CompObject.get<T>(aid));
        if (!keepold) tab_Epetra_CompObject.remove(&aid);
        break;
    case CT_Epetra_Directory_ID:
        newid = tab.alias(tab_Epetra_Directory.get<T>(aid));
        if (!keepold) tab_Epetra_Directory.remove(&aid);
        break;
    case CT_Epetra_Flops_ID:
        newid = tab.alias(tab_Epetra_Flops.get<T>(aid));
        if (!keepold) tab_Epetra_Flops.remove(&aid);
        break;
    case CT_Epetra_SrcDistObject_ID:
        newid = tab.alias(tab_Epetra_SrcDistObject.get<T>(aid));
        if (!keepold) tab_Epetra_SrcDistObject.remove(&aid);
        break;
#ifdef HAVE_MPI
    case CT_Epetra_MpiComm_ID:
        newid = tab.alias(tab_Epetra_MpiComm.get<T>(aid));
        if (!keepold) tab_Epetra_MpiComm.remove(&aid);
        break;
#endif /* HAVE_MPI */
    case CT_Epetra_CrsMatrix_ID:
        newid = tab.alias(tab_Epetra_CrsMatrix.get<T>(aid));
        if (!keepold) tab_Epetra_CrsMatrix.remove(&aid);
        break;
    case CT_Epetra_CrsGraph_ID:
        newid = tab.alias(tab_Epetra_CrsGraph.get<T>(aid));
        if (!keepold) tab_Epetra_CrsGraph.remove(&aid);
        break;
    case CT_Epetra_DistObject_ID:
        newid = tab.alias(tab_Epetra_DistObject.get<T>(aid));
        if (!keepold) tab_Epetra_DistObject.remove(&aid);
        break;
    case CT_Epetra_Vector_ID:
        newid = tab.alias(tab_Epetra_Vector.get<T>(aid));
        if (!keepold) tab_Epetra_Vector.remove(&aid);
        break;
    case CT_Epetra_Export_ID:
        newid = tab.alias(tab_Epetra_Export.get<T>(aid));
        if (!keepold) tab_Epetra_Export.remove(&aid);
        break;
    case CT_Epetra_Map_ID:
        newid = tab.alias(tab_Epetra_Map.get<T>(aid));
        if (!keepold) tab_Epetra_Map.remove(&aid);
        break;
    case CT_Epetra_BlockMap_ID:
        newid = tab.alias(tab_Epetra_BlockMap.get<T>(aid));
        if (!keepold) tab_Epetra_BlockMap.remove(&aid);
        break;
    case CT_Epetra_Import_ID:
        newid = tab.alias(tab_Epetra_Import.get<T>(aid));
        if (!keepold) tab_Epetra_Import.remove(&aid);
        break;
    case CT_Epetra_Time_ID:
        newid = tab.alias(tab_Epetra_Time.get<T>(aid));
        if (!keepold) tab_Epetra_Time.remove(&aid);
        break;
    case CT_Epetra_JadMatrix_ID:
        newid = tab.alias(tab_Epetra_JadMatrix.get<T>(aid));
        if (!keepold) tab_Epetra_JadMatrix.remove(&aid);
        break;
    case CT_Epetra_LinearProblem_ID:
        newid = tab.alias(tab_Epetra_LinearProblem.get<T>(aid));
        if (!keepold) tab_Epetra_LinearProblem.remove(&aid);
        break;
    case CT_Epetra_LAPACK_ID:
        newid = tab.alias(tab_Epetra_LAPACK.get<T>(aid));
        if (!keepold) tab_Epetra_LAPACK.remove(&aid);
        break;
    case CT_Teuchos_ParameterList_ID:
        newid = tab.alias(tab_Teuchos_ParameterList.get<T>(aid));
        if (!keepold) tab_Teuchos_ParameterList.remove(&aid);
        break;
#ifdef HAVE_CTRILINOS_AMESOS
    case CT_Amesos_BaseSolver_ID:
        newid = tab.alias(tab_Amesos_BaseSolver.get<T>(aid));
        if (!keepold) tab_Amesos_BaseSolver.remove(&aid);
        break;
#endif /* HAVE_CTRILINOS_AMESOS */
    case CT_Epetra_FECrsMatrix_ID:
        newid = tab.alias(tab_Epetra_FECrsMatrix.get<T>(aid));
        if (!keepold) tab_Epetra_FECrsMatrix.remove(&aid);
        break;
    case CT_Epetra_IntSerialDenseVector_ID:
        newid = tab.alias(tab_Epetra_IntSerialDenseVector.get<T>(aid));
        if (!keepold) tab_Epetra_IntSerialDenseVector.remove(&aid);
        break;
    case CT_Epetra_SerialDenseMatrix_ID:
        newid = tab.alias(tab_Epetra_SerialDenseMatrix.get<T>(aid));
        if (!keepold) tab_Epetra_SerialDenseMatrix.remove(&aid);
        break;
#ifdef HAVE_CTRILINOS_AZTECOO
    case CT_AztecOO_StatusTest_ID:
        newid = tab.alias(tab_AztecOO_StatusTest.get<T>(aid));
        if (!keepold) tab_AztecOO_StatusTest.remove(&aid);
        break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
    case CT_AztecOO_StatusTestCombo_ID:
        newid = tab.alias(tab_AztecOO_StatusTestCombo.get<T>(aid));
        if (!keepold) tab_AztecOO_StatusTestCombo.remove(&aid);
        break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
    case CT_AztecOO_StatusTestMaxIters_ID:
        newid = tab.alias(tab_AztecOO_StatusTestMaxIters.get<T>(aid));
        if (!keepold) tab_AztecOO_StatusTestMaxIters.remove(&aid);
        break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
    case CT_AztecOO_StatusTestResNorm_ID:
        newid = tab.alias(tab_AztecOO_StatusTestResNorm.get<T>(aid));
        if (!keepold) tab_AztecOO_StatusTestResNorm.remove(&aid);
        break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_IFPACK
    case CT_Ifpack_Preconditioner_ID:
        newid = tab.alias(tab_Ifpack_Preconditioner.get<T>(aid));
        if (!keepold) tab_Ifpack_Preconditioner.remove(&aid);
        break;
#endif /* HAVE_CTRILINOS_IFPACK */
    default:
        throw CTrilinosInvalidTypeError("invalid table id or non-polymorphic class");
    }

    return newid;
}

template <class T>
CTrilinos_Universal_ID_t TableRepos::do_alias_const(
    Table<T> &tab, CTrilinos_Universal_ID_t &aid, bool keepold)
{
    CTrilinos_Universal_ID_t newid;

    switch (aid.table) {
    case CT_Epetra_Distributor_ID:
        newid = tab.alias(tab_Epetra_Distributor.getConst<T>(aid));
        if (!keepold) tab_Epetra_Distributor.remove(&aid);
        break;
    case CT_Epetra_SerialComm_ID:
        newid = tab.alias(tab_Epetra_SerialComm.getConst<T>(aid));
        if (!keepold) tab_Epetra_SerialComm.remove(&aid);
        break;
    case CT_Epetra_BLAS_ID:
        newid = tab.alias(tab_Epetra_BLAS.getConst<T>(aid));
        if (!keepold) tab_Epetra_BLAS.remove(&aid);
        break;
    case CT_Epetra_Comm_ID:
        newid = tab.alias(tab_Epetra_Comm.getConst<T>(aid));
        if (!keepold) tab_Epetra_Comm.remove(&aid);
        break;
    case CT_Epetra_Operator_ID:
        newid = tab.alias(tab_Epetra_Operator.getConst<T>(aid));
        if (!keepold) tab_Epetra_Operator.remove(&aid);
        break;
    case CT_Epetra_MultiVector_ID:
        newid = tab.alias(tab_Epetra_MultiVector.getConst<T>(aid));
        if (!keepold) tab_Epetra_MultiVector.remove(&aid);
        break;
    case CT_Epetra_OffsetIndex_ID:
        newid = tab.alias(tab_Epetra_OffsetIndex.getConst<T>(aid));
        if (!keepold) tab_Epetra_OffsetIndex.remove(&aid);
        break;
    case CT_Epetra_Object_ID:
        newid = tab.alias(tab_Epetra_Object.getConst<T>(aid));
        if (!keepold) tab_Epetra_Object.remove(&aid);
        break;
    case CT_Epetra_RowMatrix_ID:
        newid = tab.alias(tab_Epetra_RowMatrix.getConst<T>(aid));
        if (!keepold) tab_Epetra_RowMatrix.remove(&aid);
        break;
    case CT_Epetra_CompObject_ID:
        newid = tab.alias(tab_Epetra_CompObject.getConst<T>(aid));
        if (!keepold) tab_Epetra_CompObject.remove(&aid);
        break;
    case CT_Epetra_Directory_ID:
        newid = tab.alias(tab_Epetra_Directory.getConst<T>(aid));
        if (!keepold) tab_Epetra_Directory.remove(&aid);
        break;
    case CT_Epetra_Flops_ID:
        newid = tab.alias(tab_Epetra_Flops.getConst<T>(aid));
        if (!keepold) tab_Epetra_Flops.remove(&aid);
        break;
    case CT_Epetra_SrcDistObject_ID:
        newid = tab.alias(tab_Epetra_SrcDistObject.getConst<T>(aid));
        if (!keepold) tab_Epetra_SrcDistObject.remove(&aid);
        break;
#ifdef HAVE_MPI
    case CT_Epetra_MpiComm_ID:
        newid = tab.alias(tab_Epetra_MpiComm.getConst<T>(aid));
        if (!keepold) tab_Epetra_MpiComm.remove(&aid);
        break;
#endif /* HAVE_MPI */
    case CT_Epetra_CrsMatrix_ID:
        newid = tab.alias(tab_Epetra_CrsMatrix.getConst<T>(aid));
        if (!keepold) tab_Epetra_CrsMatrix.remove(&aid);
        break;
    case CT_Epetra_CrsGraph_ID:
        newid = tab.alias(tab_Epetra_CrsGraph.getConst<T>(aid));
        if (!keepold) tab_Epetra_CrsGraph.remove(&aid);
        break;
    case CT_Epetra_DistObject_ID:
        newid = tab.alias(tab_Epetra_DistObject.getConst<T>(aid));
        if (!keepold) tab_Epetra_DistObject.remove(&aid);
        break;
    case CT_Epetra_Vector_ID:
        newid = tab.alias(tab_Epetra_Vector.getConst<T>(aid));
        if (!keepold) tab_Epetra_Vector.remove(&aid);
        break;
    case CT_Epetra_Export_ID:
        newid = tab.alias(tab_Epetra_Export.getConst<T>(aid));
        if (!keepold) tab_Epetra_Export.remove(&aid);
        break;
    case CT_Epetra_Map_ID:
        newid = tab.alias(tab_Epetra_Map.getConst<T>(aid));
        if (!keepold) tab_Epetra_Map.remove(&aid);
        break;
    case CT_Epetra_BlockMap_ID:
        newid = tab.alias(tab_Epetra_BlockMap.getConst<T>(aid));
        if (!keepold) tab_Epetra_BlockMap.remove(&aid);
        break;
    case CT_Epetra_Import_ID:
        newid = tab.alias(tab_Epetra_Import.getConst<T>(aid));
        if (!keepold) tab_Epetra_Import.remove(&aid);
        break;
    case CT_Epetra_Time_ID:
        newid = tab.alias(tab_Epetra_Time.getConst<T>(aid));
        if (!keepold) tab_Epetra_Time.remove(&aid);
        break;
    case CT_Epetra_JadMatrix_ID:
        newid = tab.alias(tab_Epetra_JadMatrix.getConst<T>(aid));
        if (!keepold) tab_Epetra_JadMatrix.remove(&aid);
        break;
    case CT_Epetra_LinearProblem_ID:
        newid = tab.alias(tab_Epetra_LinearProblem.getConst<T>(aid));
        if (!keepold) tab_Epetra_LinearProblem.remove(&aid);
        break;
    case CT_Epetra_LAPACK_ID:
        newid = tab.alias(tab_Epetra_LAPACK.getConst<T>(aid));
        if (!keepold) tab_Epetra_LAPACK.remove(&aid);
        break;
    case CT_Teuchos_ParameterList_ID:
        newid = tab.alias(tab_Teuchos_ParameterList.getConst<T>(aid));
        if (!keepold) tab_Teuchos_ParameterList.remove(&aid);
        break;
#ifdef HAVE_CTRILINOS_AMESOS
    case CT_Amesos_BaseSolver_ID:
        newid = tab.alias(tab_Amesos_BaseSolver.getConst<T>(aid));
        if (!keepold) tab_Amesos_BaseSolver.remove(&aid);
        break;
#endif /* HAVE_CTRILINOS_AMESOS */
    case CT_Epetra_FECrsMatrix_ID:
        newid = tab.alias(tab_Epetra_FECrsMatrix.getConst<T>(aid));
        if (!keepold) tab_Epetra_FECrsMatrix.remove(&aid);
        break;
    case CT_Epetra_IntSerialDenseVector_ID:
        newid = tab.alias(tab_Epetra_IntSerialDenseVector.getConst<T>(aid));
        if (!keepold) tab_Epetra_IntSerialDenseVector.remove(&aid);
        break;
    case CT_Epetra_SerialDenseMatrix_ID:
        newid = tab.alias(tab_Epetra_SerialDenseMatrix.getConst<T>(aid));
        if (!keepold) tab_Epetra_SerialDenseMatrix.remove(&aid);
        break;
#ifdef HAVE_CTRILINOS_AZTECOO
    case CT_AztecOO_StatusTest_ID:
        newid = tab.alias(tab_AztecOO_StatusTest.getConst<T>(aid));
        if (!keepold) tab_AztecOO_StatusTest.remove(&aid);
        break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
    case CT_AztecOO_StatusTestCombo_ID:
        newid = tab.alias(tab_AztecOO_StatusTestCombo.getConst<T>(aid));
        if (!keepold) tab_AztecOO_StatusTestCombo.remove(&aid);
        break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
    case CT_AztecOO_StatusTestMaxIters_ID:
        newid = tab.alias(tab_AztecOO_StatusTestMaxIters.getConst<T>(aid));
        if (!keepold) tab_AztecOO_StatusTestMaxIters.remove(&aid);
        break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
    case CT_AztecOO_StatusTestResNorm_ID:
        newid = tab.alias(tab_AztecOO_StatusTestResNorm.getConst<T>(aid));
        if (!keepold) tab_AztecOO_StatusTestResNorm.remove(&aid);
        break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_IFPACK
    case CT_Ifpack_Preconditioner_ID:
        newid = tab.alias(tab_Ifpack_Preconditioner.getConst<T>(aid));
        if (!keepold) tab_Ifpack_Preconditioner.remove(&aid);
        break;
#endif /* HAVE_CTRILINOS_IFPACK */
    default:
        throw CTrilinosInvalidTypeError("invalid table id or non-polymorphic class");
    }

    return newid;
}

template <class T>
CTrilinos_Universal_ID_t TableRepos::store(T* pobj, bool owned)
{
    Table<T> *tab = 0;
    getTable(tab);
    return tab->store(pobj, owned);
}

template <class T>
CTrilinos_Universal_ID_t TableRepos::store(const T* pobj, bool owned)
{
    Table<T> *tab = 0;
    getTable(tab);
    return tab->store(pobj, owned);
}

template <class T>
void TableRepos::purge()
{
    Table<T> *tab = 0;
    getTable(tab);
    tab->purge();
}


} // namespace CTrilinos


#endif

