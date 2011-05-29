
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


/*! @file CTrilinos_TableRepos.cpp
 * @brief Central table repository for CTrilinos. */


#include "CTrilinos_TableRepos.hpp"


namespace CTrilinos {


TableRepos & tableRepos(  )
{
    static TableRepos tr;
    return tr;
}


TableRepos::TableRepos() :
    tab_Epetra_Distributor(CT_Epetra_Distributor_ID),
    tab_Epetra_SerialComm(CT_Epetra_SerialComm_ID),
    tab_Epetra_BLAS(CT_Epetra_BLAS_ID),
    tab_Epetra_Comm(CT_Epetra_Comm_ID),
    tab_Epetra_Operator(CT_Epetra_Operator_ID),
    tab_Epetra_MultiVector(CT_Epetra_MultiVector_ID),
    tab_Epetra_OffsetIndex(CT_Epetra_OffsetIndex_ID),
    tab_Epetra_Object(CT_Epetra_Object_ID),
    tab_Epetra_RowMatrix(CT_Epetra_RowMatrix_ID),
    tab_Epetra_CompObject(CT_Epetra_CompObject_ID),
    tab_Epetra_Directory(CT_Epetra_Directory_ID),
    tab_Epetra_Flops(CT_Epetra_Flops_ID),
    tab_Epetra_SrcDistObject(CT_Epetra_SrcDistObject_ID),
#ifdef HAVE_MPI
    tab_Epetra_MpiComm(CT_Epetra_MpiComm_ID),
#endif /* HAVE_MPI */
    tab_Epetra_CrsMatrix(CT_Epetra_CrsMatrix_ID),
    tab_Epetra_CrsGraph(CT_Epetra_CrsGraph_ID),
    tab_Epetra_DistObject(CT_Epetra_DistObject_ID),
    tab_Epetra_Vector(CT_Epetra_Vector_ID),
    tab_Epetra_Export(CT_Epetra_Export_ID),
    tab_Epetra_Map(CT_Epetra_Map_ID),
    tab_Epetra_BlockMap(CT_Epetra_BlockMap_ID),
    tab_Epetra_Import(CT_Epetra_Import_ID),
    tab_Epetra_Time(CT_Epetra_Time_ID),
    tab_Epetra_JadMatrix(CT_Epetra_JadMatrix_ID),
    tab_Epetra_LinearProblem(CT_Epetra_LinearProblem_ID),
    tab_Epetra_LAPACK(CT_Epetra_LAPACK_ID),
    tab_Teuchos_ParameterList(CT_Teuchos_ParameterList_ID),
#ifdef HAVE_CTRILINOS_AMESOS
    tab_Amesos_BaseSolver(CT_Amesos_BaseSolver_ID),
#endif /* HAVE_CTRILINOS_AMESOS */
    tab_Epetra_FECrsMatrix(CT_Epetra_FECrsMatrix_ID),
    tab_Epetra_IntSerialDenseVector(CT_Epetra_IntSerialDenseVector_ID),
    tab_Epetra_SerialDenseMatrix(CT_Epetra_SerialDenseMatrix_ID),
#ifdef HAVE_CTRILINOS_AZTECOO
    tab_AztecOO_StatusTest(CT_AztecOO_StatusTest_ID),
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
    tab_AztecOO_StatusTestCombo(CT_AztecOO_StatusTestCombo_ID),
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
    tab_AztecOO_StatusTestMaxIters(CT_AztecOO_StatusTestMaxIters_ID),
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
    tab_AztecOO_StatusTestResNorm(CT_AztecOO_StatusTestResNorm_ID),
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_IFPACK
    tab_Ifpack_Preconditioner(CT_Ifpack_Preconditioner_ID),
#endif /* HAVE_CTRILINOS_IFPACK */
    call_me_lazy(true)
{
}

TableRepos::~TableRepos()
{
}

void TableRepos::getTable(Table<Epetra_Distributor> *& tab)
{
    tab = &tab_Epetra_Distributor;
}
void TableRepos::getTable(Table<Epetra_SerialComm> *& tab)
{
    tab = &tab_Epetra_SerialComm;
}
void TableRepos::getTable(Table<Epetra_BLAS> *& tab)
{
    tab = &tab_Epetra_BLAS;
}
void TableRepos::getTable(Table<Epetra_Comm> *& tab)
{
    tab = &tab_Epetra_Comm;
}
void TableRepos::getTable(Table<Epetra_Operator> *& tab)
{
    tab = &tab_Epetra_Operator;
}
void TableRepos::getTable(Table<Epetra_MultiVector> *& tab)
{
    tab = &tab_Epetra_MultiVector;
}
void TableRepos::getTable(Table<Epetra_OffsetIndex> *& tab)
{
    tab = &tab_Epetra_OffsetIndex;
}
void TableRepos::getTable(Table<Epetra_Object> *& tab)
{
    tab = &tab_Epetra_Object;
}
void TableRepos::getTable(Table<Epetra_RowMatrix> *& tab)
{
    tab = &tab_Epetra_RowMatrix;
}
void TableRepos::getTable(Table<Epetra_CompObject> *& tab)
{
    tab = &tab_Epetra_CompObject;
}
void TableRepos::getTable(Table<Epetra_Directory> *& tab)
{
    tab = &tab_Epetra_Directory;
}
void TableRepos::getTable(Table<Epetra_Flops> *& tab)
{
    tab = &tab_Epetra_Flops;
}
void TableRepos::getTable(Table<Epetra_SrcDistObject> *& tab)
{
    tab = &tab_Epetra_SrcDistObject;
}
#ifdef HAVE_MPI
void TableRepos::getTable(Table<Epetra_MpiComm> *& tab)
{
    tab = &tab_Epetra_MpiComm;
}
#endif /* HAVE_MPI */
void TableRepos::getTable(Table<Epetra_CrsMatrix> *& tab)
{
    tab = &tab_Epetra_CrsMatrix;
}
void TableRepos::getTable(Table<Epetra_CrsGraph> *& tab)
{
    tab = &tab_Epetra_CrsGraph;
}
void TableRepos::getTable(Table<Epetra_DistObject> *& tab)
{
    tab = &tab_Epetra_DistObject;
}
void TableRepos::getTable(Table<Epetra_Vector> *& tab)
{
    tab = &tab_Epetra_Vector;
}
void TableRepos::getTable(Table<Epetra_Export> *& tab)
{
    tab = &tab_Epetra_Export;
}
void TableRepos::getTable(Table<Epetra_Map> *& tab)
{
    tab = &tab_Epetra_Map;
}
void TableRepos::getTable(Table<Epetra_BlockMap> *& tab)
{
    tab = &tab_Epetra_BlockMap;
}
void TableRepos::getTable(Table<Epetra_Import> *& tab)
{
    tab = &tab_Epetra_Import;
}
void TableRepos::getTable(Table<Epetra_Time> *& tab)
{
    tab = &tab_Epetra_Time;
}
void TableRepos::getTable(Table<Epetra_JadMatrix> *& tab)
{
    tab = &tab_Epetra_JadMatrix;
}
void TableRepos::getTable(Table<Epetra_LinearProblem> *& tab)
{
    tab = &tab_Epetra_LinearProblem;
}
void TableRepos::getTable(Table<Epetra_LAPACK> *& tab)
{
    tab = &tab_Epetra_LAPACK;
}
void TableRepos::getTable(Table<Teuchos::ParameterList> *& tab)
{
    tab = &tab_Teuchos_ParameterList;
}
#ifdef HAVE_CTRILINOS_AMESOS
void TableRepos::getTable(Table<Amesos_BaseSolver> *& tab)
{
    tab = &tab_Amesos_BaseSolver;
}
#endif /* HAVE_CTRILINOS_AMESOS */
void TableRepos::getTable(Table<Epetra_FECrsMatrix> *& tab)
{
    tab = &tab_Epetra_FECrsMatrix;
}
void TableRepos::getTable(Table<Epetra_IntSerialDenseVector> *& tab)
{
    tab = &tab_Epetra_IntSerialDenseVector;
}
void TableRepos::getTable(Table<Epetra_SerialDenseMatrix> *& tab)
{
    tab = &tab_Epetra_SerialDenseMatrix;
}
#ifdef HAVE_CTRILINOS_AZTECOO
void TableRepos::getTable(Table<AztecOO_StatusTest> *& tab)
{
    tab = &tab_AztecOO_StatusTest;
}
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
void TableRepos::getTable(Table<AztecOO_StatusTestCombo> *& tab)
{
    tab = &tab_AztecOO_StatusTestCombo;
}
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
void TableRepos::getTable(Table<AztecOO_StatusTestMaxIters> *& tab)
{
    tab = &tab_AztecOO_StatusTestMaxIters;
}
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
void TableRepos::getTable(Table<AztecOO_StatusTestResNorm> *& tab)
{
    tab = &tab_AztecOO_StatusTestResNorm;
}
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_IFPACK
void TableRepos::getTable(Table<Ifpack_Preconditioner> *& tab)
{
    tab = &tab_Ifpack_Preconditioner;
}
#endif /* HAVE_CTRILINOS_IFPACK */

CTrilinos_Universal_ID_t TableRepos::alias(
    CTrilinos_Universal_ID_t aid, CTrilinos_Table_ID_t tab, bool keepold)
{
    CTrilinos_Universal_ID_t newid;

    switch (tab) {
    case CT_Epetra_Distributor_ID:
        newid = (aid.is_const ? do_alias_const(tab_Epetra_Distributor, aid, keepold)
                              : do_alias(tab_Epetra_Distributor, aid, keepold));
        break;
    case CT_Epetra_SerialComm_ID:
        newid = (aid.is_const ? do_alias_const(tab_Epetra_SerialComm, aid, keepold)
                              : do_alias(tab_Epetra_SerialComm, aid, keepold));
        break;
    case CT_Epetra_BLAS_ID:
        newid = (aid.is_const ? do_alias_const(tab_Epetra_BLAS, aid, keepold)
                              : do_alias(tab_Epetra_BLAS, aid, keepold));
        break;
    case CT_Epetra_Comm_ID:
        newid = (aid.is_const ? do_alias_const(tab_Epetra_Comm, aid, keepold)
                              : do_alias(tab_Epetra_Comm, aid, keepold));
        break;
    case CT_Epetra_Operator_ID:
        newid = (aid.is_const ? do_alias_const(tab_Epetra_Operator, aid, keepold)
                              : do_alias(tab_Epetra_Operator, aid, keepold));
        break;
    case CT_Epetra_MultiVector_ID:
        newid = (aid.is_const ? do_alias_const(tab_Epetra_MultiVector, aid, keepold)
                              : do_alias(tab_Epetra_MultiVector, aid, keepold));
        break;
    case CT_Epetra_OffsetIndex_ID:
        newid = (aid.is_const ? do_alias_const(tab_Epetra_OffsetIndex, aid, keepold)
                              : do_alias(tab_Epetra_OffsetIndex, aid, keepold));
        break;
    case CT_Epetra_Object_ID:
        newid = (aid.is_const ? do_alias_const(tab_Epetra_Object, aid, keepold)
                              : do_alias(tab_Epetra_Object, aid, keepold));
        break;
    case CT_Epetra_RowMatrix_ID:
        newid = (aid.is_const ? do_alias_const(tab_Epetra_RowMatrix, aid, keepold)
                              : do_alias(tab_Epetra_RowMatrix, aid, keepold));
        break;
    case CT_Epetra_CompObject_ID:
        newid = (aid.is_const ? do_alias_const(tab_Epetra_CompObject, aid, keepold)
                              : do_alias(tab_Epetra_CompObject, aid, keepold));
        break;
    case CT_Epetra_Directory_ID:
        newid = (aid.is_const ? do_alias_const(tab_Epetra_Directory, aid, keepold)
                              : do_alias(tab_Epetra_Directory, aid, keepold));
        break;
    case CT_Epetra_Flops_ID:
        newid = (aid.is_const ? do_alias_const(tab_Epetra_Flops, aid, keepold)
                              : do_alias(tab_Epetra_Flops, aid, keepold));
        break;
    case CT_Epetra_SrcDistObject_ID:
        newid = (aid.is_const ? do_alias_const(tab_Epetra_SrcDistObject, aid, keepold)
                              : do_alias(tab_Epetra_SrcDistObject, aid, keepold));
        break;
#ifdef HAVE_MPI
    case CT_Epetra_MpiComm_ID:
        newid = (aid.is_const ? do_alias_const(tab_Epetra_MpiComm, aid, keepold)
                              : do_alias(tab_Epetra_MpiComm, aid, keepold));
        break;
#endif /* HAVE_MPI */
    case CT_Epetra_CrsMatrix_ID:
        newid = (aid.is_const ? do_alias_const(tab_Epetra_CrsMatrix, aid, keepold)
                              : do_alias(tab_Epetra_CrsMatrix, aid, keepold));
        break;
    case CT_Epetra_CrsGraph_ID:
        newid = (aid.is_const ? do_alias_const(tab_Epetra_CrsGraph, aid, keepold)
                              : do_alias(tab_Epetra_CrsGraph, aid, keepold));
        break;
    case CT_Epetra_DistObject_ID:
        newid = (aid.is_const ? do_alias_const(tab_Epetra_DistObject, aid, keepold)
                              : do_alias(tab_Epetra_DistObject, aid, keepold));
        break;
    case CT_Epetra_Vector_ID:
        newid = (aid.is_const ? do_alias_const(tab_Epetra_Vector, aid, keepold)
                              : do_alias(tab_Epetra_Vector, aid, keepold));
        break;
    case CT_Epetra_Export_ID:
        newid = (aid.is_const ? do_alias_const(tab_Epetra_Export, aid, keepold)
                              : do_alias(tab_Epetra_Export, aid, keepold));
        break;
    case CT_Epetra_Map_ID:
        newid = (aid.is_const ? do_alias_const(tab_Epetra_Map, aid, keepold)
                              : do_alias(tab_Epetra_Map, aid, keepold));
        break;
    case CT_Epetra_BlockMap_ID:
        newid = (aid.is_const ? do_alias_const(tab_Epetra_BlockMap, aid, keepold)
                              : do_alias(tab_Epetra_BlockMap, aid, keepold));
        break;
    case CT_Epetra_Import_ID:
        newid = (aid.is_const ? do_alias_const(tab_Epetra_Import, aid, keepold)
                              : do_alias(tab_Epetra_Import, aid, keepold));
        break;
    case CT_Epetra_Time_ID:
        newid = (aid.is_const ? do_alias_const(tab_Epetra_Time, aid, keepold)
                              : do_alias(tab_Epetra_Time, aid, keepold));
        break;
    case CT_Epetra_JadMatrix_ID:
        newid = (aid.is_const ? do_alias_const(tab_Epetra_JadMatrix, aid, keepold)
                              : do_alias(tab_Epetra_JadMatrix, aid, keepold));
        break;
    case CT_Epetra_LinearProblem_ID:
        newid = (aid.is_const ? do_alias_const(tab_Epetra_LinearProblem, aid, keepold)
                              : do_alias(tab_Epetra_LinearProblem, aid, keepold));
        break;
    case CT_Epetra_LAPACK_ID:
        newid = (aid.is_const ? do_alias_const(tab_Epetra_LAPACK, aid, keepold)
                              : do_alias(tab_Epetra_LAPACK, aid, keepold));
        break;
    case CT_Teuchos_ParameterList_ID:
        newid = (aid.is_const ? do_alias_const(tab_Teuchos_ParameterList, aid, keepold)
                              : do_alias(tab_Teuchos_ParameterList, aid, keepold));
        break;
#ifdef HAVE_CTRILINOS_AMESOS
    case CT_Amesos_BaseSolver_ID:
        newid = (aid.is_const ? do_alias_const(tab_Amesos_BaseSolver, aid, keepold)
                              : do_alias(tab_Amesos_BaseSolver, aid, keepold));
        break;
#endif /* HAVE_CTRILINOS_AMESOS */
    case CT_Epetra_FECrsMatrix_ID:
        newid = (aid.is_const ? do_alias_const(tab_Epetra_FECrsMatrix, aid, keepold)
                              : do_alias(tab_Epetra_FECrsMatrix, aid, keepold));
        break;
    case CT_Epetra_IntSerialDenseVector_ID:
        newid = (aid.is_const ? do_alias_const(tab_Epetra_IntSerialDenseVector, aid, keepold)
                              : do_alias(tab_Epetra_IntSerialDenseVector, aid, keepold));
        break;
    case CT_Epetra_SerialDenseMatrix_ID:
        newid = (aid.is_const ? do_alias_const(tab_Epetra_SerialDenseMatrix, aid, keepold)
                              : do_alias(tab_Epetra_SerialDenseMatrix, aid, keepold));
        break;
#ifdef HAVE_CTRILINOS_AZTECOO
    case CT_AztecOO_StatusTest_ID:
        newid = (aid.is_const ? do_alias_const(tab_AztecOO_StatusTest, aid, keepold)
                              : do_alias(tab_AztecOO_StatusTest, aid, keepold));
        break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
    case CT_AztecOO_StatusTestCombo_ID:
        newid = (aid.is_const ? do_alias_const(tab_AztecOO_StatusTestCombo, aid, keepold)
                              : do_alias(tab_AztecOO_StatusTestCombo, aid, keepold));
        break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
    case CT_AztecOO_StatusTestMaxIters_ID:
        newid = (aid.is_const ? do_alias_const(tab_AztecOO_StatusTestMaxIters, aid, keepold)
                              : do_alias(tab_AztecOO_StatusTestMaxIters, aid, keepold));
        break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
    case CT_AztecOO_StatusTestResNorm_ID:
        newid = (aid.is_const ? do_alias_const(tab_AztecOO_StatusTestResNorm, aid, keepold)
                              : do_alias(tab_AztecOO_StatusTestResNorm, aid, keepold));
        break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_IFPACK
    case CT_Ifpack_Preconditioner_ID:
        newid = (aid.is_const ? do_alias_const(tab_Ifpack_Preconditioner, aid, keepold)
                              : do_alias(tab_Ifpack_Preconditioner, aid, keepold));
        break;
#endif /* HAVE_CTRILINOS_IFPACK */
    default:
        throw CTrilinosInvalidTypeError("invalid table id or non-polymorphic class");
    }

    return newid;
}

void TableRepos::remove(CTrilinos_Universal_ID_t * aid)
{
    switch (aid->table) {
    case CT_Epetra_Distributor_ID:
        tab_Epetra_Distributor.remove(aid);
        break;
    case CT_Epetra_SerialComm_ID:
        tab_Epetra_SerialComm.remove(aid);
        break;
    case CT_Epetra_BLAS_ID:
        tab_Epetra_BLAS.remove(aid);
        break;
    case CT_Epetra_Comm_ID:
        tab_Epetra_Comm.remove(aid);
        break;
    case CT_Epetra_Operator_ID:
        tab_Epetra_Operator.remove(aid);
        break;
    case CT_Epetra_MultiVector_ID:
        tab_Epetra_MultiVector.remove(aid);
        break;
    case CT_Epetra_OffsetIndex_ID:
        tab_Epetra_OffsetIndex.remove(aid);
        break;
    case CT_Epetra_Object_ID:
        tab_Epetra_Object.remove(aid);
        break;
    case CT_Epetra_RowMatrix_ID:
        tab_Epetra_RowMatrix.remove(aid);
        break;
    case CT_Epetra_CompObject_ID:
        tab_Epetra_CompObject.remove(aid);
        break;
    case CT_Epetra_Directory_ID:
        tab_Epetra_Directory.remove(aid);
        break;
    case CT_Epetra_Flops_ID:
        tab_Epetra_Flops.remove(aid);
        break;
    case CT_Epetra_SrcDistObject_ID:
        tab_Epetra_SrcDistObject.remove(aid);
        break;
#ifdef HAVE_MPI
    case CT_Epetra_MpiComm_ID:
        tab_Epetra_MpiComm.remove(aid);
        break;
#endif /* HAVE_MPI */
    case CT_Epetra_CrsMatrix_ID:
        tab_Epetra_CrsMatrix.remove(aid);
        break;
    case CT_Epetra_CrsGraph_ID:
        tab_Epetra_CrsGraph.remove(aid);
        break;
    case CT_Epetra_DistObject_ID:
        tab_Epetra_DistObject.remove(aid);
        break;
    case CT_Epetra_Vector_ID:
        tab_Epetra_Vector.remove(aid);
        break;
    case CT_Epetra_Export_ID:
        tab_Epetra_Export.remove(aid);
        break;
    case CT_Epetra_Map_ID:
        tab_Epetra_Map.remove(aid);
        break;
    case CT_Epetra_BlockMap_ID:
        tab_Epetra_BlockMap.remove(aid);
        break;
    case CT_Epetra_Import_ID:
        tab_Epetra_Import.remove(aid);
        break;
    case CT_Epetra_Time_ID:
        tab_Epetra_Time.remove(aid);
        break;
    case CT_Epetra_JadMatrix_ID:
        tab_Epetra_JadMatrix.remove(aid);
        break;
    case CT_Epetra_LinearProblem_ID:
        tab_Epetra_LinearProblem.remove(aid);
        break;
    case CT_Epetra_LAPACK_ID:
        tab_Epetra_LAPACK.remove(aid);
        break;
    case CT_Teuchos_ParameterList_ID:
        tab_Teuchos_ParameterList.remove(aid);
        break;
#ifdef HAVE_CTRILINOS_AMESOS
    case CT_Amesos_BaseSolver_ID:
        tab_Amesos_BaseSolver.remove(aid);
        break;
#endif /* HAVE_CTRILINOS_AMESOS */
    case CT_Epetra_FECrsMatrix_ID:
        tab_Epetra_FECrsMatrix.remove(aid);
        break;
    case CT_Epetra_IntSerialDenseVector_ID:
        tab_Epetra_IntSerialDenseVector.remove(aid);
        break;
    case CT_Epetra_SerialDenseMatrix_ID:
        tab_Epetra_SerialDenseMatrix.remove(aid);
        break;
#ifdef HAVE_CTRILINOS_AZTECOO
    case CT_AztecOO_StatusTest_ID:
        tab_AztecOO_StatusTest.remove(aid);
        break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
    case CT_AztecOO_StatusTestCombo_ID:
        tab_AztecOO_StatusTestCombo.remove(aid);
        break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
    case CT_AztecOO_StatusTestMaxIters_ID:
        tab_AztecOO_StatusTestMaxIters.remove(aid);
        break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
    case CT_AztecOO_StatusTestResNorm_ID:
        tab_AztecOO_StatusTestResNorm.remove(aid);
        break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_IFPACK
    case CT_Ifpack_Preconditioner_ID:
        tab_Ifpack_Preconditioner.remove(aid);
        break;
#endif /* HAVE_CTRILINOS_IFPACK */
    default:
        throw CTrilinosInvalidTypeError("invalid table id");
    }
}

void TableRepos::purgeAll()
{
    tab_Epetra_Distributor.purge();
    tab_Epetra_SerialComm.purge();
    tab_Epetra_BLAS.purge();
    tab_Epetra_Comm.purge();
    tab_Epetra_Operator.purge();
    tab_Epetra_MultiVector.purge();
    tab_Epetra_OffsetIndex.purge();
    tab_Epetra_Object.purge();
    tab_Epetra_RowMatrix.purge();
    tab_Epetra_CompObject.purge();
    tab_Epetra_Directory.purge();
    tab_Epetra_Flops.purge();
    tab_Epetra_SrcDistObject.purge();
#ifdef HAVE_MPI
    tab_Epetra_MpiComm.purge();
#endif /* HAVE_MPI */
    tab_Epetra_CrsMatrix.purge();
    tab_Epetra_CrsGraph.purge();
    tab_Epetra_DistObject.purge();
    tab_Epetra_Vector.purge();
    tab_Epetra_Export.purge();
    tab_Epetra_Map.purge();
    tab_Epetra_BlockMap.purge();
    tab_Epetra_Import.purge();
    tab_Epetra_Time.purge();
    tab_Epetra_JadMatrix.purge();
    tab_Epetra_LinearProblem.purge();
    tab_Epetra_LAPACK.purge();
    tab_Teuchos_ParameterList.purge();
#ifdef HAVE_CTRILINOS_AMESOS
    tab_Amesos_BaseSolver.purge();
#endif /* HAVE_CTRILINOS_AMESOS */
    tab_Epetra_FECrsMatrix.purge();
    tab_Epetra_IntSerialDenseVector.purge();
    tab_Epetra_SerialDenseMatrix.purge();
#ifdef HAVE_CTRILINOS_AZTECOO
    tab_AztecOO_StatusTest.purge();
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
    tab_AztecOO_StatusTestCombo.purge();
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
    tab_AztecOO_StatusTestMaxIters.purge();
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
    tab_AztecOO_StatusTestResNorm.purge();
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_IFPACK
    tab_Ifpack_Preconditioner.purge();
#endif /* HAVE_CTRILINOS_IFPACK */
}

bool TableRepos::typeCheck(CTrilinos_Universal_ID_t aid, CTrilinos_Table_ID_t type)
{
    if (aid.table == type)
        return true;

    try {
        if (aid.is_const) {
            switch (type) {
            case CT_Epetra_Distributor_ID:
                getConst<Epetra_Distributor>(aid);
                break;
            case CT_Epetra_SerialComm_ID:
                getConst<Epetra_SerialComm>(aid);
                break;
            case CT_Epetra_BLAS_ID:
                getConst<Epetra_BLAS>(aid);
                break;
            case CT_Epetra_Comm_ID:
                getConst<Epetra_Comm>(aid);
                break;
            case CT_Epetra_Operator_ID:
                getConst<Epetra_Operator>(aid);
                break;
            case CT_Epetra_MultiVector_ID:
                getConst<Epetra_MultiVector>(aid);
                break;
            case CT_Epetra_OffsetIndex_ID:
                getConst<Epetra_OffsetIndex>(aid);
                break;
            case CT_Epetra_Object_ID:
                getConst<Epetra_Object>(aid);
                break;
            case CT_Epetra_RowMatrix_ID:
                getConst<Epetra_RowMatrix>(aid);
                break;
            case CT_Epetra_CompObject_ID:
                getConst<Epetra_CompObject>(aid);
                break;
            case CT_Epetra_Directory_ID:
                getConst<Epetra_Directory>(aid);
                break;
            case CT_Epetra_Flops_ID:
                getConst<Epetra_Flops>(aid);
                break;
            case CT_Epetra_SrcDistObject_ID:
                getConst<Epetra_SrcDistObject>(aid);
                break;
#ifdef HAVE_MPI
            case CT_Epetra_MpiComm_ID:
                getConst<Epetra_MpiComm>(aid);
                break;
#endif /* HAVE_MPI */
            case CT_Epetra_CrsMatrix_ID:
                getConst<Epetra_CrsMatrix>(aid);
                break;
            case CT_Epetra_CrsGraph_ID:
                getConst<Epetra_CrsGraph>(aid);
                break;
            case CT_Epetra_DistObject_ID:
                getConst<Epetra_DistObject>(aid);
                break;
            case CT_Epetra_Vector_ID:
                getConst<Epetra_Vector>(aid);
                break;
            case CT_Epetra_Export_ID:
                getConst<Epetra_Export>(aid);
                break;
            case CT_Epetra_Map_ID:
                getConst<Epetra_Map>(aid);
                break;
            case CT_Epetra_BlockMap_ID:
                getConst<Epetra_BlockMap>(aid);
                break;
            case CT_Epetra_Import_ID:
                getConst<Epetra_Import>(aid);
                break;
            case CT_Epetra_Time_ID:
                getConst<Epetra_Time>(aid);
                break;
            case CT_Epetra_JadMatrix_ID:
                getConst<Epetra_JadMatrix>(aid);
                break;
            case CT_Epetra_LinearProblem_ID:
                getConst<Epetra_LinearProblem>(aid);
                break;
            case CT_Epetra_LAPACK_ID:
                getConst<Epetra_LAPACK>(aid);
                break;
            case CT_Teuchos_ParameterList_ID:
                getConst<Teuchos::ParameterList>(aid);
                break;
#ifdef HAVE_CTRILINOS_AMESOS
            case CT_Amesos_BaseSolver_ID:
                getConst<Amesos_BaseSolver>(aid);
                break;
#endif /* HAVE_CTRILINOS_AMESOS */
            case CT_Epetra_FECrsMatrix_ID:
                getConst<Epetra_FECrsMatrix>(aid);
                break;
            case CT_Epetra_IntSerialDenseVector_ID:
                getConst<Epetra_IntSerialDenseVector>(aid);
                break;
            case CT_Epetra_SerialDenseMatrix_ID:
                getConst<Epetra_SerialDenseMatrix>(aid);
                break;
#ifdef HAVE_CTRILINOS_AZTECOO
            case CT_AztecOO_StatusTest_ID:
                getConst<AztecOO_StatusTest>(aid);
                break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
            case CT_AztecOO_StatusTestCombo_ID:
                getConst<AztecOO_StatusTestCombo>(aid);
                break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
            case CT_AztecOO_StatusTestMaxIters_ID:
                getConst<AztecOO_StatusTestMaxIters>(aid);
                break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
            case CT_AztecOO_StatusTestResNorm_ID:
                getConst<AztecOO_StatusTestResNorm>(aid);
                break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_IFPACK
            case CT_Ifpack_Preconditioner_ID:
                getConst<Ifpack_Preconditioner>(aid);
                break;
#endif /* HAVE_CTRILINOS_IFPACK */
            default:
                return false;
            }
        } else {
            switch (type) {
            case CT_Epetra_Distributor_ID:
                get<Epetra_Distributor>(aid);
                break;
            case CT_Epetra_SerialComm_ID:
                get<Epetra_SerialComm>(aid);
                break;
            case CT_Epetra_BLAS_ID:
                get<Epetra_BLAS>(aid);
                break;
            case CT_Epetra_Comm_ID:
                get<Epetra_Comm>(aid);
                break;
            case CT_Epetra_Operator_ID:
                get<Epetra_Operator>(aid);
                break;
            case CT_Epetra_MultiVector_ID:
                get<Epetra_MultiVector>(aid);
                break;
            case CT_Epetra_OffsetIndex_ID:
                get<Epetra_OffsetIndex>(aid);
                break;
            case CT_Epetra_Object_ID:
                get<Epetra_Object>(aid);
                break;
            case CT_Epetra_RowMatrix_ID:
                get<Epetra_RowMatrix>(aid);
                break;
            case CT_Epetra_CompObject_ID:
                get<Epetra_CompObject>(aid);
                break;
            case CT_Epetra_Directory_ID:
                get<Epetra_Directory>(aid);
                break;
            case CT_Epetra_Flops_ID:
                get<Epetra_Flops>(aid);
                break;
            case CT_Epetra_SrcDistObject_ID:
                get<Epetra_SrcDistObject>(aid);
                break;
#ifdef HAVE_MPI
            case CT_Epetra_MpiComm_ID:
                get<Epetra_MpiComm>(aid);
                break;
#endif /* HAVE_MPI */
            case CT_Epetra_CrsMatrix_ID:
                get<Epetra_CrsMatrix>(aid);
                break;
            case CT_Epetra_CrsGraph_ID:
                get<Epetra_CrsGraph>(aid);
                break;
            case CT_Epetra_DistObject_ID:
                get<Epetra_DistObject>(aid);
                break;
            case CT_Epetra_Vector_ID:
                get<Epetra_Vector>(aid);
                break;
            case CT_Epetra_Export_ID:
                get<Epetra_Export>(aid);
                break;
            case CT_Epetra_Map_ID:
                get<Epetra_Map>(aid);
                break;
            case CT_Epetra_BlockMap_ID:
                get<Epetra_BlockMap>(aid);
                break;
            case CT_Epetra_Import_ID:
                get<Epetra_Import>(aid);
                break;
            case CT_Epetra_Time_ID:
                get<Epetra_Time>(aid);
                break;
            case CT_Epetra_JadMatrix_ID:
                get<Epetra_JadMatrix>(aid);
                break;
            case CT_Epetra_LinearProblem_ID:
                get<Epetra_LinearProblem>(aid);
                break;
            case CT_Epetra_LAPACK_ID:
                get<Epetra_LAPACK>(aid);
                break;
            case CT_Teuchos_ParameterList_ID:
                get<Teuchos::ParameterList>(aid);
                break;
#ifdef HAVE_CTRILINOS_AMESOS
            case CT_Amesos_BaseSolver_ID:
                get<Amesos_BaseSolver>(aid);
                break;
#endif /* HAVE_CTRILINOS_AMESOS */
            case CT_Epetra_FECrsMatrix_ID:
                get<Epetra_FECrsMatrix>(aid);
                break;
            case CT_Epetra_IntSerialDenseVector_ID:
                get<Epetra_IntSerialDenseVector>(aid);
                break;
            case CT_Epetra_SerialDenseMatrix_ID:
                get<Epetra_SerialDenseMatrix>(aid);
                break;
#ifdef HAVE_CTRILINOS_AZTECOO
            case CT_AztecOO_StatusTest_ID:
                get<AztecOO_StatusTest>(aid);
                break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
            case CT_AztecOO_StatusTestCombo_ID:
                get<AztecOO_StatusTestCombo>(aid);
                break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
            case CT_AztecOO_StatusTestMaxIters_ID:
                get<AztecOO_StatusTestMaxIters>(aid);
                break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
            case CT_AztecOO_StatusTestResNorm_ID:
                get<AztecOO_StatusTestResNorm>(aid);
                break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_IFPACK
            case CT_Ifpack_Preconditioner_ID:
                get<Ifpack_Preconditioner>(aid);
                break;
#endif /* HAVE_CTRILINOS_IFPACK */
            default:
                return false;
            }
        }
    } catch (...) {
        return false;
    }

    return true;
}


} // namespace CTrilinos


