
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


/*! @file CTrilinos_flex_enums.h
 * @brief Defines structs and enums needed for CTrilinos. */


#ifndef CTRILINOS_FLEX_ENUMS_H
#define CTRILINOS_FLEX_ENUMS_H


#include "CTrilinos_config.h"
#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif


typedef union {
    CTrilinos_Universal_ID_t universal;
    CT_Epetra_Distributor_ID_t Epetra_Distributor;
} CT_Epetra_Distributor_ID_Flex_t;

typedef union {
    CTrilinos_Universal_ID_t universal;
    CT_Epetra_SerialComm_ID_t Epetra_SerialComm;
    CT_Epetra_Comm_ID_t Epetra_Comm;
    CT_Epetra_Object_ID_t Epetra_Object;
} CT_Epetra_SerialComm_ID_Flex_t;

typedef union {
    CTrilinos_Universal_ID_t universal;
    CT_Epetra_BLAS_ID_t Epetra_BLAS;
} CT_Epetra_BLAS_ID_Flex_t;

typedef union {
    CTrilinos_Universal_ID_t universal;
    CT_Epetra_Comm_ID_t Epetra_Comm;
} CT_Epetra_Comm_ID_Flex_t;

typedef union {
    CTrilinos_Universal_ID_t universal;
    CT_Epetra_Operator_ID_t Epetra_Operator;
} CT_Epetra_Operator_ID_Flex_t;

typedef union {
    CTrilinos_Universal_ID_t universal;
    CT_Epetra_MultiVector_ID_t Epetra_MultiVector;
    CT_Epetra_BLAS_ID_t Epetra_BLAS;
    CT_Epetra_CompObject_ID_t Epetra_CompObject;
    CT_Epetra_DistObject_ID_t Epetra_DistObject;
    CT_Epetra_Object_ID_t Epetra_Object;
    CT_Epetra_SrcDistObject_ID_t Epetra_SrcDistObject;
} CT_Epetra_MultiVector_ID_Flex_t;

typedef union {
    CTrilinos_Universal_ID_t universal;
    CT_Epetra_OffsetIndex_ID_t Epetra_OffsetIndex;
    CT_Epetra_Object_ID_t Epetra_Object;
} CT_Epetra_OffsetIndex_ID_Flex_t;

typedef union {
    CTrilinos_Universal_ID_t universal;
    CT_Epetra_Object_ID_t Epetra_Object;
} CT_Epetra_Object_ID_Flex_t;

typedef union {
    CTrilinos_Universal_ID_t universal;
    CT_Epetra_RowMatrix_ID_t Epetra_RowMatrix;
    CT_Epetra_Operator_ID_t Epetra_Operator;
    CT_Epetra_SrcDistObject_ID_t Epetra_SrcDistObject;
} CT_Epetra_RowMatrix_ID_Flex_t;

typedef union {
    CTrilinos_Universal_ID_t universal;
    CT_Epetra_CompObject_ID_t Epetra_CompObject;
} CT_Epetra_CompObject_ID_Flex_t;

typedef union {
    CTrilinos_Universal_ID_t universal;
    CT_Epetra_Directory_ID_t Epetra_Directory;
} CT_Epetra_Directory_ID_Flex_t;

typedef union {
    CTrilinos_Universal_ID_t universal;
    CT_Epetra_Flops_ID_t Epetra_Flops;
} CT_Epetra_Flops_ID_Flex_t;

typedef union {
    CTrilinos_Universal_ID_t universal;
    CT_Epetra_SrcDistObject_ID_t Epetra_SrcDistObject;
} CT_Epetra_SrcDistObject_ID_Flex_t;

#ifdef HAVE_MPI
typedef union {
    CTrilinos_Universal_ID_t universal;
    CT_Epetra_MpiComm_ID_t Epetra_MpiComm;
#ifdef HAVE_MPI
    CT_Epetra_Comm_ID_t Epetra_Comm;
#endif /* HAVE_MPI */
#ifdef HAVE_MPI
    CT_Epetra_Object_ID_t Epetra_Object;
#endif /* HAVE_MPI */
} CT_Epetra_MpiComm_ID_Flex_t;
#endif /* HAVE_MPI */

typedef union {
    CTrilinos_Universal_ID_t universal;
    CT_Epetra_CrsMatrix_ID_t Epetra_CrsMatrix;
    CT_Epetra_BLAS_ID_t Epetra_BLAS;
    CT_Epetra_CompObject_ID_t Epetra_CompObject;
    CT_Epetra_DistObject_ID_t Epetra_DistObject;
    CT_Epetra_Object_ID_t Epetra_Object;
    CT_Epetra_Operator_ID_t Epetra_Operator;
    CT_Epetra_RowMatrix_ID_t Epetra_RowMatrix;
    CT_Epetra_SrcDistObject_ID_t Epetra_SrcDistObject;
} CT_Epetra_CrsMatrix_ID_Flex_t;

typedef union {
    CTrilinos_Universal_ID_t universal;
    CT_Epetra_CrsGraph_ID_t Epetra_CrsGraph;
    CT_Epetra_DistObject_ID_t Epetra_DistObject;
    CT_Epetra_Object_ID_t Epetra_Object;
    CT_Epetra_SrcDistObject_ID_t Epetra_SrcDistObject;
} CT_Epetra_CrsGraph_ID_Flex_t;

typedef union {
    CTrilinos_Universal_ID_t universal;
    CT_Epetra_DistObject_ID_t Epetra_DistObject;
    CT_Epetra_Object_ID_t Epetra_Object;
    CT_Epetra_SrcDistObject_ID_t Epetra_SrcDistObject;
} CT_Epetra_DistObject_ID_Flex_t;

typedef union {
    CTrilinos_Universal_ID_t universal;
    CT_Epetra_Vector_ID_t Epetra_Vector;
    CT_Epetra_BLAS_ID_t Epetra_BLAS;
    CT_Epetra_CompObject_ID_t Epetra_CompObject;
    CT_Epetra_DistObject_ID_t Epetra_DistObject;
    CT_Epetra_MultiVector_ID_t Epetra_MultiVector;
    CT_Epetra_Object_ID_t Epetra_Object;
    CT_Epetra_SrcDistObject_ID_t Epetra_SrcDistObject;
} CT_Epetra_Vector_ID_Flex_t;

typedef union {
    CTrilinos_Universal_ID_t universal;
    CT_Epetra_Export_ID_t Epetra_Export;
    CT_Epetra_Object_ID_t Epetra_Object;
} CT_Epetra_Export_ID_Flex_t;

typedef union {
    CTrilinos_Universal_ID_t universal;
    CT_Epetra_Map_ID_t Epetra_Map;
    CT_Epetra_BlockMap_ID_t Epetra_BlockMap;
    CT_Epetra_Object_ID_t Epetra_Object;
} CT_Epetra_Map_ID_Flex_t;

typedef union {
    CTrilinos_Universal_ID_t universal;
    CT_Epetra_BlockMap_ID_t Epetra_BlockMap;
    CT_Epetra_Object_ID_t Epetra_Object;
} CT_Epetra_BlockMap_ID_Flex_t;

typedef union {
    CTrilinos_Universal_ID_t universal;
    CT_Epetra_Import_ID_t Epetra_Import;
    CT_Epetra_Object_ID_t Epetra_Object;
} CT_Epetra_Import_ID_Flex_t;

typedef union {
    CTrilinos_Universal_ID_t universal;
    CT_Epetra_Time_ID_t Epetra_Time;
    CT_Epetra_Object_ID_t Epetra_Object;
} CT_Epetra_Time_ID_Flex_t;

typedef union {
    CTrilinos_Universal_ID_t universal;
    CT_Epetra_JadMatrix_ID_t Epetra_JadMatrix;
    CT_Epetra_CompObject_ID_t Epetra_CompObject;
    CT_Epetra_Object_ID_t Epetra_Object;
    CT_Epetra_Operator_ID_t Epetra_Operator;
    CT_Epetra_RowMatrix_ID_t Epetra_RowMatrix;
    CT_Epetra_SrcDistObject_ID_t Epetra_SrcDistObject;
} CT_Epetra_JadMatrix_ID_Flex_t;

typedef union {
    CTrilinos_Universal_ID_t universal;
    CT_Epetra_LinearProblem_ID_t Epetra_LinearProblem;
} CT_Epetra_LinearProblem_ID_Flex_t;

typedef union {
    CTrilinos_Universal_ID_t universal;
    CT_Epetra_LAPACK_ID_t Epetra_LAPACK;
} CT_Epetra_LAPACK_ID_Flex_t;

typedef union {
    CTrilinos_Universal_ID_t universal;
    CT_Teuchos_CommandLineProcessor_ID_t Teuchos_CommandLineProcessor;
} CT_Teuchos_CommandLineProcessor_ID_Flex_t;

typedef union {
    CTrilinos_Universal_ID_t universal;
    CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList;
} CT_Teuchos_ParameterList_ID_Flex_t;

typedef union {
    CTrilinos_Universal_ID_t universal;
    CT_Teuchos_ParameterEntry_ID_t Teuchos_ParameterEntry;
} CT_Teuchos_ParameterEntry_ID_Flex_t;

typedef union {
    CTrilinos_Universal_ID_t universal;
    CT_Teuchos_any_ID_t Teuchos_any;
} CT_Teuchos_any_ID_Flex_t;

#ifdef HAVE_CTRILINOS_AMESOS
typedef union {
    CTrilinos_Universal_ID_t universal;
    CT_Amesos_BaseSolver_ID_t Amesos_BaseSolver;
} CT_Amesos_BaseSolver_ID_Flex_t;
#endif /* HAVE_CTRILINOS_AMESOS */

#ifdef HAVE_CTRILINOS_AMESOS
typedef union {
    CTrilinos_Universal_ID_t universal;
    CT_Amesos_ID_t Amesos;
} CT_Amesos_ID_Flex_t;
#endif /* HAVE_CTRILINOS_AMESOS */

typedef union {
    CTrilinos_Universal_ID_t universal;
    CT_Epetra_FECrsMatrix_ID_t Epetra_FECrsMatrix;
    CT_Epetra_BLAS_ID_t Epetra_BLAS;
    CT_Epetra_CompObject_ID_t Epetra_CompObject;
    CT_Epetra_CrsMatrix_ID_t Epetra_CrsMatrix;
    CT_Epetra_DistObject_ID_t Epetra_DistObject;
    CT_Epetra_Object_ID_t Epetra_Object;
    CT_Epetra_Operator_ID_t Epetra_Operator;
    CT_Epetra_RowMatrix_ID_t Epetra_RowMatrix;
    CT_Epetra_SrcDistObject_ID_t Epetra_SrcDistObject;
} CT_Epetra_FECrsMatrix_ID_Flex_t;

typedef union {
    CTrilinos_Universal_ID_t universal;
    CT_Epetra_IntSerialDenseVector_ID_t Epetra_IntSerialDenseVector;
    CT_Epetra_Object_ID_t Epetra_Object;
} CT_Epetra_IntSerialDenseVector_ID_Flex_t;

typedef union {
    CTrilinos_Universal_ID_t universal;
    CT_Epetra_SerialDenseMatrix_ID_t Epetra_SerialDenseMatrix;
    CT_Epetra_BLAS_ID_t Epetra_BLAS;
    CT_Epetra_CompObject_ID_t Epetra_CompObject;
    CT_Epetra_Object_ID_t Epetra_Object;
} CT_Epetra_SerialDenseMatrix_ID_Flex_t;

#ifdef HAVE_CTRILINOS_AZTECOO
typedef union {
    CTrilinos_Universal_ID_t universal;
    CT_AztecOO_ID_t AztecOO;
} CT_AztecOO_ID_Flex_t;
#endif /* HAVE_CTRILINOS_AZTECOO */

#ifdef HAVE_CTRILINOS_AZTECOO
typedef union {
    CTrilinos_Universal_ID_t universal;
    CT_AztecOO_StatusTest_ID_t AztecOO_StatusTest;
} CT_AztecOO_StatusTest_ID_Flex_t;
#endif /* HAVE_CTRILINOS_AZTECOO */

#ifdef HAVE_CTRILINOS_AZTECOO
typedef union {
    CTrilinos_Universal_ID_t universal;
    CT_AztecOO_StatusTestCombo_ID_t AztecOO_StatusTestCombo;
#ifdef HAVE_CTRILINOS_AZTECOO
    CT_AztecOO_StatusTest_ID_t AztecOO_StatusTest;
#endif /* HAVE_CTRILINOS_AZTECOO */
} CT_AztecOO_StatusTestCombo_ID_Flex_t;
#endif /* HAVE_CTRILINOS_AZTECOO */

#ifdef HAVE_CTRILINOS_AZTECOO
typedef union {
    CTrilinos_Universal_ID_t universal;
    CT_AztecOO_StatusTestMaxIters_ID_t AztecOO_StatusTestMaxIters;
#ifdef HAVE_CTRILINOS_AZTECOO
    CT_AztecOO_StatusTest_ID_t AztecOO_StatusTest;
#endif /* HAVE_CTRILINOS_AZTECOO */
} CT_AztecOO_StatusTestMaxIters_ID_Flex_t;
#endif /* HAVE_CTRILINOS_AZTECOO */

#ifdef HAVE_CTRILINOS_AZTECOO
typedef union {
    CTrilinos_Universal_ID_t universal;
    CT_AztecOO_StatusTestResNorm_ID_t AztecOO_StatusTestResNorm;
#ifdef HAVE_CTRILINOS_AZTECOO
    CT_AztecOO_StatusTest_ID_t AztecOO_StatusTest;
#endif /* HAVE_CTRILINOS_AZTECOO */
} CT_AztecOO_StatusTestResNorm_ID_Flex_t;
#endif /* HAVE_CTRILINOS_AZTECOO */

#ifdef HAVE_CTRILINOS_IFPACK
typedef union {
    CTrilinos_Universal_ID_t universal;
    CT_Ifpack_ID_t Ifpack;
} CT_Ifpack_ID_Flex_t;
#endif /* HAVE_CTRILINOS_IFPACK */

#ifdef HAVE_CTRILINOS_IFPACK
typedef union {
    CTrilinos_Universal_ID_t universal;
    CT_Ifpack_Preconditioner_ID_t Ifpack_Preconditioner;
#ifdef HAVE_CTRILINOS_IFPACK
    CT_Epetra_Operator_ID_t Epetra_Operator;
#endif /* HAVE_CTRILINOS_IFPACK */
} CT_Ifpack_Preconditioner_ID_Flex_t;
#endif /* HAVE_CTRILINOS_IFPACK */


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif
