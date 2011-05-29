
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


/*! @file CTrilinos_enum_wrappers.h
 * @brief Defines enums used by CTrilinos to wrap other enums. */


#ifndef CTRILINOS_ENUM_WRAPPERS_H
#define CTRILINOS_ENUM_WRAPPERS_H


#include "CTrilinos_config.h"


#ifdef __cplusplus
extern "C" {
#endif


/*! Wrapper for Epetra_DataAccess */
typedef enum {
    CT_Epetra_DataAccess_E_Copy,
    CT_Epetra_DataAccess_E_View
} CT_Epetra_DataAccess_E_t;

/*! Wrapper for Epetra_CombineMode */
typedef enum {
    CT_Epetra_CombineMode_E_Add,
    CT_Epetra_CombineMode_E_Zero,
    CT_Epetra_CombineMode_E_Insert,
    CT_Epetra_CombineMode_E_InsertAdd,
    CT_Epetra_CombineMode_E_Average,
    CT_Epetra_CombineMode_E_AbsMax
} CT_Epetra_CombineMode_E_t;

/*! Wrapper for ProblemDifficultyLevel */
typedef enum {
    CT_ProblemDifficultyLevel_E_easy,
    CT_ProblemDifficultyLevel_E_moderate,
    CT_ProblemDifficultyLevel_E_hard,
    CT_ProblemDifficultyLevel_E_unsure
} CT_ProblemDifficultyLevel_E_t;

/*! Wrapper for Teuchos::EValidateUsed */
typedef enum {
    CT_EValidateUsed_E_VALIDATE_USED_ENABLED,
    CT_EValidateUsed_E_VALIDATE_USED_DISABLED
} CT_EValidateUsed_E_t;

/*! Wrapper for Teuchos::EValidateDefaults */
typedef enum {
    CT_EValidateDefaults_E_VALIDATE_DEFAULTS_ENABLED,
    CT_EValidateDefaults_E_VALIDATE_DEFAULTS_DISABLED
} CT_EValidateDefaults_E_t;

/*! Wrapper for Teuchos::CommandLineProcessor::EParseCommandLineReturn */
typedef enum {
    CT_EParseCommandLineReturn_E_PARSE_SUCCESSFUL = 0,
    CT_EParseCommandLineReturn_E_PARSE_HELP_PRINTED = 1,
    CT_EParseCommandLineReturn_E_PARSE_UNRECOGNIZED_OPTION = 2
} CT_EParseCommandLineReturn_E_t;

/*! Wrapper for Teuchos::CommandLineProcessor::EOptType */
typedef enum {
    CT_EOptType_E_OPT_NONE,
    CT_EOptType_E_OPT_BOOL_TRUE,
    CT_EOptType_E_OPT_BOOL_FALSE,
    CT_EOptType_E_OPT_INT,
    CT_EOptType_E_OPT_DOUBLE,
    CT_EOptType_E_OPT_STRING,
    CT_EOptType_E_OPT_ENUM_INT
} CT_EOptType_E_t;

#ifdef HAVE_CTRILINOS_AZTECOO
/*! Wrapper for AztecOO_StatusType */
typedef enum {
    CT_AztecOO_StatusType_E_Unchecked = 2,
    CT_AztecOO_StatusType_E_Unconverged = 1,
    CT_AztecOO_StatusType_E_Converged = 0,
    CT_AztecOO_StatusType_E_Failed = -1,
    CT_AztecOO_StatusType_E_NaN = -2,
    CT_AztecOO_StatusType_E_PartialFailed = -3
} CT_AztecOO_StatusType_E_t;
#endif /* HAVE_CTRILINOS_AZTECOO */

#ifdef HAVE_CTRILINOS_AZTECOO
/*! Wrapper for AztecOO_StatusTestCombo::ComboType */
typedef enum {
    CT_ComboType_E_AND,
    CT_ComboType_E_OR,
    CT_ComboType_E_SEQ
} CT_ComboType_E_t;
#endif /* HAVE_CTRILINOS_AZTECOO */

#ifdef HAVE_CTRILINOS_AZTECOO
/*! Wrapper for AztecOO_StatusTestResNorm::ResType */
typedef enum {
    CT_ResType_E_Implicit,
    CT_ResType_E_Explicit
} CT_ResType_E_t;
#endif /* HAVE_CTRILINOS_AZTECOO */

#ifdef HAVE_CTRILINOS_AZTECOO
/*! Wrapper for AztecOO_StatusTestResNorm::NormType */
typedef enum {
    CT_NormType_E_OneNorm,
    CT_NormType_E_TwoNorm,
    CT_NormType_E_InfNorm
} CT_NormType_E_t;
#endif /* HAVE_CTRILINOS_AZTECOO */

#ifdef HAVE_CTRILINOS_AZTECOO
/*! Wrapper for AztecOO_StatusTestResNorm::ScaleType */
typedef enum {
    CT_ScaleType_E_NormOfRHS,
    CT_ScaleType_E_NormOfInitRes,
    CT_ScaleType_E_None,
    CT_ScaleType_E_UserProvided
} CT_ScaleType_E_t;
#endif /* HAVE_CTRILINOS_AZTECOO */

#ifdef HAVE_CTRILINOS_IFPACK
/*! Wrapper for Ifpack::EPrecType */
typedef enum {
    CT_EPrecType_E_POINT_RELAXATION,
    CT_EPrecType_E_POINT_RELAXATION_STAND_ALONE,
    CT_EPrecType_E_BLOCK_RELAXATION,
    CT_EPrecType_E_BLOCK_RELAXATION_STAND_ALONE,
    CT_EPrecType_E_BLOCK_RELAXATION_STAND_ALONE_ILU,
    CT_EPrecType_E_BLOCK_RELAXATION_STAND_ALONE_AMESOS,
    CT_EPrecType_E_BLOCK_RELAXATION_AMESOS,
    CT_EPrecType_E_AMESOS,
    CT_EPrecType_E_AMESOS_STAND_ALONE,
    CT_EPrecType_E_IC,
    CT_EPrecType_E_IC_STAND_ALONE,
    CT_EPrecType_E_ICT,
    CT_EPrecType_E_ICT_STAND_ALONE,
    CT_EPrecType_E_ILU,
    CT_EPrecType_E_ILU_STAND_ALONE,
    CT_EPrecType_E_ILUT,
    CT_EPrecType_E_ILUT_STAND_ALONE,
    CT_EPrecType_E_SPARSKIT,
    CT_EPrecType_E_HIPS,
    CT_EPrecType_E_HYPRE,
    CT_EPrecType_E_CHEBYSHEV
} CT_EPrecType_E_t;
#endif /* HAVE_CTRILINOS_IFPACK */


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif
