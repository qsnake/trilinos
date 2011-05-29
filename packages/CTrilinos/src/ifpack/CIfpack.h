
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


#include "CTrilinos_config.h"


#ifdef HAVE_CTRILINOS_IFPACK



/*! @file CIfpack.h
 * @brief Wrappers for Ifpack */

/* True C header file! */


#ifndef CIFPACK_H
#define CIFPACK_H


#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif



/*! @name Ifpack constructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   Ifpack::Ifpack()
*/
CT_Ifpack_ID_t Ifpack_Create (  );

/*@}*/

/*! @name Ifpack destructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   Ifpack::~Ifpack()
*/
void Ifpack_Destroy ( CT_Ifpack_ID_t * selfID );

/*@}*/

/*! @name Ifpack member wrappers */
/*@{*/

/*! @brief Wrapper for 
   Ifpack_Preconditioner* Ifpack::Create(const string PrecType, Epetra_RowMatrix* Matrix, const int overlap = 0)
*/
CT_Ifpack_Preconditioner_ID_t Ifpack_CreatePreconditioner_UsingName ( 
  CT_Ifpack_ID_t selfID, const char PrecType[], 
  CT_Epetra_RowMatrix_ID_t MatrixID, const int overlap );

/*! @brief Wrapper for 
   int Ifpack::SetParameters(int argc, char* argv[], Teuchos::ParameterList& List, string& PrecType, int& Overlap)
*/
int Ifpack_SetParameters ( 
  CT_Ifpack_ID_t selfID, int argc, char * argv[], 
  CT_Teuchos_ParameterList_ID_t ListID, char * PrecType[], 
  int * Overlap );

/*@}*/

/*! @name Ifpack static function wrappers */
/*@{*/

/*! @brief Wrapper for 
   static const char* Ifpack::toString(const EPrecType precType)
*/
const char * Ifpack_toString ( const CT_EPrecType_E_t precType );

/*! @brief Wrapper for 
   static Ifpack_Preconditioner* Ifpack::Create( EPrecType PrecType, Epetra_RowMatrix* Matrix, const int overlap = 0 )
*/
CT_Ifpack_Preconditioner_ID_t Ifpack_CreatePreconditioner_UsingType ( 
  CT_EPrecType_E_t PrecType, CT_Epetra_RowMatrix_ID_t MatrixID, 
  const int overlap );

/*@}*/


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CIFPACK_H */

#endif /* HAVE_CTRILINOS_IFPACK */


