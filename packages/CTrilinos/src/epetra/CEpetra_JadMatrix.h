
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


/*! @file CEpetra_JadMatrix.h
 * @brief Wrappers for Epetra_JadMatrix */

/* True C header file! */


#ifndef CEPETRA_JADMATRIX_H
#define CEPETRA_JADMATRIX_H


#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif



/*! @name ID struct conversion functions */
/*@{*/

/*! @brief Changes the ID struct from the universal
   (generalized) struct type to the class-specific one.
*/
CT_Epetra_JadMatrix_ID_t Epetra_JadMatrix_Degeneralize ( 
  CTrilinos_Universal_ID_t id );

/*! @brief Changes the ID struct from the class-specific
   struct type to the universal (generalized) one.
*/
CTrilinos_Universal_ID_t Epetra_JadMatrix_Generalize ( 
  CT_Epetra_JadMatrix_ID_t id );

/*@}*/

/*! @name Epetra_JadMatrix constructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   Epetra_JadMatrix::Epetra_JadMatrix(const Epetra_RowMatrix & Matrix)
*/
CT_Epetra_JadMatrix_ID_t Epetra_JadMatrix_Create ( 
  CT_Epetra_RowMatrix_ID_t MatrixID );

/*@}*/

/*! @name Epetra_JadMatrix destructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   virtual Epetra_JadMatrix::~Epetra_JadMatrix()
*/
void Epetra_JadMatrix_Destroy ( CT_Epetra_JadMatrix_ID_t * selfID );

/*@}*/

/*! @name Epetra_JadMatrix member wrappers */
/*@{*/

/*! @brief Wrapper for 
   int Epetra_JadMatrix::UpdateValues(const Epetra_RowMatrix & Matrix, bool CheckStructure = false)
*/
int Epetra_JadMatrix_UpdateValues ( 
  CT_Epetra_JadMatrix_ID_t selfID, 
  CT_Epetra_RowMatrix_ID_t MatrixID, boolean CheckStructure );

/*! @brief Wrapper for 
   int Epetra_JadMatrix::ExtractMyRowCopy(int MyRow, int Length, int & NumEntries, double *Values, int * Indices) const
*/
int Epetra_JadMatrix_ExtractMyRowCopy ( 
  CT_Epetra_JadMatrix_ID_t selfID, int MyRow, int Length, 
  int * NumEntries, double * Values, int * Indices );

/*! @brief Wrapper for 
   int Epetra_JadMatrix::ExtractMyEntryView(int CurEntry, double * &Value, int & RowIndex, int & ColIndex)
*/
int Epetra_JadMatrix_ExtractMyEntryView ( 
  CT_Epetra_JadMatrix_ID_t selfID, int CurEntry, double * * Value, 
  int * RowIndex, int * ColIndex );

/*! @brief Wrapper for 
   int Epetra_JadMatrix::ExtractMyEntryView(int CurEntry, double const * & Value, int & RowIndex, int & ColIndex) const
*/
int Epetra_JadMatrix_ExtractMyEntryView_Const ( 
  CT_Epetra_JadMatrix_ID_t selfID, int CurEntry, 
  double const ** Value, int * RowIndex, int * ColIndex );

/*! @brief Wrapper for 
   int Epetra_JadMatrix::NumMyRowEntries(int MyRow, int & NumEntries) const
*/
int Epetra_JadMatrix_NumMyRowEntries ( 
  CT_Epetra_JadMatrix_ID_t selfID, int MyRow, int * NumEntries );

/*! @brief Wrapper for 
   int Epetra_JadMatrix::Multiply(bool TransA, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
*/
int Epetra_JadMatrix_Multiply ( 
  CT_Epetra_JadMatrix_ID_t selfID, boolean TransA, 
  CT_Epetra_MultiVector_ID_t XID, CT_Epetra_MultiVector_ID_t YID );

/*! @brief Wrapper for 
   int Epetra_JadMatrix::Solve(bool Upper, bool Trans, bool UnitDiagonal, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
*/
int Epetra_JadMatrix_Solve ( 
  CT_Epetra_JadMatrix_ID_t selfID, boolean Upper, boolean Trans, 
  boolean UnitDiagonal, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t YID );

/*@}*/


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CEPETRA_JADMATRIX_H */

