
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


/*! @file CEpetra_SerialDenseMatrix.h
 * @brief Wrappers for Epetra_SerialDenseMatrix */

/* True C header file! */


#ifndef CEPETRA_SERIALDENSEMATRIX_H
#define CEPETRA_SERIALDENSEMATRIX_H


#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif



/*! @name ID struct conversion functions */
/*@{*/

/*! @brief Changes the ID struct from the universal
   (generalized) struct type to the class-specific one.
*/
CT_Epetra_SerialDenseMatrix_ID_t Epetra_SerialDenseMatrix_Degeneralize ( 
  CTrilinos_Universal_ID_t id );

/*! @brief Changes the ID struct from the class-specific
   struct type to the universal (generalized) one.
*/
CTrilinos_Universal_ID_t Epetra_SerialDenseMatrix_Generalize ( 
  CT_Epetra_SerialDenseMatrix_ID_t id );

/*@}*/

/*! @name Epetra_SerialDenseMatrix constructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   Epetra_SerialDenseMatrix::Epetra_SerialDenseMatrix(bool set_object_label=true)
*/
CT_Epetra_SerialDenseMatrix_ID_t Epetra_SerialDenseMatrix_Create_Empty ( 
  boolean set_object_label );

/*! @brief Wrapper for 
   Epetra_SerialDenseMatrix::Epetra_SerialDenseMatrix(int NumRows, int NumCols, bool set_object_label=true)
*/
CT_Epetra_SerialDenseMatrix_ID_t Epetra_SerialDenseMatrix_Create ( 
  int NumRows, int NumCols, boolean set_object_label );

/*! @brief Wrapper for 
   Epetra_SerialDenseMatrix::Epetra_SerialDenseMatrix(Epetra_DataAccess CV, double* A_in, int LDA_in, int NumRows, int NumCols, bool set_object_label=true)
*/
CT_Epetra_SerialDenseMatrix_ID_t Epetra_SerialDenseMatrix_Create_FromArray ( 
  CT_Epetra_DataAccess_E_t CV, double * A_in, int LDA_in, 
  int NumRows, int NumCols, boolean set_object_label );

/*! @brief Wrapper for 
   Epetra_SerialDenseMatrix::Epetra_SerialDenseMatrix(const Epetra_SerialDenseMatrix& Source)
*/
CT_Epetra_SerialDenseMatrix_ID_t Epetra_SerialDenseMatrix_Duplicate ( 
  CT_Epetra_SerialDenseMatrix_ID_t SourceID );

/*@}*/

/*! @name Epetra_SerialDenseMatrix destructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   virtual Epetra_SerialDenseMatrix::~Epetra_SerialDenseMatrix()
*/
void Epetra_SerialDenseMatrix_Destroy ( 
  CT_Epetra_SerialDenseMatrix_ID_t * selfID );

/*@}*/

/*! @name Epetra_SerialDenseMatrix member wrappers */
/*@{*/

/*! @brief Wrapper for 
   int Epetra_SerialDenseMatrix::Shape(int NumRows, int NumCols)
*/
int Epetra_SerialDenseMatrix_Shape ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID, int NumRows, int NumCols );

/*! @brief Wrapper for 
   int Epetra_SerialDenseMatrix::Reshape(int NumRows, int NumCols)
*/
int Epetra_SerialDenseMatrix_Reshape ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID, int NumRows, int NumCols );

/*! @brief Wrapper for 
   int Epetra_SerialDenseMatrix::Multiply(char TransA, char TransB, double ScalarAB, const Epetra_SerialDenseMatrix& A, const Epetra_SerialDenseMatrix& B, double ScalarThis)
*/
int Epetra_SerialDenseMatrix_Multiply_Matrix ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID, char TransA, char TransB, 
  double ScalarAB, CT_Epetra_SerialDenseMatrix_ID_t AID, 
  CT_Epetra_SerialDenseMatrix_ID_t BID, double ScalarThis );

/*! @brief Wrapper for 
   int Epetra_SerialDenseMatrix::Multiply(bool transA, const Epetra_SerialDenseMatrix& x, Epetra_SerialDenseMatrix& y)
*/
int Epetra_SerialDenseMatrix_Multiply_Vector ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID, boolean transA, 
  CT_Epetra_SerialDenseMatrix_ID_t xID, 
  CT_Epetra_SerialDenseMatrix_ID_t yID );

/*! @brief Wrapper for 
   int Epetra_SerialDenseMatrix::Scale(double ScalarA)
*/
int Epetra_SerialDenseMatrix_Scale ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID, double ScalarA );

/*! @brief Wrapper for 
   virtual double Epetra_SerialDenseMatrix::NormOne() const
*/
double Epetra_SerialDenseMatrix_NormOne ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID );

/*! @brief Wrapper for 
   virtual double Epetra_SerialDenseMatrix::NormInf() const
*/
double Epetra_SerialDenseMatrix_NormInf ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_SerialDenseMatrix::Random()
*/
int Epetra_SerialDenseMatrix_Random ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_SerialDenseMatrix::M() const
*/
int Epetra_SerialDenseMatrix_M ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_SerialDenseMatrix::N() const
*/
int Epetra_SerialDenseMatrix_N ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID );

/*! @brief Wrapper for 
   double* Epetra_SerialDenseMatrix::A() const
*/
double * Epetra_SerialDenseMatrix_A_Const ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID );

/*! @brief Wrapper for 
   double* Epetra_SerialDenseMatrix::A()
*/
double * Epetra_SerialDenseMatrix_A ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_SerialDenseMatrix::LDA() const
*/
int Epetra_SerialDenseMatrix_LDA ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID );

/*! @brief Wrapper for 
   Epetra_DataAccess Epetra_SerialDenseMatrix::CV() const
*/
CT_Epetra_DataAccess_E_t Epetra_SerialDenseMatrix_CV ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID );

/*! @brief Wrapper for 
   virtual double Epetra_SerialDenseMatrix::OneNorm() const
*/
double Epetra_SerialDenseMatrix_OneNorm ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID );

/*! @brief Wrapper for 
   virtual double Epetra_SerialDenseMatrix::InfNorm() const
*/
double Epetra_SerialDenseMatrix_InfNorm ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID );

/*! @brief Wrapper for 
   virtual int Epetra_SerialDenseMatrix::SetUseTranspose(bool UseTranspose_in)
*/
int Epetra_SerialDenseMatrix_SetUseTranspose ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID, boolean UseTranspose_in );

/*! @brief Wrapper for 
   virtual int Epetra_SerialDenseMatrix::Apply(const Epetra_SerialDenseMatrix& X, Epetra_SerialDenseMatrix& Y)
*/
int Epetra_SerialDenseMatrix_Apply ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID, 
  CT_Epetra_SerialDenseMatrix_ID_t XID, 
  CT_Epetra_SerialDenseMatrix_ID_t YID );

/*! @brief Wrapper for 
   virtual int Epetra_SerialDenseMatrix::ApplyInverse(const Epetra_SerialDenseMatrix & X, Epetra_SerialDenseMatrix & Y)
*/
int Epetra_SerialDenseMatrix_ApplyInverse ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID, 
  CT_Epetra_SerialDenseMatrix_ID_t XID, 
  CT_Epetra_SerialDenseMatrix_ID_t YID );

/*! @brief Wrapper for 
   virtual const char * Epetra_SerialDenseMatrix::Label() const
*/
const char * Epetra_SerialDenseMatrix_Label ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID );

/*! @brief Wrapper for 
   virtual bool Epetra_SerialDenseMatrix::UseTranspose() const
*/
boolean Epetra_SerialDenseMatrix_UseTranspose ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID );

/*! @brief Wrapper for 
   virtual bool Epetra_SerialDenseMatrix::HasNormInf() const
*/
boolean Epetra_SerialDenseMatrix_HasNormInf ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID );

/*! @brief Wrapper for 
   virtual int Epetra_SerialDenseMatrix::RowDim() const
*/
int Epetra_SerialDenseMatrix_RowDim ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID );

/*! @brief Wrapper for 
   virtual int Epetra_SerialDenseMatrix::ColDim() const
*/
int Epetra_SerialDenseMatrix_ColDim ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID );

/*@}*/

/*! @name Epetra_SerialDenseMatrix operator wrappers */
/*@{*/

/*! @brief Wrapper for 
   Epetra_SerialDenseMatrix & Epetra_SerialDenseMatrix::operator= (const Epetra_SerialDenseMatrix& Source)
*/
void Epetra_SerialDenseMatrix_Assign ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID, 
  CT_Epetra_SerialDenseMatrix_ID_t SourceID );

/*! @brief Wrapper for 
   bool Epetra_SerialDenseMatrix::operator==(const Epetra_SerialDenseMatrix& rhs) const
*/
boolean Epetra_SerialDenseMatrix_IsEqual ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID, 
  CT_Epetra_SerialDenseMatrix_ID_t rhsID );

/*! @brief Wrapper for 
   bool Epetra_SerialDenseMatrix::operator!=(const Epetra_SerialDenseMatrix& rhs) const
*/
boolean Epetra_SerialDenseMatrix_NotEqual ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID, 
  CT_Epetra_SerialDenseMatrix_ID_t rhsID );

/*! @brief Wrapper for 
   Epetra_SerialDenseMatrix & Epetra_SerialDenseMatrix::operator+= (const Epetra_SerialDenseMatrix& Source)
*/
void Epetra_SerialDenseMatrix_AddTo ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID, 
  CT_Epetra_SerialDenseMatrix_ID_t SourceID );

/*! @brief Wrapper for 
   double& Epetra_SerialDenseMatrix::operator() (int RowIndex, int ColIndex)
*/
void Epetra_SerialDenseMatrix_setElement ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID, int RowIndex, 
  int ColIndex, double * value );

/*! @brief Wrapper for 
   const double& Epetra_SerialDenseMatrix::operator() (int RowIndex, int ColIndex) const
*/
double Epetra_SerialDenseMatrix_getElement ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID, int RowIndex, 
  int ColIndex );

/*! @brief Wrapper for 
   const double* Epetra_SerialDenseMatrix::operator[] (int ColIndex) const
*/
const double * Epetra_SerialDenseMatrix_getColumn ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID, int ColIndex );

/*@}*/


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CEPETRA_SERIALDENSEMATRIX_H */

