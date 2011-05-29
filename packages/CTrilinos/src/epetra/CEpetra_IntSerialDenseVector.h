
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


/*! @file CEpetra_IntSerialDenseVector.h
 * @brief Wrappers for Epetra_IntSerialDenseVector */

/* True C header file! */


#ifndef CEPETRA_INTSERIALDENSEVECTOR_H
#define CEPETRA_INTSERIALDENSEVECTOR_H


#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif



/*! @name ID struct conversion functions */
/*@{*/

/*! @brief Changes the ID struct from the universal
   (generalized) struct type to the class-specific one.
*/
CT_Epetra_IntSerialDenseVector_ID_t Epetra_IntSerialDenseVector_Degeneralize ( 
  CTrilinos_Universal_ID_t id );

/*! @brief Changes the ID struct from the class-specific
   struct type to the universal (generalized) one.
*/
CTrilinos_Universal_ID_t Epetra_IntSerialDenseVector_Generalize ( 
  CT_Epetra_IntSerialDenseVector_ID_t id );

/*@}*/

/*! @name Epetra_IntSerialDenseVector constructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   Epetra_IntSerialDenseVector::Epetra_IntSerialDenseVector()
*/
CT_Epetra_IntSerialDenseVector_ID_t Epetra_IntSerialDenseVector_Create_Empty ( 
   );

/*! @brief Wrapper for 
   Epetra_IntSerialDenseVector::Epetra_IntSerialDenseVector(int Length_in)
*/
CT_Epetra_IntSerialDenseVector_ID_t Epetra_IntSerialDenseVector_Create ( 
  int Length_in );

/*! @brief Wrapper for 
   Epetra_IntSerialDenseVector::Epetra_IntSerialDenseVector(Epetra_DataAccess CV_in, int* Values_in, int Length_in)
*/
CT_Epetra_IntSerialDenseVector_ID_t Epetra_IntSerialDenseVector_Create_FromArray ( 
  CT_Epetra_DataAccess_E_t CV_in, int * Values_in, int Length_in );

/*! @brief Wrapper for 
   Epetra_IntSerialDenseVector::Epetra_IntSerialDenseVector(const Epetra_IntSerialDenseVector& Source)
*/
CT_Epetra_IntSerialDenseVector_ID_t Epetra_IntSerialDenseVector_Duplicate ( 
  CT_Epetra_IntSerialDenseVector_ID_t SourceID );

/*@}*/

/*! @name Epetra_IntSerialDenseVector destructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   virtual Epetra_IntSerialDenseVector::~Epetra_IntSerialDenseVector()
*/
void Epetra_IntSerialDenseVector_Destroy ( 
  CT_Epetra_IntSerialDenseVector_ID_t * selfID );

/*@}*/

/*! @name Epetra_IntSerialDenseVector member wrappers */
/*@{*/

/*! @brief Wrapper for 
   int Epetra_IntSerialDenseVector::Size(int Length_in)
*/
int Epetra_IntSerialDenseVector_Size ( 
  CT_Epetra_IntSerialDenseVector_ID_t selfID, int Length_in );

/*! @brief Wrapper for 
   int Epetra_IntSerialDenseVector::Resize(int Length_in)
*/
int Epetra_IntSerialDenseVector_Resize ( 
  CT_Epetra_IntSerialDenseVector_ID_t selfID, int Length_in );

/*! @brief Wrapper for 
   int Epetra_IntSerialDenseVector::Random()
*/
int Epetra_IntSerialDenseVector_Random ( 
  CT_Epetra_IntSerialDenseVector_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_IntSerialDenseVector::Length() const
*/
int Epetra_IntSerialDenseVector_Length ( 
  CT_Epetra_IntSerialDenseVector_ID_t selfID );

/*! @brief Wrapper for 
   int* Epetra_IntSerialDenseVector::Values()
*/
int * Epetra_IntSerialDenseVector_Values ( 
  CT_Epetra_IntSerialDenseVector_ID_t selfID );

/*! @brief Wrapper for 
   const int* Epetra_IntSerialDenseVector::Values() const
*/
const int * Epetra_IntSerialDenseVector_Values_Const ( 
  CT_Epetra_IntSerialDenseVector_ID_t selfID );

/*! @brief Wrapper for 
   Epetra_DataAccess Epetra_IntSerialDenseVector::CV() const
*/
CT_Epetra_DataAccess_E_t Epetra_IntSerialDenseVector_CV ( 
  CT_Epetra_IntSerialDenseVector_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_IntSerialDenseVector::MakeViewOf(const Epetra_IntSerialDenseVector& Source)
*/
int Epetra_IntSerialDenseVector_MakeViewOf ( 
  CT_Epetra_IntSerialDenseVector_ID_t selfID, 
  CT_Epetra_IntSerialDenseVector_ID_t SourceID );

/*@}*/

/*! @name Epetra_IntSerialDenseVector operator wrappers */
/*@{*/

/*! @brief Wrapper for 
   int& Epetra_IntSerialDenseVector::operator() (int Index)
*/
void Epetra_IntSerialDenseVector_setElement ( 
  CT_Epetra_IntSerialDenseVector_ID_t selfID, int Index, 
  int * value );

/*! @brief Wrapper for 
   const int& Epetra_IntSerialDenseVector::operator() (int Index) const
*/
int Epetra_IntSerialDenseVector_getElement ( 
  CT_Epetra_IntSerialDenseVector_ID_t selfID, int Index );

/*! @brief Wrapper for 
   int& Epetra_IntSerialDenseVector::operator[] (int Index)
*/
void Epetra_IntSerialDenseVector_setElement_Bracket ( 
  CT_Epetra_IntSerialDenseVector_ID_t selfID, int Index, 
  int * value );

/*! @brief Wrapper for 
   const int& Epetra_IntSerialDenseVector::operator[] (int Index) const
*/
int Epetra_IntSerialDenseVector_getElement_Bracket ( 
  CT_Epetra_IntSerialDenseVector_ID_t selfID, int Index );

/*! @brief Wrapper for 
   Epetra_IntSerialDenseVector& Epetra_IntSerialDenseVector::operator= (const Epetra_IntSerialDenseVector& Source)
*/
void Epetra_IntSerialDenseVector_Assign ( 
  CT_Epetra_IntSerialDenseVector_ID_t selfID, 
  CT_Epetra_IntSerialDenseVector_ID_t SourceID );

/*@}*/


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CEPETRA_INTSERIALDENSEVECTOR_H */

