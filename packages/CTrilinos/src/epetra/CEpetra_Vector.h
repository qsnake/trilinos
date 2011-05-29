
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


/*! @file CEpetra_Vector.h
 * @brief Wrappers for Epetra_Vector */

/* True C header file! */


#ifndef CEPETRA_VECTOR_H
#define CEPETRA_VECTOR_H


#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif



/*! @name ID struct conversion functions */
/*@{*/

/*! @brief Changes the ID struct from the universal
   (generalized) struct type to the class-specific one.
*/
CT_Epetra_Vector_ID_t Epetra_Vector_Degeneralize ( 
  CTrilinos_Universal_ID_t id );

/*! @brief Changes the ID struct from the class-specific
   struct type to the universal (generalized) one.
*/
CTrilinos_Universal_ID_t Epetra_Vector_Generalize ( 
  CT_Epetra_Vector_ID_t id );

/*@}*/

/*! @name Epetra_Vector constructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   Epetra_Vector::Epetra_Vector(const Epetra_BlockMap& Map, bool zeroOut = true)
*/
CT_Epetra_Vector_ID_t Epetra_Vector_Create ( 
  CT_Epetra_BlockMap_ID_t MapID, boolean zeroOut );

/*! @brief Wrapper for 
   Epetra_Vector::Epetra_Vector(const Epetra_Vector& Source)
*/
CT_Epetra_Vector_ID_t Epetra_Vector_Duplicate ( 
  CT_Epetra_Vector_ID_t SourceID );

/*! @brief Wrapper for 
   Epetra_Vector::Epetra_Vector(Epetra_DataAccess CV, const Epetra_BlockMap& Map, double *V)
*/
CT_Epetra_Vector_ID_t Epetra_Vector_Create_FromArray ( 
  CT_Epetra_DataAccess_E_t CV, CT_Epetra_BlockMap_ID_t MapID, 
  double * V );

/*! @brief Wrapper for 
   Epetra_Vector::Epetra_Vector(Epetra_DataAccess CV, const Epetra_MultiVector& Source, int Index)
*/
CT_Epetra_Vector_ID_t Epetra_Vector_FromSource ( 
  CT_Epetra_DataAccess_E_t CV, CT_Epetra_MultiVector_ID_t SourceID, 
  int Index );

/*@}*/

/*! @name Epetra_Vector destructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   virtual Epetra_Vector::~Epetra_Vector()
*/
void Epetra_Vector_Destroy ( CT_Epetra_Vector_ID_t * selfID );

/*@}*/

/*! @name Epetra_Vector member wrappers */
/*@{*/

/*! @brief Wrapper for 
   int Epetra_Vector::ReplaceGlobalValues(int NumEntries, double * Values, int * Indices)
*/
int Epetra_Vector_ReplaceGlobalValues ( 
  CT_Epetra_Vector_ID_t selfID, int NumEntries, double * Values, 
  int * Indices );

/*! @brief Wrapper for 
   int Epetra_Vector::ReplaceMyValues(int NumEntries, double * Values, int * Indices)
*/
int Epetra_Vector_ReplaceMyValues ( 
  CT_Epetra_Vector_ID_t selfID, int NumEntries, double * Values, 
  int * Indices );

/*! @brief Wrapper for 
   int Epetra_Vector::SumIntoGlobalValues(int NumEntries, double * Values, int * Indices)
*/
int Epetra_Vector_SumIntoGlobalValues ( 
  CT_Epetra_Vector_ID_t selfID, int NumEntries, double * Values, 
  int * Indices );

/*! @brief Wrapper for 
   int Epetra_Vector::SumIntoMyValues(int NumEntries, double * Values, int * Indices)
*/
int Epetra_Vector_SumIntoMyValues ( 
  CT_Epetra_Vector_ID_t selfID, int NumEntries, double * Values, 
  int * Indices );

/*! @brief Wrapper for 
   int Epetra_Vector::ReplaceGlobalValues(int NumEntries, int BlockOffset, double * Values, int * Indices)
*/
int Epetra_Vector_ReplaceGlobalValues_BlockPos ( 
  CT_Epetra_Vector_ID_t selfID, int NumEntries, int BlockOffset, 
  double * Values, int * Indices );

/*! @brief Wrapper for 
   int Epetra_Vector::ReplaceMyValues(int NumEntries, int BlockOffset, double * Values, int * Indices)
*/
int Epetra_Vector_ReplaceMyValues_BlockPos ( 
  CT_Epetra_Vector_ID_t selfID, int NumEntries, int BlockOffset, 
  double * Values, int * Indices );

/*! @brief Wrapper for 
   int Epetra_Vector::SumIntoGlobalValues(int NumEntries, int BlockOffset, double * Values, int * Indices)
*/
int Epetra_Vector_SumIntoGlobalValues_BlockPos ( 
  CT_Epetra_Vector_ID_t selfID, int NumEntries, int BlockOffset, 
  double * Values, int * Indices );

/*! @brief Wrapper for 
   int Epetra_Vector::SumIntoMyValues(int NumEntries, int BlockOffset, double * Values, int * Indices)
*/
int Epetra_Vector_SumIntoMyValues_BlockPos ( 
  CT_Epetra_Vector_ID_t selfID, int NumEntries, int BlockOffset, 
  double * Values, int * Indices );

/*! @brief Wrapper for 
   int Epetra_Vector::ExtractCopy(double *V) const
*/
int Epetra_Vector_ExtractCopy ( 
  CT_Epetra_Vector_ID_t selfID, double * V );

/*! @brief Wrapper for 
   int Epetra_Vector::ExtractView(double **V) const
*/
int Epetra_Vector_ExtractView ( 
  CT_Epetra_Vector_ID_t selfID, double ** V );

/*@}*/

/*! @name Epetra_Vector operator wrappers */
/*@{*/

/*! @brief Wrapper for 
   const double& Epetra_Vector::operator[] (int index) const
*/
double Epetra_Vector_getElement ( 
  CT_Epetra_Vector_ID_t selfID, int index );

/*@}*/


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CEPETRA_VECTOR_H */

