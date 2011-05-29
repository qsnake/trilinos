
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



/*! @file CIfpack_Preconditioner.h
 * @brief Wrappers for Ifpack_Preconditioner */

/* True C header file! */


#ifndef CIFPACK_PRECONDITIONER_H
#define CIFPACK_PRECONDITIONER_H


#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif



/*! @name ID struct conversion functions */
/*@{*/

/*! @brief Changes the ID struct from the universal
   (generalized) struct type to the class-specific one.
*/
CT_Ifpack_Preconditioner_ID_t Ifpack_Preconditioner_Degeneralize ( 
  CTrilinos_Universal_ID_t id );

/*! @brief Changes the ID struct from the class-specific
   struct type to the universal (generalized) one.
*/
CTrilinos_Universal_ID_t Ifpack_Preconditioner_Generalize ( 
  CT_Ifpack_Preconditioner_ID_t id );

/*@}*/

/*! @name Ifpack_Preconditioner member wrappers */
/*@{*/

/*! @brief Wrapper for 
   virtual int Ifpack_Preconditioner::SetParameters(Teuchos::ParameterList& List) = 0
*/
int Ifpack_Preconditioner_SetParameters ( 
  CT_Ifpack_Preconditioner_ID_t selfID, 
  CT_Teuchos_ParameterList_ID_t ListID );

/*! @brief Wrapper for 
   virtual int Ifpack_Preconditioner::Initialize() = 0
*/
int Ifpack_Preconditioner_Initialize ( 
  CT_Ifpack_Preconditioner_ID_t selfID );

/*! @brief Wrapper for 
   virtual bool Ifpack_Preconditioner::IsInitialized() const = 0
*/
boolean Ifpack_Preconditioner_IsInitialized ( 
  CT_Ifpack_Preconditioner_ID_t selfID );

/*! @brief Wrapper for 
   virtual int Ifpack_Preconditioner::Compute() = 0
*/
int Ifpack_Preconditioner_Compute ( 
  CT_Ifpack_Preconditioner_ID_t selfID );

/*! @brief Wrapper for 
   virtual bool Ifpack_Preconditioner::IsComputed() const = 0
*/
boolean Ifpack_Preconditioner_IsComputed ( 
  CT_Ifpack_Preconditioner_ID_t selfID );

/*! @brief Wrapper for 
   virtual double Ifpack_Preconditioner::Condest() const = 0
*/
double Ifpack_Preconditioner_Condest ( 
  CT_Ifpack_Preconditioner_ID_t selfID );

/*! @brief Wrapper for 
   virtual int Ifpack_Preconditioner::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const = 0
*/
int Ifpack_Preconditioner_ApplyInverse ( 
  CT_Ifpack_Preconditioner_ID_t selfID, 
  CT_Epetra_MultiVector_ID_t XID, CT_Epetra_MultiVector_ID_t YID );

/*! @brief Wrapper for 
   virtual const Epetra_RowMatrix& Ifpack_Preconditioner::Matrix() const = 0
*/
CT_Epetra_RowMatrix_ID_t Ifpack_Preconditioner_Matrix ( 
  CT_Ifpack_Preconditioner_ID_t selfID );

/*! @brief Wrapper for 
   virtual int Ifpack_Preconditioner::NumInitialize() const = 0
*/
int Ifpack_Preconditioner_NumInitialize ( 
  CT_Ifpack_Preconditioner_ID_t selfID );

/*! @brief Wrapper for 
   virtual int Ifpack_Preconditioner::NumCompute() const = 0
*/
int Ifpack_Preconditioner_NumCompute ( 
  CT_Ifpack_Preconditioner_ID_t selfID );

/*! @brief Wrapper for 
   virtual int Ifpack_Preconditioner::NumApplyInverse() const = 0
*/
int Ifpack_Preconditioner_NumApplyInverse ( 
  CT_Ifpack_Preconditioner_ID_t selfID );

/*! @brief Wrapper for 
   virtual double Ifpack_Preconditioner::InitializeTime() const = 0
*/
double Ifpack_Preconditioner_InitializeTime ( 
  CT_Ifpack_Preconditioner_ID_t selfID );

/*! @brief Wrapper for 
   virtual double Ifpack_Preconditioner::ComputeTime() const = 0
*/
double Ifpack_Preconditioner_ComputeTime ( 
  CT_Ifpack_Preconditioner_ID_t selfID );

/*! @brief Wrapper for 
   virtual double Ifpack_Preconditioner::ApplyInverseTime() const = 0
*/
double Ifpack_Preconditioner_ApplyInverseTime ( 
  CT_Ifpack_Preconditioner_ID_t selfID );

/*! @brief Wrapper for 
   virtual double Ifpack_Preconditioner::InitializeFlops() const = 0
*/
double Ifpack_Preconditioner_InitializeFlops ( 
  CT_Ifpack_Preconditioner_ID_t selfID );

/*! @brief Wrapper for 
   virtual double Ifpack_Preconditioner::ComputeFlops() const = 0
*/
double Ifpack_Preconditioner_ComputeFlops ( 
  CT_Ifpack_Preconditioner_ID_t selfID );

/*! @brief Wrapper for 
   virtual double Ifpack_Preconditioner::ApplyInverseFlops() const = 0
*/
double Ifpack_Preconditioner_ApplyInverseFlops ( 
  CT_Ifpack_Preconditioner_ID_t selfID );

/*@}*/


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CIFPACK_PRECONDITIONER_H */

#endif /* HAVE_CTRILINOS_IFPACK */


