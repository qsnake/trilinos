
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


/*! @file CEpetra_LinearProblem.h
 * @brief Wrappers for Epetra_LinearProblem */

/* True C header file! */


#ifndef CEPETRA_LINEARPROBLEM_H
#define CEPETRA_LINEARPROBLEM_H


#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif



/*! @name ID struct conversion functions */
/*@{*/

/*! @brief Changes the ID struct from the universal
   (generalized) struct type to the class-specific one.
*/
CT_Epetra_LinearProblem_ID_t Epetra_LinearProblem_Degeneralize ( 
  CTrilinos_Universal_ID_t id );

/*! @brief Changes the ID struct from the class-specific
   struct type to the universal (generalized) one.
*/
CTrilinos_Universal_ID_t Epetra_LinearProblem_Generalize ( 
  CT_Epetra_LinearProblem_ID_t id );

/*@}*/

/*! @name Epetra_LinearProblem constructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   Epetra_LinearProblem::Epetra_LinearProblem(void)
*/
CT_Epetra_LinearProblem_ID_t Epetra_LinearProblem_Create (  );

/*! @brief Wrapper for 
   Epetra_LinearProblem::Epetra_LinearProblem(Epetra_RowMatrix * A, Epetra_MultiVector * X, Epetra_MultiVector * B)
*/
CT_Epetra_LinearProblem_ID_t Epetra_LinearProblem_Create_FromMatrix ( 
  CT_Epetra_RowMatrix_ID_t AID, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t BID );

/*! @brief Wrapper for 
   Epetra_LinearProblem::Epetra_LinearProblem(Epetra_Operator * A, Epetra_MultiVector * X, Epetra_MultiVector * B)
*/
CT_Epetra_LinearProblem_ID_t Epetra_LinearProblem_Create_FromOperator ( 
  CT_Epetra_Operator_ID_t AID, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t BID );

/*! @brief Wrapper for 
   Epetra_LinearProblem::Epetra_LinearProblem(const Epetra_LinearProblem& Problem)
*/
CT_Epetra_LinearProblem_ID_t Epetra_LinearProblem_Duplicate ( 
  CT_Epetra_LinearProblem_ID_t ProblemID );

/*@}*/

/*! @name Epetra_LinearProblem destructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   virtual Epetra_LinearProblem::~Epetra_LinearProblem(void)
*/
void Epetra_LinearProblem_Destroy ( 
  CT_Epetra_LinearProblem_ID_t * selfID );

/*@}*/

/*! @name Epetra_LinearProblem member wrappers */
/*@{*/

/*! @brief Wrapper for 
   int Epetra_LinearProblem::CheckInput() const
*/
int Epetra_LinearProblem_CheckInput ( 
  CT_Epetra_LinearProblem_ID_t selfID );

/*! @brief Wrapper for 
   void Epetra_LinearProblem::AssertSymmetric()
*/
void Epetra_LinearProblem_AssertSymmetric ( 
  CT_Epetra_LinearProblem_ID_t selfID );

/*! @brief Wrapper for 
   void Epetra_LinearProblem::SetPDL(ProblemDifficultyLevel PDL)
*/
void Epetra_LinearProblem_SetPDL ( 
  CT_Epetra_LinearProblem_ID_t selfID, 
  CT_ProblemDifficultyLevel_E_t PDL );

/*! @brief Wrapper for 
   void Epetra_LinearProblem::SetOperator(Epetra_RowMatrix * A)
*/
void Epetra_LinearProblem_SetOperator_Matrix ( 
  CT_Epetra_LinearProblem_ID_t selfID, CT_Epetra_RowMatrix_ID_t AID );

/*! @brief Wrapper for 
   void Epetra_LinearProblem::SetOperator(Epetra_Operator * A)
*/
void Epetra_LinearProblem_SetOperator ( 
  CT_Epetra_LinearProblem_ID_t selfID, CT_Epetra_Operator_ID_t AID );

/*! @brief Wrapper for 
   void Epetra_LinearProblem::SetLHS(Epetra_MultiVector * X)
*/
void Epetra_LinearProblem_SetLHS ( 
  CT_Epetra_LinearProblem_ID_t selfID, 
  CT_Epetra_MultiVector_ID_t XID );

/*! @brief Wrapper for 
   void Epetra_LinearProblem::SetRHS(Epetra_MultiVector * B)
*/
void Epetra_LinearProblem_SetRHS ( 
  CT_Epetra_LinearProblem_ID_t selfID, 
  CT_Epetra_MultiVector_ID_t BID );

/*! @brief Wrapper for 
   int Epetra_LinearProblem::LeftScale(const Epetra_Vector & D)
*/
int Epetra_LinearProblem_LeftScale ( 
  CT_Epetra_LinearProblem_ID_t selfID, CT_Epetra_Vector_ID_t DID );

/*! @brief Wrapper for 
   int Epetra_LinearProblem::RightScale(const Epetra_Vector & D)
*/
int Epetra_LinearProblem_RightScale ( 
  CT_Epetra_LinearProblem_ID_t selfID, CT_Epetra_Vector_ID_t DID );

/*! @brief Wrapper for 
   Epetra_Operator * Epetra_LinearProblem::GetOperator() const
*/
CT_Epetra_Operator_ID_t Epetra_LinearProblem_GetOperator ( 
  CT_Epetra_LinearProblem_ID_t selfID );

/*! @brief Wrapper for 
   Epetra_RowMatrix * Epetra_LinearProblem::GetMatrix() const
*/
CT_Epetra_RowMatrix_ID_t Epetra_LinearProblem_GetMatrix ( 
  CT_Epetra_LinearProblem_ID_t selfID );

/*! @brief Wrapper for 
   Epetra_MultiVector * Epetra_LinearProblem::GetLHS() const
*/
CT_Epetra_MultiVector_ID_t Epetra_LinearProblem_GetLHS ( 
  CT_Epetra_LinearProblem_ID_t selfID );

/*! @brief Wrapper for 
   Epetra_MultiVector * Epetra_LinearProblem::GetRHS() const
*/
CT_Epetra_MultiVector_ID_t Epetra_LinearProblem_GetRHS ( 
  CT_Epetra_LinearProblem_ID_t selfID );

/*! @brief Wrapper for 
   ProblemDifficultyLevel Epetra_LinearProblem::GetPDL() const
*/
CT_ProblemDifficultyLevel_E_t Epetra_LinearProblem_GetPDL ( 
  CT_Epetra_LinearProblem_ID_t selfID );

/*! @brief Wrapper for 
   bool Epetra_LinearProblem::IsOperatorSymmetric() const
*/
boolean Epetra_LinearProblem_IsOperatorSymmetric ( 
  CT_Epetra_LinearProblem_ID_t selfID );

/*@}*/


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CEPETRA_LINEARPROBLEM_H */

