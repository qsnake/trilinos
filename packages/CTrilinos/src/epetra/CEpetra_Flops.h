
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


/*! @file CEpetra_Flops.h
 * @brief Wrappers for Epetra_Flops */

/* True C header file! */


#ifndef CEPETRA_FLOPS_H
#define CEPETRA_FLOPS_H


#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif



/*! @name ID struct conversion functions */
/*@{*/

/*! @brief Changes the ID struct from the universal
   (generalized) struct type to the class-specific one.
*/
CT_Epetra_Flops_ID_t Epetra_Flops_Degeneralize ( 
  CTrilinos_Universal_ID_t id );

/*! @brief Changes the ID struct from the class-specific
   struct type to the universal (generalized) one.
*/
CTrilinos_Universal_ID_t Epetra_Flops_Generalize ( 
  CT_Epetra_Flops_ID_t id );

/*@}*/

/*! @name Epetra_Flops constructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   Epetra_Flops::Epetra_Flops(void)
*/
CT_Epetra_Flops_ID_t Epetra_Flops_Create (  );

/*! @brief Wrapper for 
   Epetra_Flops::Epetra_Flops(const Epetra_Flops& Flops_in)
*/
CT_Epetra_Flops_ID_t Epetra_Flops_Duplicate ( 
  CT_Epetra_Flops_ID_t Flops_inID );

/*@}*/

/*! @name Epetra_Flops destructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   virtual Epetra_Flops::~Epetra_Flops(void)
*/
void Epetra_Flops_Destroy ( CT_Epetra_Flops_ID_t * selfID );

/*@}*/

/*! @name Epetra_Flops member wrappers */
/*@{*/

/*! @brief Wrapper for 
   double Epetra_Flops::Flops() const
*/
double Epetra_Flops_Flops ( CT_Epetra_Flops_ID_t selfID );

/*! @brief Wrapper for 
   void Epetra_Flops::ResetFlops()
*/
void Epetra_Flops_ResetFlops ( CT_Epetra_Flops_ID_t selfID );

/*@}*/

/*! @name Epetra_Flops operator wrappers */
/*@{*/

/*! @brief Wrapper for 
   Epetra_Flops& Epetra_Flops::operator=(const Epetra_Flops& src)
*/
void Epetra_Flops_Assign ( 
  CT_Epetra_Flops_ID_t selfID, CT_Epetra_Flops_ID_t srcID );

/*@}*/


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CEPETRA_FLOPS_H */

