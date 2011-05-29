
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


/*! @file CTeuchos_any.h
 * @brief Wrappers for Teuchos::any */

/* True C header file! */


#ifndef CTEUCHOS_ANY_H
#define CTEUCHOS_ANY_H


#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif



/*! @name any constructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   Teuchos::any::any()
*/
CT_Teuchos_any_ID_t Teuchos_any_Create (  );

/*! @brief Wrapper for 
   template<typename ValueType> explicit Teuchos::any::any(const ValueType & value)
*/
CT_Teuchos_any_ID_t Teuchos_any_Create_double ( double value );

/*! @brief Wrapper for 
   template<typename ValueType> explicit Teuchos::any::any(const ValueType & value)
*/
CT_Teuchos_any_ID_t Teuchos_any_Create_int ( int value );

/*! @brief Wrapper for 
   Teuchos::any::any(const any & other)
*/
CT_Teuchos_any_ID_t Teuchos_any_Duplicate ( 
  CT_Teuchos_any_ID_t otherID );

/*@}*/

/*! @name any destructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   Teuchos::any::~any()
*/
void Teuchos_any_Destroy ( CT_Teuchos_any_ID_t * selfID );

/*@}*/

/*! @name any member wrappers */
/*@{*/

/*! @brief Wrapper for 
   any & Teuchos::any::swap(any & rhs)
*/
CT_Teuchos_any_ID_t Teuchos_any_swap ( 
  CT_Teuchos_any_ID_t selfID, CT_Teuchos_any_ID_t rhsID );

/*! @brief Wrapper for 
   bool Teuchos::any::empty() const
*/
boolean Teuchos_any_empty ( CT_Teuchos_any_ID_t selfID );

/*! @brief Wrapper for 
   std::string Teuchos::any::typeName() const
*/
const char * Teuchos_any_typeName ( CT_Teuchos_any_ID_t selfID );

/*! @brief Wrapper for 
   bool Teuchos::any::same( const any &other ) const
*/
boolean Teuchos_any_same ( 
  CT_Teuchos_any_ID_t selfID, CT_Teuchos_any_ID_t otherID );

/*@}*/

/*! @name any operator wrappers */
/*@{*/

/*! @brief Wrapper for 
   any & Teuchos::any::operator=(const any & rhs)
*/
void Teuchos_any_Assign ( 
  CT_Teuchos_any_ID_t selfID, CT_Teuchos_any_ID_t rhsID );

/*@}*/


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CTEUCHOS_ANY_H */

