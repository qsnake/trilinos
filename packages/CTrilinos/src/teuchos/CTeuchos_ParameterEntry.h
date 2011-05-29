
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


/*! @file CTeuchos_ParameterEntry.h
 * @brief Wrappers for Teuchos::ParameterEntry */

/* True C header file! */


#ifndef CTEUCHOS_PARAMETERENTRY_H
#define CTEUCHOS_PARAMETERENTRY_H


#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif



/*! @name ParameterEntry constructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   Teuchos::ParameterEntry::ParameterEntry()
*/
CT_Teuchos_ParameterEntry_ID_t Teuchos_ParameterEntry_Create (  );

/*! @brief Wrapper for 
   Teuchos::ParameterEntry::ParameterEntry(const ParameterEntry& source)
*/
CT_Teuchos_ParameterEntry_ID_t Teuchos_ParameterEntry_Duplicate ( 
  CT_Teuchos_ParameterEntry_ID_t sourceID );

/*@}*/

/*! @name ParameterEntry destructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   Teuchos::ParameterEntry::~ParameterEntry()
*/
void Teuchos_ParameterEntry_Destroy ( 
  CT_Teuchos_ParameterEntry_ID_t * selfID );

/*@}*/

/*! @name ParameterEntry member wrappers */
/*@{*/

/*! @brief Wrapper for 
   void Teuchos::ParameterEntry::setAnyValue( const any &value, bool isDefault = false )
*/
void Teuchos_ParameterEntry_setAnyValue ( 
  CT_Teuchos_ParameterEntry_ID_t selfID, 
  CT_Teuchos_any_ID_t valueID, boolean isDefault );

/*! @brief Wrapper for 
   void Teuchos::ParameterEntry::setDocString(const std::string &docString)
*/
void Teuchos_ParameterEntry_setDocString ( 
  CT_Teuchos_ParameterEntry_ID_t selfID, const char docString[] );

/*! @brief Wrapper for 
   ParameterList& Teuchos::ParameterEntry::setList( bool isDefault = false, const std::string &docString = "" )
*/
CT_Teuchos_ParameterList_ID_t Teuchos_ParameterEntry_setList ( 
  CT_Teuchos_ParameterEntry_ID_t selfID, boolean isDefault, 
  const char docString[] );

/*! @brief Wrapper for 
   template<typename T> inline T& Teuchos::ParameterEntry::getValue(T *ptr) const
*/
double Teuchos_ParameterEntry_getValue_double ( 
  CT_Teuchos_ParameterEntry_ID_t selfID, double * ptr );

/*! @brief Wrapper for 
   template<typename T> inline T& Teuchos::ParameterEntry::getValue(T *ptr) const
*/
int Teuchos_ParameterEntry_getValue_int ( 
  CT_Teuchos_ParameterEntry_ID_t selfID, int * ptr );

/*! @brief Wrapper for 
   inline any& Teuchos::ParameterEntry::getAny(bool activeQry = true)
*/
CT_Teuchos_any_ID_t Teuchos_ParameterEntry_getAny ( 
  CT_Teuchos_ParameterEntry_ID_t selfID, boolean activeQry );

/*! @brief Wrapper for 
   inline const any& Teuchos::ParameterEntry::getAny(bool activeQry = true) const
*/
CT_Teuchos_any_ID_t Teuchos_ParameterEntry_getAny_const ( 
  CT_Teuchos_ParameterEntry_ID_t selfID, boolean activeQry );

/*! @brief Wrapper for 
   inline bool Teuchos::ParameterEntry::isUsed() const
*/
boolean Teuchos_ParameterEntry_isUsed ( 
  CT_Teuchos_ParameterEntry_ID_t selfID );

/*! @brief Wrapper for 
   bool Teuchos::ParameterEntry::isList() const
*/
boolean Teuchos_ParameterEntry_isList ( 
  CT_Teuchos_ParameterEntry_ID_t selfID );

/*! @brief Wrapper for 
   template <typename T> inline bool Teuchos::ParameterEntry::isType() const
*/
boolean Teuchos_ParameterEntry_isType_double ( 
  CT_Teuchos_ParameterEntry_ID_t selfID );

/*! @brief Wrapper for 
   template <typename T> inline bool Teuchos::ParameterEntry::isType() const
*/
boolean Teuchos_ParameterEntry_isType_int ( 
  CT_Teuchos_ParameterEntry_ID_t selfID );

/*! @brief Wrapper for 
   inline bool Teuchos::ParameterEntry::isDefault() const
*/
boolean Teuchos_ParameterEntry_isDefault ( 
  CT_Teuchos_ParameterEntry_ID_t selfID );

/*! @brief Wrapper for 
   inline std::string Teuchos::ParameterEntry::docString() const
*/
const char * Teuchos_ParameterEntry_docString ( 
  CT_Teuchos_ParameterEntry_ID_t selfID );

/*@}*/

/*! @name ParameterEntry operator wrappers */
/*@{*/

/*! @brief Wrapper for 
   ParameterEntry& Teuchos::ParameterEntry::operator=(const ParameterEntry& source)
*/
void Teuchos_ParameterEntry_Assign ( 
  CT_Teuchos_ParameterEntry_ID_t selfID, 
  CT_Teuchos_ParameterEntry_ID_t sourceID );

/*@}*/


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CTEUCHOS_PARAMETERENTRY_H */

