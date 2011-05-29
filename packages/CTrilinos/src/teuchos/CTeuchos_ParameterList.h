
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


/*! @file CTeuchos_ParameterList.h
 * @brief Wrappers for Teuchos::ParameterList */

/* True C header file! */


#ifndef CTEUCHOS_PARAMETERLIST_H
#define CTEUCHOS_PARAMETERLIST_H


#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif



/*! @name ID struct conversion functions */
/*@{*/

/*! @brief Changes the ID struct from the universal
   (generalized) struct type to the class-specific one.
*/
CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_Degeneralize ( 
  CTrilinos_Universal_ID_t id );

/*! @brief Changes the ID struct from the class-specific
   struct type to the universal (generalized) one.
*/
CTrilinos_Universal_ID_t Teuchos_ParameterList_Generalize ( 
  CT_Teuchos_ParameterList_ID_t id );

/*@}*/

/*! @name ParameterList constructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   Teuchos::ParameterList::ParameterList()
*/
CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_Create (  );

/*! @brief Wrapper for 
   Teuchos::ParameterList::ParameterList(const std::string &name)
*/
CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_Create_WithName ( 
  const char name[] );

/*! @brief Wrapper for 
   Teuchos::ParameterList::ParameterList(const ParameterList& source)
*/
CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_Create_FromSource ( 
  CT_Teuchos_ParameterList_ID_t sourceID );

/*@}*/

/*! @name ParameterList destructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   virtual Teuchos::ParameterList::~ParameterList()
*/
void Teuchos_ParameterList_Destroy ( 
  CT_Teuchos_ParameterList_ID_t * selfID );

/*@}*/

/*! @name ParameterList member wrappers */
/*@{*/

/*! @brief Wrapper for 
   ParameterList& Teuchos::ParameterList::setName( const std::string &name )
*/
CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_setName ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] );

/*! @brief Wrapper for 
   ParameterList& Teuchos::ParameterList::setParameters(const ParameterList& source)
*/
CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_setParameters ( 
  CT_Teuchos_ParameterList_ID_t selfID, 
  CT_Teuchos_ParameterList_ID_t sourceID );

/*! @brief Wrapper for 
   ParameterList& Teuchos::ParameterList::setParametersNotAlreadySet(const ParameterList& source)
*/
CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_setParametersNotAlreadySet ( 
  CT_Teuchos_ParameterList_ID_t selfID, 
  CT_Teuchos_ParameterList_ID_t sourceID );

/*! @brief Wrapper for 
   ParameterList& Teuchos::ParameterList::disableRecursiveValidation()
*/
CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_disableRecursiveValidation ( 
  CT_Teuchos_ParameterList_ID_t selfID );

/*! @brief Wrapper for 
   template<typename T> ParameterList& Teuchos::ParameterList::set( std::string const& name, T const& value, std::string const& docString = "" ,RCP<const ParameterEntryValidator> const& validator = Teuchos::null )
*/
CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_set_double ( 
  CT_Teuchos_ParameterList_ID_t selfID, char const name[], 
  double value, char const docString[] );

/*! @brief Wrapper for 
   template<typename T> ParameterList& Teuchos::ParameterList::set( std::string const& name, T const& value, std::string const& docString = "" ,RCP<const ParameterEntryValidator> const& validator = Teuchos::null )
*/
CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_set_int ( 
  CT_Teuchos_ParameterList_ID_t selfID, char const name[], 
  int value, char const docString[] );

/*! @brief Wrapper for 
   ParameterList& Teuchos::ParameterList::set( std::string const& name, char value[], std::string const& docString = "" ,RCP<const ParameterEntryValidator> const& validator = Teuchos::null )
*/
CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_set_str ( 
  CT_Teuchos_ParameterList_ID_t selfID, char const name[], 
  char value[], char const docString[] );

/*! @brief Wrapper for 
   ParameterList& Teuchos::ParameterList::set( std::string const& name, ParameterList const& value, std::string const& docString = "" )
*/
CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_set ( 
  CT_Teuchos_ParameterList_ID_t selfID, char const name[], 
  CT_Teuchos_ParameterList_ID_t valueID, char const docString[] );

/*! @brief Wrapper for 
   ParameterList& Teuchos::ParameterList::setEntry(const std::string& name, const ParameterEntry& entry)
*/
CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_setEntry ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[], 
  CT_Teuchos_ParameterEntry_ID_t entryID );

/*! @brief Wrapper for 
   template<typename T> T& Teuchos::ParameterList::get(const std::string& name, T def_value)
*/
double Teuchos_ParameterList_get_double_def ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[], 
  double def_value );

/*! @brief Wrapper for 
   template<typename T> T& Teuchos::ParameterList::get(const std::string& name, T def_value)
*/
int Teuchos_ParameterList_get_int_def ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[], 
  int def_value );

/*! @brief Wrapper for 
   std::string& Teuchos::ParameterList::get(const std::string& name, char def_value[])
*/
const char * Teuchos_ParameterList_get_char_def ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[], 
  char def_value[] );

/*! @brief Wrapper for 
   std::string& Teuchos::ParameterList::get(const std::string& name, const char def_value[])
*/
const char * Teuchos_ParameterList_get_const_char_def ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[], 
  const char def_value[] );

/*! @brief Wrapper for 
   template<typename T> T& Teuchos::ParameterList::get(const std::string& name)
*/
double Teuchos_ParameterList_get_double ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] );

/*! @brief Wrapper for 
   template<typename T> T& Teuchos::ParameterList::get(const std::string& name)
*/
int Teuchos_ParameterList_get_int ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] );

/*! @brief Wrapper for 
   template<typename T> const T& Teuchos::ParameterList::get(const std::string& name) const
*/
double Teuchos_ParameterList_get_double_const ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] );

/*! @brief Wrapper for 
   template<typename T> const T& Teuchos::ParameterList::get(const std::string& name) const
*/
int Teuchos_ParameterList_get_int_const ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] );

/*! @brief Wrapper for 
   template<typename T> inline T* Teuchos::ParameterList::getPtr(const std::string& name)
*/
double * Teuchos_ParameterList_getPtr_double ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] );

/*! @brief Wrapper for 
   template<typename T> inline T* Teuchos::ParameterList::getPtr(const std::string& name)
*/
int * Teuchos_ParameterList_getPtr_int ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] );

/*! @brief Wrapper for 
   template<typename T> inline const T* Teuchos::ParameterList::getPtr(const std::string& name) const
*/
const double * Teuchos_ParameterList_getPtr_double_const ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] );

/*! @brief Wrapper for 
   template<typename T> inline const T* Teuchos::ParameterList::getPtr(const std::string& name) const
*/
const int * Teuchos_ParameterList_getPtr_int_const ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] );

/*! @brief Wrapper for 
   ParameterEntry& Teuchos::ParameterList::getEntry(const std::string& name)
*/
CT_Teuchos_ParameterEntry_ID_t Teuchos_ParameterList_getEntry ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] );

/*! @brief Wrapper for 
   inline const ParameterEntry& Teuchos::ParameterList::getEntry(const std::string& name) const
*/
CT_Teuchos_ParameterEntry_ID_t Teuchos_ParameterList_getEntry_const ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] );

/*! @brief Wrapper for 
   inline ParameterEntry* Teuchos::ParameterList::getEntryPtr(const std::string& name)
*/
CT_Teuchos_ParameterEntry_ID_t Teuchos_ParameterList_getEntryPtr ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] );

/*! @brief Wrapper for 
   inline const ParameterEntry* Teuchos::ParameterList::getEntryPtr(const std::string& name) const
*/
CT_Teuchos_ParameterEntry_ID_t Teuchos_ParameterList_getEntryPtr_const ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] );

/*! @brief Wrapper for 
   bool Teuchos::ParameterList::remove( std::string const& name, bool throwIfNotExists = true )
*/
boolean Teuchos_ParameterList_remove ( 
  CT_Teuchos_ParameterList_ID_t selfID, char const name[], 
  boolean throwIfNotExists );

/*! @brief Wrapper for 
   ParameterList& Teuchos::ParameterList::sublist( const std::string& name, bool mustAlreadyExist = false ,const std::string& docString = "" )
*/
CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_sublist ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[], 
  boolean mustAlreadyExist, const char docString[] );

/*! @brief Wrapper for 
   const ParameterList& Teuchos::ParameterList::sublist(const std::string& name) const
*/
CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_sublist_existing ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] );

/*! @brief Wrapper for 
   const std::string& Teuchos::ParameterList::name() const
*/
const char * Teuchos_ParameterList_name_it ( 
  CT_Teuchos_ParameterList_ID_t selfID );

/*! @brief Wrapper for 
   bool Teuchos::ParameterList::isParameter(const std::string& name) const
*/
boolean Teuchos_ParameterList_isParameter ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] );

/*! @brief Wrapper for 
   bool Teuchos::ParameterList::isSublist(const std::string& name) const
*/
boolean Teuchos_ParameterList_isSublist ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] );

/*! @brief Wrapper for 
   template<typename T> bool Teuchos::ParameterList::isType(const std::string& name) const
*/
boolean Teuchos_ParameterList_isType_double ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] );

/*! @brief Wrapper for 
   template<typename T> bool Teuchos::ParameterList::isType(const std::string& name) const
*/
boolean Teuchos_ParameterList_isType_int ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] );

/*! @brief Wrapper for 
   template<typename T> bool Teuchos::ParameterList::isType(const std::string& name, T* ptr) const
*/
boolean Teuchos_ParameterList_isType_double_type ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[], 
  double * ptr );

/*! @brief Wrapper for 
   template<typename T> bool Teuchos::ParameterList::isType(const std::string& name, T* ptr) const
*/
boolean Teuchos_ParameterList_isType_int_type ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[], 
  int * ptr );

/*! @brief Wrapper for 
   std::string Teuchos::ParameterList::currentParametersString() const
*/
const char * Teuchos_ParameterList_currentParametersString ( 
  CT_Teuchos_ParameterList_ID_t selfID );

/*! @brief Wrapper for 
   void Teuchos::ParameterList::validateParameters( ParameterList const& validParamList, int const depth = 1000, EValidateUsed const validateUsed = VALIDATE_USED_ENABLED, EValidateDefaults const validateDefaults = VALIDATE_DEFAULTS_ENABLED ) const
*/
void Teuchos_ParameterList_validateParameters ( 
  CT_Teuchos_ParameterList_ID_t selfID, 
  CT_Teuchos_ParameterList_ID_t validParamListID, int const depth, 
  const CT_EValidateUsed_E_t validateUsed, 
  const CT_EValidateDefaults_E_t validateDefaults );

/*! @brief Wrapper for 
   void Teuchos::ParameterList::validateParametersAndSetDefaults( ParameterList const& validParamList, int const depth = 1000 )
*/
void Teuchos_ParameterList_validateParametersAndSetDefaults ( 
  CT_Teuchos_ParameterList_ID_t selfID, 
  CT_Teuchos_ParameterList_ID_t validParamListID, int const depth );

/*@}*/

/*! @name ParameterList operator wrappers */
/*@{*/

/*! @brief Wrapper for 
   ParameterList& Teuchos::ParameterList::operator=(const ParameterList& source)
*/
void Teuchos_ParameterList_Assign ( 
  CT_Teuchos_ParameterList_ID_t selfID, 
  CT_Teuchos_ParameterList_ID_t sourceID );

/*@}*/


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CTEUCHOS_PARAMETERLIST_H */

