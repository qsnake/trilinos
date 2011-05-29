
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

#include "CTrilinos_enums.h"
#include "CTeuchos_CommandLineProcessor.h"
#include "CTeuchos_CommandLineProcessor_Cpp.hpp"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_Table.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


/* table to hold objects of type Teuchos::CommandLineProcessor */
Table<Teuchos::CommandLineProcessor>& tableOfCommandLineProcessors()
{
    static Table<Teuchos::CommandLineProcessor> loc_tableOfCommandLineProcessors(CT_Teuchos_CommandLineProcessor_ID);
    return loc_tableOfCommandLineProcessors;
}


} // namespace


//
// Definitions from CTeuchos_CommandLineProcessor.h
//


extern "C" {


CT_Teuchos_CommandLineProcessor_ID_t Teuchos_CommandLineProcessor_Create ( 
  boolean throwExceptions, boolean recogniseAllOptions, 
  boolean addOutputSetupOptions )
{
    return CTeuchos::storeNewCommandLineProcessor(
        new Teuchos::CommandLineProcessor(((throwExceptions) != 
        FALSE ? true : false), ((recogniseAllOptions) != 
        FALSE ? true : false), ((addOutputSetupOptions) != 
        FALSE ? true : false)));
}

void Teuchos_CommandLineProcessor_Destroy ( 
  CT_Teuchos_CommandLineProcessor_ID_t * selfID )
{
    CTeuchos::removeCommandLineProcessor(selfID);
}

void Teuchos_CommandLineProcessor_throwExceptions_set ( 
  CT_Teuchos_CommandLineProcessor_ID_t selfID, 
  const boolean throwExceptions )
{
    CTeuchos::getCommandLineProcessor(selfID)->throwExceptions(
        ((throwExceptions) != FALSE ? true : false));
}

boolean Teuchos_CommandLineProcessor_throwExceptions_get ( 
  CT_Teuchos_CommandLineProcessor_ID_t selfID )
{
    return ((CTeuchos::getConstCommandLineProcessor(
        selfID)->throwExceptions()) ? TRUE : FALSE);
}

void Teuchos_CommandLineProcessor_recogniseAllOptions_set ( 
  CT_Teuchos_CommandLineProcessor_ID_t selfID, 
  const boolean recogniseAllOptions )
{
    CTeuchos::getCommandLineProcessor(selfID)->recogniseAllOptions(
        ((recogniseAllOptions) != FALSE ? true : false));
}

boolean Teuchos_CommandLineProcessor_recogniseAllOptions_get ( 
  CT_Teuchos_CommandLineProcessor_ID_t selfID )
{
    return ((CTeuchos::getConstCommandLineProcessor(
        selfID)->recogniseAllOptions()) ? TRUE : FALSE);
}

void Teuchos_CommandLineProcessor_addOutputSetupOptions_set ( 
  CT_Teuchos_CommandLineProcessor_ID_t selfID, 
  const boolean addOutputSetupOptions )
{
    CTeuchos::getCommandLineProcessor(selfID)->addOutputSetupOptions(
        ((addOutputSetupOptions) != FALSE ? true : false));
}

boolean Teuchos_CommandLineProcessor_addOutputSetupOptions_get ( 
  CT_Teuchos_CommandLineProcessor_ID_t selfID )
{
    return ((CTeuchos::getConstCommandLineProcessor(
        selfID)->addOutputSetupOptions()) ? TRUE : FALSE);
}

void Teuchos_CommandLineProcessor_setDocString ( 
  CT_Teuchos_CommandLineProcessor_ID_t selfID, 
  const char doc_string[] )
{
    CTeuchos::getCommandLineProcessor(selfID)->setDocString(doc_string);
}

void Teuchos_CommandLineProcessor_setOption_bool ( 
  CT_Teuchos_CommandLineProcessor_ID_t selfID, 
  const char option_true[], const char option_false[], 
  boolean * option_val, const char documentation[] )
{
    bool *tmp_option_val = NULL;
    CTrilinos::pass_bool_in(option_val, tmp_option_val);
    CTeuchos::getCommandLineProcessor(selfID)->setOption(option_true, 
        option_false, tmp_option_val, documentation);
    CTrilinos::pass_bool_out(tmp_option_val, option_val);
    delete tmp_option_val;
}

void Teuchos_CommandLineProcessor_setOption_int ( 
  CT_Teuchos_CommandLineProcessor_ID_t selfID, 
  const char option_name[], int * option_val, 
  const char documentation[], const boolean required )
{
    CTeuchos::getCommandLineProcessor(selfID)->setOption(option_name, 
        option_val, documentation, ((required) != FALSE ? true : false));
}

void Teuchos_CommandLineProcessor_setOption_double ( 
  CT_Teuchos_CommandLineProcessor_ID_t selfID, 
  const char option_name[], double * option_val, 
  const char documentation[], const boolean required )
{
    CTeuchos::getCommandLineProcessor(selfID)->setOption(option_name, 
        option_val, documentation, ((required) != FALSE ? true : false));
}

void Teuchos_CommandLineProcessor_setOption_str ( 
  CT_Teuchos_CommandLineProcessor_ID_t selfID, 
  const char option_name[], char * option_val[], 
  const char documentation[], const boolean required )
{
    std::string *tmp_option_val = NULL;
    CTrilinos::pass_string_in(option_val, tmp_option_val);
    CTeuchos::getCommandLineProcessor(selfID)->setOption(option_name, 
        tmp_option_val, documentation, ((required) != FALSE ? true : false));
    CTrilinos::pass_string_out(tmp_option_val, option_val);
    delete tmp_option_val;
}


} // extern "C"


//
// Definitions from CTeuchos_CommandLineProcessor_Cpp.hpp
//


/* get Teuchos::CommandLineProcessor from non-const table using CT_Teuchos_CommandLineProcessor_ID */
const Teuchos::RCP<Teuchos::CommandLineProcessor>
CTeuchos::getCommandLineProcessor( CT_Teuchos_CommandLineProcessor_ID_t id )
{
    return tableOfCommandLineProcessors().get(
        CTrilinos::abstractType<CT_Teuchos_CommandLineProcessor_ID_t>(id));
}

/* get Teuchos::CommandLineProcessor from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Teuchos::CommandLineProcessor>
CTeuchos::getCommandLineProcessor( CTrilinos_Universal_ID_t id )
{
    return tableOfCommandLineProcessors().get(id);
}

/* get const Teuchos::CommandLineProcessor from either the const or non-const table
 * using CT_Teuchos_CommandLineProcessor_ID */
const Teuchos::RCP<const Teuchos::CommandLineProcessor>
CTeuchos::getConstCommandLineProcessor( CT_Teuchos_CommandLineProcessor_ID_t id )
{
    return tableOfCommandLineProcessors().getConst(
        CTrilinos::abstractType<CT_Teuchos_CommandLineProcessor_ID_t>(id));
}

/* get const Teuchos::CommandLineProcessor from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Teuchos::CommandLineProcessor>
CTeuchos::getConstCommandLineProcessor( CTrilinos_Universal_ID_t id )
{
    return tableOfCommandLineProcessors().getConst(id);
}

/* store Teuchos::CommandLineProcessor (owned) in non-const table */
CT_Teuchos_CommandLineProcessor_ID_t
CTeuchos::storeNewCommandLineProcessor( Teuchos::CommandLineProcessor *pobj )
{
    return CTrilinos::concreteType<CT_Teuchos_CommandLineProcessor_ID_t>(
        tableOfCommandLineProcessors().store(pobj, true));
}

/* store Teuchos::CommandLineProcessor in non-const table */
CT_Teuchos_CommandLineProcessor_ID_t
CTeuchos::storeCommandLineProcessor( Teuchos::CommandLineProcessor *pobj )
{
    return CTrilinos::concreteType<CT_Teuchos_CommandLineProcessor_ID_t>(
        tableOfCommandLineProcessors().store(pobj, false));
}

/* store const Teuchos::CommandLineProcessor in const table */
CT_Teuchos_CommandLineProcessor_ID_t
CTeuchos::storeConstCommandLineProcessor( const Teuchos::CommandLineProcessor *pobj )
{
    return CTrilinos::concreteType<CT_Teuchos_CommandLineProcessor_ID_t>(
        tableOfCommandLineProcessors().store(pobj, false));
}

/* remove Teuchos::CommandLineProcessor from table using CT_Teuchos_CommandLineProcessor_ID */
void
CTeuchos::removeCommandLineProcessor( CT_Teuchos_CommandLineProcessor_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Teuchos_CommandLineProcessor_ID_t>(*id);
    tableOfCommandLineProcessors().remove(&aid);
    *id = CTrilinos::concreteType<CT_Teuchos_CommandLineProcessor_ID_t>(aid);
}

/* purge Teuchos::CommandLineProcessor table */
void
CTeuchos::purgeCommandLineProcessor(  )
{
    tableOfCommandLineProcessors().purge();
}



