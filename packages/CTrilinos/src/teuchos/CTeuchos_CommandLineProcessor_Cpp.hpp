
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

#ifndef CTEUCHOS_COMMANDLINEPROCESSOR_CPP_HPP
#define CTEUCHOS_COMMANDLINEPROCESSOR_CPP_HPP


#include "CTrilinos_enums.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_CommandLineProcessor.hpp"


namespace CTeuchos {


using Teuchos::RCP;


/*! get Teuchos::CommandLineProcessor from non-const table using CT_Teuchos_CommandLineProcessor_ID */
const RCP<Teuchos::CommandLineProcessor>
getCommandLineProcessor( CT_Teuchos_CommandLineProcessor_ID_t id );

/*! get Teuchos::CommandLineProcessor from non-const table using CTrilinos_Universal_ID_t */
const RCP<Teuchos::CommandLineProcessor>
getCommandLineProcessor( CTrilinos_Universal_ID_t id );

/*! get const Teuchos::CommandLineProcessor from either the const or non-const table
 * using CT_Teuchos_CommandLineProcessor_ID */
const RCP<const Teuchos::CommandLineProcessor>
getConstCommandLineProcessor( CT_Teuchos_CommandLineProcessor_ID_t id );

/*! get const Teuchos::CommandLineProcessor from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const RCP<const Teuchos::CommandLineProcessor>
getConstCommandLineProcessor( CTrilinos_Universal_ID_t id );

/*! store Teuchos::CommandLineProcessor (owned) in non-const table */
CT_Teuchos_CommandLineProcessor_ID_t
storeNewCommandLineProcessor( Teuchos::CommandLineProcessor *pobj );

/*! store Teuchos::CommandLineProcessor in non-const table */
CT_Teuchos_CommandLineProcessor_ID_t
storeCommandLineProcessor( Teuchos::CommandLineProcessor *pobj );

/*! store const Teuchos::CommandLineProcessor in const table */
CT_Teuchos_CommandLineProcessor_ID_t
storeConstCommandLineProcessor( const Teuchos::CommandLineProcessor *pobj );

/* remove Teuchos::CommandLineProcessor from table using CT_Teuchos_CommandLineProcessor_ID */
void
removeCommandLineProcessor( CT_Teuchos_CommandLineProcessor_ID_t *id );

/* purge Teuchos::CommandLineProcessor table */
void
purgeCommandLineProcessor(  );

} // namespace CTeuchos


#endif // CTEUCHOS_COMMANDLINEPROCESSOR_CPP_HPP


