
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

#ifndef CTEUCHOS_PARAMETERENTRY_CPP_HPP
#define CTEUCHOS_PARAMETERENTRY_CPP_HPP


#include "CTrilinos_enums.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterEntry.hpp"


namespace CTeuchos {


using Teuchos::RCP;


/*! get Teuchos::ParameterEntry from non-const table using CT_Teuchos_ParameterEntry_ID */
const RCP<Teuchos::ParameterEntry>
getParameterEntry( CT_Teuchos_ParameterEntry_ID_t id );

/*! get Teuchos::ParameterEntry from non-const table using CTrilinos_Universal_ID_t */
const RCP<Teuchos::ParameterEntry>
getParameterEntry( CTrilinos_Universal_ID_t id );

/*! get const Teuchos::ParameterEntry from either the const or non-const table
 * using CT_Teuchos_ParameterEntry_ID */
const RCP<const Teuchos::ParameterEntry>
getConstParameterEntry( CT_Teuchos_ParameterEntry_ID_t id );

/*! get const Teuchos::ParameterEntry from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const RCP<const Teuchos::ParameterEntry>
getConstParameterEntry( CTrilinos_Universal_ID_t id );

/*! store Teuchos::ParameterEntry (owned) in non-const table */
CT_Teuchos_ParameterEntry_ID_t
storeNewParameterEntry( Teuchos::ParameterEntry *pobj );

/*! store Teuchos::ParameterEntry in non-const table */
CT_Teuchos_ParameterEntry_ID_t
storeParameterEntry( Teuchos::ParameterEntry *pobj );

/*! store const Teuchos::ParameterEntry in const table */
CT_Teuchos_ParameterEntry_ID_t
storeConstParameterEntry( const Teuchos::ParameterEntry *pobj );

/* remove Teuchos::ParameterEntry from table using CT_Teuchos_ParameterEntry_ID */
void
removeParameterEntry( CT_Teuchos_ParameterEntry_ID_t *id );

/* purge Teuchos::ParameterEntry table */
void
purgeParameterEntry(  );

} // namespace CTeuchos


#endif // CTEUCHOS_PARAMETERENTRY_CPP_HPP


