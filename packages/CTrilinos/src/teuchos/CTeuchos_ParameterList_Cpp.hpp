
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

#ifndef CTEUCHOS_PARAMETERLIST_CPP_HPP
#define CTEUCHOS_PARAMETERLIST_CPP_HPP


#include "CTrilinos_enums.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"


namespace CTeuchos {


using Teuchos::RCP;


/*! get Teuchos::ParameterList from non-const table using CT_Teuchos_ParameterList_ID */
const RCP<Teuchos::ParameterList>
getParameterList( CT_Teuchos_ParameterList_ID_t id );

/*! get Teuchos::ParameterList from non-const table using CTrilinos_Universal_ID_t */
const RCP<Teuchos::ParameterList>
getParameterList( CTrilinos_Universal_ID_t id );

/*! get const Teuchos::ParameterList from either the const or non-const table
 * using CT_Teuchos_ParameterList_ID */
const RCP<const Teuchos::ParameterList>
getConstParameterList( CT_Teuchos_ParameterList_ID_t id );

/*! get const Teuchos::ParameterList from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const RCP<const Teuchos::ParameterList>
getConstParameterList( CTrilinos_Universal_ID_t id );

/*! store Teuchos::ParameterList (owned) in non-const table */
CT_Teuchos_ParameterList_ID_t
storeNewParameterList( Teuchos::ParameterList *pobj );

/*! store Teuchos::ParameterList in non-const table */
CT_Teuchos_ParameterList_ID_t
storeParameterList( Teuchos::ParameterList *pobj );

/*! store const Teuchos::ParameterList in const table */
CT_Teuchos_ParameterList_ID_t
storeConstParameterList( const Teuchos::ParameterList *pobj );

/* remove Teuchos::ParameterList from table using CT_Teuchos_ParameterList_ID */
void
removeParameterList( CT_Teuchos_ParameterList_ID_t *id );

/* purge Teuchos::ParameterList table */
void
purgeParameterList(  );

} // namespace CTeuchos


#endif // CTEUCHOS_PARAMETERLIST_CPP_HPP


