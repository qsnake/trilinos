
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

#ifndef CTEUCHOS_ANY_CPP_HPP
#define CTEUCHOS_ANY_CPP_HPP


#include "CTrilinos_enums.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_any.hpp"


namespace CTeuchos {


using Teuchos::RCP;


/*! get Teuchos::any from non-const table using CT_Teuchos_any_ID */
const RCP<Teuchos::any>
getany( CT_Teuchos_any_ID_t id );

/*! get Teuchos::any from non-const table using CTrilinos_Universal_ID_t */
const RCP<Teuchos::any>
getany( CTrilinos_Universal_ID_t id );

/*! get const Teuchos::any from either the const or non-const table
 * using CT_Teuchos_any_ID */
const RCP<const Teuchos::any>
getConstany( CT_Teuchos_any_ID_t id );

/*! get const Teuchos::any from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const RCP<const Teuchos::any>
getConstany( CTrilinos_Universal_ID_t id );

/*! store Teuchos::any (owned) in non-const table */
CT_Teuchos_any_ID_t
storeNewany( Teuchos::any *pobj );

/*! store Teuchos::any in non-const table */
CT_Teuchos_any_ID_t
storeany( Teuchos::any *pobj );

/*! store const Teuchos::any in const table */
CT_Teuchos_any_ID_t
storeConstany( const Teuchos::any *pobj );

/* remove Teuchos::any from table using CT_Teuchos_any_ID */
void
removeany( CT_Teuchos_any_ID_t *id );

/* purge Teuchos::any table */
void
purgeany(  );

} // namespace CTeuchos


#endif // CTEUCHOS_ANY_CPP_HPP


