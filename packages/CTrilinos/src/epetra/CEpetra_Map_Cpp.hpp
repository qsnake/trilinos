
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

#ifndef CEPETRA_MAP_CPP_HPP
#define CEPETRA_MAP_CPP_HPP


#include "CTrilinos_enums.h"
#include "Teuchos_RCP.hpp"
#include "Epetra_Map.h"


namespace CEpetra {


using Teuchos::RCP;


/*! get Epetra_Map from non-const table using CT_Epetra_Map_ID */
const RCP<Epetra_Map>
getMap( CT_Epetra_Map_ID_t id );

/*! get Epetra_Map from non-const table using CTrilinos_Universal_ID_t */
const RCP<Epetra_Map>
getMap( CTrilinos_Universal_ID_t id );

/*! get const Epetra_Map from either the const or non-const table
 * using CT_Epetra_Map_ID */
const RCP<const Epetra_Map>
getConstMap( CT_Epetra_Map_ID_t id );

/*! get const Epetra_Map from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const RCP<const Epetra_Map>
getConstMap( CTrilinos_Universal_ID_t id );

/*! store Epetra_Map (owned) in non-const table */
CT_Epetra_Map_ID_t
storeNewMap( Epetra_Map *pobj );

/*! store Epetra_Map in non-const table */
CT_Epetra_Map_ID_t
storeMap( Epetra_Map *pobj );

/*! store const Epetra_Map in const table */
CT_Epetra_Map_ID_t
storeConstMap( const Epetra_Map *pobj );

/* remove Epetra_Map from table using CT_Epetra_Map_ID */
void
removeMap( CT_Epetra_Map_ID_t *id );

/* purge Epetra_Map table */
void
purgeMap(  );

} // namespace CEpetra


#endif // CEPETRA_MAP_CPP_HPP


