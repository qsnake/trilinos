
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

#ifndef CEPETRA_BLOCKMAP_CPP_HPP
#define CEPETRA_BLOCKMAP_CPP_HPP


#include "CTrilinos_enums.h"
#include "Teuchos_RCP.hpp"
#include "Epetra_BlockMap.h"


namespace CEpetra {


using Teuchos::RCP;


/*! get Epetra_BlockMap from non-const table using CT_Epetra_BlockMap_ID */
const RCP<Epetra_BlockMap>
getBlockMap( CT_Epetra_BlockMap_ID_t id );

/*! get Epetra_BlockMap from non-const table using CTrilinos_Universal_ID_t */
const RCP<Epetra_BlockMap>
getBlockMap( CTrilinos_Universal_ID_t id );

/*! get const Epetra_BlockMap from either the const or non-const table
 * using CT_Epetra_BlockMap_ID */
const RCP<const Epetra_BlockMap>
getConstBlockMap( CT_Epetra_BlockMap_ID_t id );

/*! get const Epetra_BlockMap from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const RCP<const Epetra_BlockMap>
getConstBlockMap( CTrilinos_Universal_ID_t id );

/*! store Epetra_BlockMap (owned) in non-const table */
CT_Epetra_BlockMap_ID_t
storeNewBlockMap( Epetra_BlockMap *pobj );

/*! store Epetra_BlockMap in non-const table */
CT_Epetra_BlockMap_ID_t
storeBlockMap( Epetra_BlockMap *pobj );

/*! store const Epetra_BlockMap in const table */
CT_Epetra_BlockMap_ID_t
storeConstBlockMap( const Epetra_BlockMap *pobj );

/* remove Epetra_BlockMap from table using CT_Epetra_BlockMap_ID */
void
removeBlockMap( CT_Epetra_BlockMap_ID_t *id );

/* purge Epetra_BlockMap table */
void
purgeBlockMap(  );

} // namespace CEpetra


#endif // CEPETRA_BLOCKMAP_CPP_HPP


