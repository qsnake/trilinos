
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

#ifndef CEPETRA_EXPORT_CPP_HPP
#define CEPETRA_EXPORT_CPP_HPP


#include "CTrilinos_enums.h"
#include "Teuchos_RCP.hpp"
#include "Epetra_Export.h"


namespace CEpetra {


using Teuchos::RCP;


/*! get Epetra_Export from non-const table using CT_Epetra_Export_ID */
const RCP<Epetra_Export>
getExport( CT_Epetra_Export_ID_t id );

/*! get Epetra_Export from non-const table using CTrilinos_Universal_ID_t */
const RCP<Epetra_Export>
getExport( CTrilinos_Universal_ID_t id );

/*! get const Epetra_Export from either the const or non-const table
 * using CT_Epetra_Export_ID */
const RCP<const Epetra_Export>
getConstExport( CT_Epetra_Export_ID_t id );

/*! get const Epetra_Export from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const RCP<const Epetra_Export>
getConstExport( CTrilinos_Universal_ID_t id );

/*! store Epetra_Export (owned) in non-const table */
CT_Epetra_Export_ID_t
storeNewExport( Epetra_Export *pobj );

/*! store Epetra_Export in non-const table */
CT_Epetra_Export_ID_t
storeExport( Epetra_Export *pobj );

/*! store const Epetra_Export in const table */
CT_Epetra_Export_ID_t
storeConstExport( const Epetra_Export *pobj );

/* remove Epetra_Export from table using CT_Epetra_Export_ID */
void
removeExport( CT_Epetra_Export_ID_t *id );

/* purge Epetra_Export table */
void
purgeExport(  );

} // namespace CEpetra


#endif // CEPETRA_EXPORT_CPP_HPP


