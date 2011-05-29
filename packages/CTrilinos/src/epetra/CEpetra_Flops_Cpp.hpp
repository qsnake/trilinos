
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

#ifndef CEPETRA_FLOPS_CPP_HPP
#define CEPETRA_FLOPS_CPP_HPP


#include "CTrilinos_enums.h"
#include "Teuchos_RCP.hpp"
#include "Epetra_Flops.h"


namespace CEpetra {


using Teuchos::RCP;


/*! get Epetra_Flops from non-const table using CT_Epetra_Flops_ID */
const RCP<Epetra_Flops>
getFlops( CT_Epetra_Flops_ID_t id );

/*! get Epetra_Flops from non-const table using CTrilinos_Universal_ID_t */
const RCP<Epetra_Flops>
getFlops( CTrilinos_Universal_ID_t id );

/*! get const Epetra_Flops from either the const or non-const table
 * using CT_Epetra_Flops_ID */
const RCP<const Epetra_Flops>
getConstFlops( CT_Epetra_Flops_ID_t id );

/*! get const Epetra_Flops from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const RCP<const Epetra_Flops>
getConstFlops( CTrilinos_Universal_ID_t id );

/*! store Epetra_Flops (owned) in non-const table */
CT_Epetra_Flops_ID_t
storeNewFlops( Epetra_Flops *pobj );

/*! store Epetra_Flops in non-const table */
CT_Epetra_Flops_ID_t
storeFlops( Epetra_Flops *pobj );

/*! store const Epetra_Flops in const table */
CT_Epetra_Flops_ID_t
storeConstFlops( const Epetra_Flops *pobj );

/* remove Epetra_Flops from table using CT_Epetra_Flops_ID */
void
removeFlops( CT_Epetra_Flops_ID_t *id );

/* purge Epetra_Flops table */
void
purgeFlops(  );

} // namespace CEpetra


#endif // CEPETRA_FLOPS_CPP_HPP


