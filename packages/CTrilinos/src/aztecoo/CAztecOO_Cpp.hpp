
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


#ifdef HAVE_CTRILINOS_AZTECOO


#ifndef CAZTECOO_CPP_HPP
#define CAZTECOO_CPP_HPP


#include "CTrilinos_enums.h"
#include "Teuchos_RCP.hpp"
#include "AztecOO.h"


namespace CAztecOO {


using Teuchos::RCP;


/*! get AztecOO from non-const table using CT_AztecOO_ID */
const RCP<AztecOO>
getAztecOO( CT_AztecOO_ID_t id );

/*! get AztecOO from non-const table using CTrilinos_Universal_ID_t */
const RCP<AztecOO>
getAztecOO( CTrilinos_Universal_ID_t id );

/*! get const AztecOO from either the const or non-const table
 * using CT_AztecOO_ID */
const RCP<const AztecOO>
getConstAztecOO( CT_AztecOO_ID_t id );

/*! get const AztecOO from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const RCP<const AztecOO>
getConstAztecOO( CTrilinos_Universal_ID_t id );

/*! store AztecOO (owned) in non-const table */
CT_AztecOO_ID_t
storeNewAztecOO( AztecOO *pobj );

/*! store AztecOO in non-const table */
CT_AztecOO_ID_t
storeAztecOO( AztecOO *pobj );

/*! store const AztecOO in const table */
CT_AztecOO_ID_t
storeConstAztecOO( const AztecOO *pobj );

/* remove AztecOO from table using CT_AztecOO_ID */
void
removeAztecOO( CT_AztecOO_ID_t *id );

/* purge AztecOO table */
void
purgeAztecOO(  );

} // namespace CAztecOO


#endif // CAZTECOO_CPP_HPP


#endif /* HAVE_CTRILINOS_AZTECOO */


