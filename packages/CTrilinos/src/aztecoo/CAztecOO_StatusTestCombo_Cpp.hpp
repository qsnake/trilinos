
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


#ifndef CAZTECOO_STATUSTESTCOMBO_CPP_HPP
#define CAZTECOO_STATUSTESTCOMBO_CPP_HPP


#include "CTrilinos_enums.h"
#include "Teuchos_RCP.hpp"
#include "AztecOO_StatusTestCombo.h"


namespace CAztecOO {


using Teuchos::RCP;


/*! get AztecOO_StatusTestCombo from non-const table using CT_AztecOO_StatusTestCombo_ID */
const RCP<AztecOO_StatusTestCombo>
getStatusTestCombo( CT_AztecOO_StatusTestCombo_ID_t id );

/*! get AztecOO_StatusTestCombo from non-const table using CTrilinos_Universal_ID_t */
const RCP<AztecOO_StatusTestCombo>
getStatusTestCombo( CTrilinos_Universal_ID_t id );

/*! get const AztecOO_StatusTestCombo from either the const or non-const table
 * using CT_AztecOO_StatusTestCombo_ID */
const RCP<const AztecOO_StatusTestCombo>
getConstStatusTestCombo( CT_AztecOO_StatusTestCombo_ID_t id );

/*! get const AztecOO_StatusTestCombo from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const RCP<const AztecOO_StatusTestCombo>
getConstStatusTestCombo( CTrilinos_Universal_ID_t id );

/*! store AztecOO_StatusTestCombo (owned) in non-const table */
CT_AztecOO_StatusTestCombo_ID_t
storeNewStatusTestCombo( AztecOO_StatusTestCombo *pobj );

/*! store AztecOO_StatusTestCombo in non-const table */
CT_AztecOO_StatusTestCombo_ID_t
storeStatusTestCombo( AztecOO_StatusTestCombo *pobj );

/*! store const AztecOO_StatusTestCombo in const table */
CT_AztecOO_StatusTestCombo_ID_t
storeConstStatusTestCombo( const AztecOO_StatusTestCombo *pobj );

/* remove AztecOO_StatusTestCombo from table using CT_AztecOO_StatusTestCombo_ID */
void
removeStatusTestCombo( CT_AztecOO_StatusTestCombo_ID_t *id );

/* purge AztecOO_StatusTestCombo table */
void
purgeStatusTestCombo(  );

} // namespace CAztecOO


#endif // CAZTECOO_STATUSTESTCOMBO_CPP_HPP


#endif /* HAVE_CTRILINOS_AZTECOO */


