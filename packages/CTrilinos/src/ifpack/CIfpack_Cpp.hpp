
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


#ifdef HAVE_CTRILINOS_IFPACK


#ifndef CIFPACK_CPP_HPP
#define CIFPACK_CPP_HPP


#include "CTrilinos_enums.h"
#include "Teuchos_RCP.hpp"
#include "Ifpack.h"


namespace CIfpack {


using Teuchos::RCP;


/*! get Ifpack from non-const table using CT_Ifpack_ID */
const RCP<Ifpack>
getIfpack( CT_Ifpack_ID_t id );

/*! get Ifpack from non-const table using CTrilinos_Universal_ID_t */
const RCP<Ifpack>
getIfpack( CTrilinos_Universal_ID_t id );

/*! get const Ifpack from either the const or non-const table
 * using CT_Ifpack_ID */
const RCP<const Ifpack>
getConstIfpack( CT_Ifpack_ID_t id );

/*! get const Ifpack from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const RCP<const Ifpack>
getConstIfpack( CTrilinos_Universal_ID_t id );

/*! store Ifpack (owned) in non-const table */
CT_Ifpack_ID_t
storeNewIfpack( Ifpack *pobj );

/*! store Ifpack in non-const table */
CT_Ifpack_ID_t
storeIfpack( Ifpack *pobj );

/*! store const Ifpack in const table */
CT_Ifpack_ID_t
storeConstIfpack( const Ifpack *pobj );

/* remove Ifpack from table using CT_Ifpack_ID */
void
removeIfpack( CT_Ifpack_ID_t *id );

/* purge Ifpack table */
void
purgeIfpack(  );

} // namespace CIfpack


#endif // CIFPACK_CPP_HPP


#endif /* HAVE_CTRILINOS_IFPACK */


