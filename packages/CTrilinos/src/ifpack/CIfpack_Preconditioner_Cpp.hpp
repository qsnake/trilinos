
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


#ifndef CIFPACK_PRECONDITIONER_CPP_HPP
#define CIFPACK_PRECONDITIONER_CPP_HPP


#include "CTrilinos_enums.h"
#include "Teuchos_RCP.hpp"
#include "Ifpack_Preconditioner.h"


namespace CIfpack {


using Teuchos::RCP;


/*! get Ifpack_Preconditioner from non-const table using CT_Ifpack_Preconditioner_ID */
const RCP<Ifpack_Preconditioner>
getPreconditioner( CT_Ifpack_Preconditioner_ID_t id );

/*! get Ifpack_Preconditioner from non-const table using CTrilinos_Universal_ID_t */
const RCP<Ifpack_Preconditioner>
getPreconditioner( CTrilinos_Universal_ID_t id );

/*! get const Ifpack_Preconditioner from either the const or non-const table
 * using CT_Ifpack_Preconditioner_ID */
const RCP<const Ifpack_Preconditioner>
getConstPreconditioner( CT_Ifpack_Preconditioner_ID_t id );

/*! get const Ifpack_Preconditioner from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const RCP<const Ifpack_Preconditioner>
getConstPreconditioner( CTrilinos_Universal_ID_t id );

/*! store Ifpack_Preconditioner (owned) in non-const table */
CT_Ifpack_Preconditioner_ID_t
storeNewPreconditioner( Ifpack_Preconditioner *pobj );

/*! store Ifpack_Preconditioner in non-const table */
CT_Ifpack_Preconditioner_ID_t
storePreconditioner( Ifpack_Preconditioner *pobj );

/*! store const Ifpack_Preconditioner in const table */
CT_Ifpack_Preconditioner_ID_t
storeConstPreconditioner( const Ifpack_Preconditioner *pobj );

/* remove Ifpack_Preconditioner from table using CT_Ifpack_Preconditioner_ID */
void
removePreconditioner( CT_Ifpack_Preconditioner_ID_t *id );

/* purge Ifpack_Preconditioner table */
void
purgePreconditioner(  );

} // namespace CIfpack


#endif // CIFPACK_PRECONDITIONER_CPP_HPP


#endif /* HAVE_CTRILINOS_IFPACK */


