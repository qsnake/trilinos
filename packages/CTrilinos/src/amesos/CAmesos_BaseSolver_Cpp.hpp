
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


#ifdef HAVE_CTRILINOS_AMESOS


#ifndef CAMESOS_BASESOLVER_CPP_HPP
#define CAMESOS_BASESOLVER_CPP_HPP


#include "CTrilinos_enums.h"
#include "Teuchos_RCP.hpp"
#include "Amesos_BaseSolver.h"


namespace CAmesos {


using Teuchos::RCP;


/*! get Amesos_BaseSolver from non-const table using CT_Amesos_BaseSolver_ID */
const RCP<Amesos_BaseSolver>
getBaseSolver( CT_Amesos_BaseSolver_ID_t id );

/*! get Amesos_BaseSolver from non-const table using CTrilinos_Universal_ID_t */
const RCP<Amesos_BaseSolver>
getBaseSolver( CTrilinos_Universal_ID_t id );

/*! get const Amesos_BaseSolver from either the const or non-const table
 * using CT_Amesos_BaseSolver_ID */
const RCP<const Amesos_BaseSolver>
getConstBaseSolver( CT_Amesos_BaseSolver_ID_t id );

/*! get const Amesos_BaseSolver from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const RCP<const Amesos_BaseSolver>
getConstBaseSolver( CTrilinos_Universal_ID_t id );

/*! store Amesos_BaseSolver (owned) in non-const table */
CT_Amesos_BaseSolver_ID_t
storeNewBaseSolver( Amesos_BaseSolver *pobj );

/*! store Amesos_BaseSolver in non-const table */
CT_Amesos_BaseSolver_ID_t
storeBaseSolver( Amesos_BaseSolver *pobj );

/*! store const Amesos_BaseSolver in const table */
CT_Amesos_BaseSolver_ID_t
storeConstBaseSolver( const Amesos_BaseSolver *pobj );

/* remove Amesos_BaseSolver from table using CT_Amesos_BaseSolver_ID */
void
removeBaseSolver( CT_Amesos_BaseSolver_ID_t *id );

/* purge Amesos_BaseSolver table */
void
purgeBaseSolver(  );

} // namespace CAmesos


#endif // CAMESOS_BASESOLVER_CPP_HPP


#endif /* HAVE_CTRILINOS_AMESOS */


