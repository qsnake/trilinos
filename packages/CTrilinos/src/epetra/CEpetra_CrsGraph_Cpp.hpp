
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

#ifndef CEPETRA_CRSGRAPH_CPP_HPP
#define CEPETRA_CRSGRAPH_CPP_HPP


#include "CTrilinos_enums.h"
#include "Teuchos_RCP.hpp"
#include "Epetra_CrsGraph.h"


namespace CEpetra {


using Teuchos::RCP;


/*! get Epetra_CrsGraph from non-const table using CT_Epetra_CrsGraph_ID */
const RCP<Epetra_CrsGraph>
getCrsGraph( CT_Epetra_CrsGraph_ID_t id );

/*! get Epetra_CrsGraph from non-const table using CTrilinos_Universal_ID_t */
const RCP<Epetra_CrsGraph>
getCrsGraph( CTrilinos_Universal_ID_t id );

/*! get const Epetra_CrsGraph from either the const or non-const table
 * using CT_Epetra_CrsGraph_ID */
const RCP<const Epetra_CrsGraph>
getConstCrsGraph( CT_Epetra_CrsGraph_ID_t id );

/*! get const Epetra_CrsGraph from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const RCP<const Epetra_CrsGraph>
getConstCrsGraph( CTrilinos_Universal_ID_t id );

/*! store Epetra_CrsGraph (owned) in non-const table */
CT_Epetra_CrsGraph_ID_t
storeNewCrsGraph( Epetra_CrsGraph *pobj );

/*! store Epetra_CrsGraph in non-const table */
CT_Epetra_CrsGraph_ID_t
storeCrsGraph( Epetra_CrsGraph *pobj );

/*! store const Epetra_CrsGraph in const table */
CT_Epetra_CrsGraph_ID_t
storeConstCrsGraph( const Epetra_CrsGraph *pobj );

/* remove Epetra_CrsGraph from table using CT_Epetra_CrsGraph_ID */
void
removeCrsGraph( CT_Epetra_CrsGraph_ID_t *id );

/* purge Epetra_CrsGraph table */
void
purgeCrsGraph(  );

} // namespace CEpetra


#endif // CEPETRA_CRSGRAPH_CPP_HPP


