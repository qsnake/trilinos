
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

#ifndef CEPETRA_ROWMATRIX_CPP_HPP
#define CEPETRA_ROWMATRIX_CPP_HPP


#include "CTrilinos_enums.h"
#include "Teuchos_RCP.hpp"
#include "Epetra_RowMatrix.h"


namespace CEpetra {


using Teuchos::RCP;


/*! get Epetra_RowMatrix from non-const table using CT_Epetra_RowMatrix_ID */
const RCP<Epetra_RowMatrix>
getRowMatrix( CT_Epetra_RowMatrix_ID_t id );

/*! get Epetra_RowMatrix from non-const table using CTrilinos_Universal_ID_t */
const RCP<Epetra_RowMatrix>
getRowMatrix( CTrilinos_Universal_ID_t id );

/*! get const Epetra_RowMatrix from either the const or non-const table
 * using CT_Epetra_RowMatrix_ID */
const RCP<const Epetra_RowMatrix>
getConstRowMatrix( CT_Epetra_RowMatrix_ID_t id );

/*! get const Epetra_RowMatrix from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const RCP<const Epetra_RowMatrix>
getConstRowMatrix( CTrilinos_Universal_ID_t id );

/*! store Epetra_RowMatrix (owned) in non-const table */
CT_Epetra_RowMatrix_ID_t
storeNewRowMatrix( Epetra_RowMatrix *pobj );

/*! store Epetra_RowMatrix in non-const table */
CT_Epetra_RowMatrix_ID_t
storeRowMatrix( Epetra_RowMatrix *pobj );

/*! store const Epetra_RowMatrix in const table */
CT_Epetra_RowMatrix_ID_t
storeConstRowMatrix( const Epetra_RowMatrix *pobj );

/* remove Epetra_RowMatrix from table using CT_Epetra_RowMatrix_ID */
void
removeRowMatrix( CT_Epetra_RowMatrix_ID_t *id );

/* purge Epetra_RowMatrix table */
void
purgeRowMatrix(  );

} // namespace CEpetra


#endif // CEPETRA_ROWMATRIX_CPP_HPP


