
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


/*! @file CTrilinos_table_man.h
 * @brief Provides interface for CTrilinos table micro-management. */


#ifndef CTRILINOS_TABLE_MAN_H
#define CTRILINOS_TABLE_MAN_H


#include "CTrilinos_config.h"
#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif


/*! Copies the RCP from one table into a second table. The new ID
 *  will be returned from the function. Both the old and the new
 *  IDs will need to be removed from the tables in order to destroy
 *  the object. */
CTrilinos_Universal_ID_t CT_Alias(CTrilinos_Universal_ID_t aid, CTrilinos_Table_ID_t new_table);

/*! Removes the RCP from one table and puts it in another. *aid will
 *  hold the new struct value afterward. Only the new RCP will need
 *  to be removed in order to destroy the object. */
void CT_Migrate(CTrilinos_Universal_ID_t *aid, CTrilinos_Table_ID_t new_table);

/*! Checks to see if the underlying object referenced by a table
 *  entry is dynamic_cast'able to a given type (can be used to
 *  distinguish, e.g., an Epetra_SerialComm from an Epetra_MpiComm). */
boolean CT_TypeCheck(CTrilinos_Universal_ID_t aid, CTrilinos_Table_ID_t type);


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif
