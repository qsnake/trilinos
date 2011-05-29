
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


/*! @file CTrilinos_external_utils.h
 * @brief Provides interface for CTrilinos utilities needed mostly from Fortran. */


#ifndef CTRILINOS_EXTERNAL_UTILS_H
#define CTRILINOS_EXTERNAL_UTILS_H


#include "CTrilinos_config.h"
#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif


#ifdef HAVE_MPI

/*! Create an Epetra_MpiComm from Fortran */
CT_Epetra_MpiComm_ID_t Epetra_MpiComm_Fortran_Create ( int fcomm );

#endif /* HAVE_MPI */


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif
