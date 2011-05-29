
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



/*! @file CAmesos.h
 * @brief Wrappers for Amesos */

/* True C header file! */


#ifndef CAMESOS_H
#define CAMESOS_H


#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif



/*! @name Amesos constructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   Amesos::Amesos()
*/
CT_Amesos_ID_t Amesos_Create (  );

/*@}*/

/*! @name Amesos destructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   Amesos::~Amesos()
*/
void Amesos_Destroy ( CT_Amesos_ID_t * selfID );

/*@}*/

/*! @name Amesos member wrappers */
/*@{*/

/*! @brief Wrapper for 
   Amesos_BaseSolver *Amesos::Create(const char *ClassType, const Epetra_LinearProblem& LinearProblem )
*/
CT_Amesos_BaseSolver_ID_t Amesos_CreateSolver ( 
  CT_Amesos_ID_t selfID, const char * ClassType, 
  CT_Epetra_LinearProblem_ID_t LinearProblemID );

/*! @brief Wrapper for 
   bool Amesos::Query(const char * ClassType)
*/
boolean Amesos_Query ( 
  CT_Amesos_ID_t selfID, const char * ClassType );

/*@}*/

/*! @name Amesos static function wrappers */
/*@{*/

/*! @brief Wrapper for 
   static Teuchos::ParameterList Amesos::GetValidParameters()
*/
CT_Teuchos_ParameterList_ID_t Amesos_GetValidParameters (  );

/*@}*/


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CAMESOS_H */

#endif /* HAVE_CTRILINOS_AMESOS */


