
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


#include "CTrilinos_enums.h"
#include "CAmesos_BaseSolver.h"
#include "CAmesos_BaseSolver_Cpp.hpp"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTeuchos_ParameterList_Cpp.hpp"
#include "CEpetra_LinearProblem_Cpp.hpp"
#include "CEpetra_Comm_Cpp.hpp"


//
// Definitions from CAmesos_BaseSolver.h
//


extern "C" {


CT_Amesos_BaseSolver_ID_t Amesos_BaseSolver_Degeneralize ( 
  CTrilinos_Universal_ID_t id )
{
    return CTrilinos::concreteType<CT_Amesos_BaseSolver_ID_t>(id);
}

CTrilinos_Universal_ID_t Amesos_BaseSolver_Generalize ( 
  CT_Amesos_BaseSolver_ID_t id )
{
    return CTrilinos::abstractType<CT_Amesos_BaseSolver_ID_t>(id);
}

void Amesos_BaseSolver_Destroy ( CT_Amesos_BaseSolver_ID_t * selfID )
{
    CAmesos::removeBaseSolver(selfID);
}

int Amesos_BaseSolver_SymbolicFactorization ( 
  CT_Amesos_BaseSolver_ID_t selfID )
{
    return CAmesos::getBaseSolver(selfID)->SymbolicFactorization();
}

int Amesos_BaseSolver_NumericFactorization ( 
  CT_Amesos_BaseSolver_ID_t selfID )
{
    return CAmesos::getBaseSolver(selfID)->NumericFactorization();
}

int Amesos_BaseSolver_Solve ( CT_Amesos_BaseSolver_ID_t selfID )
{
    return CAmesos::getBaseSolver(selfID)->Solve();
}

int Amesos_BaseSolver_SetUseTranspose ( 
  CT_Amesos_BaseSolver_ID_t selfID, boolean UseTranspose )
{
    return CAmesos::getBaseSolver(selfID)->SetUseTranspose(((UseTranspose) != 
        FALSE ? true : false));
}

boolean Amesos_BaseSolver_UseTranspose ( 
  CT_Amesos_BaseSolver_ID_t selfID )
{
    return ((CAmesos::getConstBaseSolver(
        selfID)->UseTranspose()) ? TRUE : FALSE);
}

int Amesos_BaseSolver_SetParameters ( 
  CT_Amesos_BaseSolver_ID_t selfID, 
  CT_Teuchos_ParameterList_ID_t ParameterListID )
{
    const Teuchos::RCP<Teuchos::ParameterList> ParameterList = 
        CTeuchos::getParameterList(ParameterListID);
    return CAmesos::getBaseSolver(selfID)->SetParameters(*ParameterList);
}

CT_Epetra_LinearProblem_ID_t Amesos_BaseSolver_GetProblem ( 
  CT_Amesos_BaseSolver_ID_t selfID )
{
    return CEpetra::storeConstLinearProblem(CAmesos::getConstBaseSolver(
        selfID)->GetProblem());
}

boolean Amesos_BaseSolver_MatrixShapeOK ( 
  CT_Amesos_BaseSolver_ID_t selfID )
{
    return ((CAmesos::getConstBaseSolver(
        selfID)->MatrixShapeOK()) ? TRUE : FALSE);
}

CT_Epetra_Comm_ID_t Amesos_BaseSolver_Comm ( 
  CT_Amesos_BaseSolver_ID_t selfID )
{
    return CEpetra::storeConstComm(&( CAmesos::getConstBaseSolver(
        selfID)->Comm() ));
}

int Amesos_BaseSolver_NumSymbolicFact ( 
  CT_Amesos_BaseSolver_ID_t selfID )
{
    return CAmesos::getConstBaseSolver(selfID)->NumSymbolicFact();
}

int Amesos_BaseSolver_NumNumericFact ( 
  CT_Amesos_BaseSolver_ID_t selfID )
{
    return CAmesos::getConstBaseSolver(selfID)->NumNumericFact();
}

int Amesos_BaseSolver_NumSolve ( CT_Amesos_BaseSolver_ID_t selfID )
{
    return CAmesos::getConstBaseSolver(selfID)->NumSolve();
}

void Amesos_BaseSolver_PrintStatus ( 
  CT_Amesos_BaseSolver_ID_t selfID )
{
    CAmesos::getConstBaseSolver(selfID)->PrintStatus();
}

void Amesos_BaseSolver_PrintTiming ( 
  CT_Amesos_BaseSolver_ID_t selfID )
{
    CAmesos::getConstBaseSolver(selfID)->PrintTiming();
}

void Amesos_BaseSolver_setParameterList ( 
  CT_Amesos_BaseSolver_ID_t selfID, 
  CT_Teuchos_ParameterList_ID_t paramListID )
{
    const Teuchos::RCP<Teuchos::ParameterList> paramList = 
        CTeuchos::getParameterList(paramListID);
    CAmesos::getBaseSolver(selfID)->setParameterList(paramList);
}

CT_Teuchos_ParameterList_ID_t Amesos_BaseSolver_getNonconstParameterList ( 
  CT_Amesos_BaseSolver_ID_t selfID )
{
    return CTeuchos::storeParameterList(CAmesos::getBaseSolver(
        selfID)->getNonconstParameterList().getRawPtr());
}

CT_Teuchos_ParameterList_ID_t Amesos_BaseSolver_unsetParameterList ( 
  CT_Amesos_BaseSolver_ID_t selfID )
{
    return CTeuchos::storeParameterList(CAmesos::getBaseSolver(
        selfID)->unsetParameterList().getRawPtr());
}

void Amesos_BaseSolver_GetTiming ( 
  CT_Amesos_BaseSolver_ID_t selfID, 
  CT_Teuchos_ParameterList_ID_t TimingParameterListID )
{
    const Teuchos::RCP<Teuchos::ParameterList> TimingParameterList = 
        CTeuchos::getParameterList(TimingParameterListID);
    CAmesos::getConstBaseSolver(selfID)->GetTiming(*TimingParameterList);
}


} // extern "C"




#endif /* HAVE_CTRILINOS_AMESOS */


