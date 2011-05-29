
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


#include "CTrilinos_enums.h"
#include "CIfpack_Preconditioner.h"
#include "CIfpack_Preconditioner_Cpp.hpp"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTeuchos_ParameterList_Cpp.hpp"
#include "CEpetra_MultiVector_Cpp.hpp"
#include "CEpetra_RowMatrix_Cpp.hpp"


//
// Definitions from CIfpack_Preconditioner.h
//


extern "C" {


CT_Ifpack_Preconditioner_ID_t Ifpack_Preconditioner_Degeneralize ( 
  CTrilinos_Universal_ID_t id )
{
    return CTrilinos::concreteType<CT_Ifpack_Preconditioner_ID_t>(id);
}

CTrilinos_Universal_ID_t Ifpack_Preconditioner_Generalize ( 
  CT_Ifpack_Preconditioner_ID_t id )
{
    return CTrilinos::abstractType<CT_Ifpack_Preconditioner_ID_t>(id);
}

int Ifpack_Preconditioner_SetParameters ( 
  CT_Ifpack_Preconditioner_ID_t selfID, 
  CT_Teuchos_ParameterList_ID_t ListID )
{
    const Teuchos::RCP<Teuchos::ParameterList> List = 
        CTeuchos::getParameterList(ListID);
    return CIfpack::getPreconditioner(selfID)->SetParameters(*List);
}

int Ifpack_Preconditioner_Initialize ( 
  CT_Ifpack_Preconditioner_ID_t selfID )
{
    return CIfpack::getPreconditioner(selfID)->Initialize();
}

boolean Ifpack_Preconditioner_IsInitialized ( 
  CT_Ifpack_Preconditioner_ID_t selfID )
{
    return ((CIfpack::getConstPreconditioner(
        selfID)->IsInitialized()) ? TRUE : FALSE);
}

int Ifpack_Preconditioner_Compute ( 
  CT_Ifpack_Preconditioner_ID_t selfID )
{
    return CIfpack::getPreconditioner(selfID)->Compute();
}

boolean Ifpack_Preconditioner_IsComputed ( 
  CT_Ifpack_Preconditioner_ID_t selfID )
{
    return ((CIfpack::getConstPreconditioner(
        selfID)->IsComputed()) ? TRUE : FALSE);
}

double Ifpack_Preconditioner_Condest ( 
  CT_Ifpack_Preconditioner_ID_t selfID )
{
    return CIfpack::getConstPreconditioner(selfID)->Condest();
}

int Ifpack_Preconditioner_ApplyInverse ( 
  CT_Ifpack_Preconditioner_ID_t selfID, 
  CT_Epetra_MultiVector_ID_t XID, CT_Epetra_MultiVector_ID_t YID )
{
    const Teuchos::RCP<const Epetra_MultiVector> X = 
        CEpetra::getConstMultiVector(XID);
    const Teuchos::RCP<Epetra_MultiVector> Y = CEpetra::getMultiVector(YID);
    return CIfpack::getConstPreconditioner(selfID)->ApplyInverse(*X, *Y);
}

CT_Epetra_RowMatrix_ID_t Ifpack_Preconditioner_Matrix ( 
  CT_Ifpack_Preconditioner_ID_t selfID )
{
    return CEpetra::storeConstRowMatrix(&( CIfpack::getConstPreconditioner(
        selfID)->Matrix() ));
}

int Ifpack_Preconditioner_NumInitialize ( 
  CT_Ifpack_Preconditioner_ID_t selfID )
{
    return CIfpack::getConstPreconditioner(selfID)->NumInitialize();
}

int Ifpack_Preconditioner_NumCompute ( 
  CT_Ifpack_Preconditioner_ID_t selfID )
{
    return CIfpack::getConstPreconditioner(selfID)->NumCompute();
}

int Ifpack_Preconditioner_NumApplyInverse ( 
  CT_Ifpack_Preconditioner_ID_t selfID )
{
    return CIfpack::getConstPreconditioner(selfID)->NumApplyInverse();
}

double Ifpack_Preconditioner_InitializeTime ( 
  CT_Ifpack_Preconditioner_ID_t selfID )
{
    return CIfpack::getConstPreconditioner(selfID)->InitializeTime();
}

double Ifpack_Preconditioner_ComputeTime ( 
  CT_Ifpack_Preconditioner_ID_t selfID )
{
    return CIfpack::getConstPreconditioner(selfID)->ComputeTime();
}

double Ifpack_Preconditioner_ApplyInverseTime ( 
  CT_Ifpack_Preconditioner_ID_t selfID )
{
    return CIfpack::getConstPreconditioner(selfID)->ApplyInverseTime();
}

double Ifpack_Preconditioner_InitializeFlops ( 
  CT_Ifpack_Preconditioner_ID_t selfID )
{
    return CIfpack::getConstPreconditioner(selfID)->InitializeFlops();
}

double Ifpack_Preconditioner_ComputeFlops ( 
  CT_Ifpack_Preconditioner_ID_t selfID )
{
    return CIfpack::getConstPreconditioner(selfID)->ComputeFlops();
}

double Ifpack_Preconditioner_ApplyInverseFlops ( 
  CT_Ifpack_Preconditioner_ID_t selfID )
{
    return CIfpack::getConstPreconditioner(selfID)->ApplyInverseFlops();
}


} // extern "C"




#endif /* HAVE_CTRILINOS_IFPACK */


