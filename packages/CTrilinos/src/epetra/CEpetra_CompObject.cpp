
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

#include "CTrilinos_enums.h"
#include "CEpetra_CompObject.h"
#include "CEpetra_CompObject_Cpp.hpp"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CEpetra_Flops_Cpp.hpp"


//
// Definitions from CEpetra_CompObject.h
//


extern "C" {


CT_Epetra_CompObject_ID_t Epetra_CompObject_Degeneralize ( 
  CTrilinos_Universal_ID_t id )
{
    return CTrilinos::concreteType<CT_Epetra_CompObject_ID_t>(id);
}

CTrilinos_Universal_ID_t Epetra_CompObject_Generalize ( 
  CT_Epetra_CompObject_ID_t id )
{
    return CTrilinos::abstractType<CT_Epetra_CompObject_ID_t>(id);
}

CT_Epetra_CompObject_ID_t Epetra_CompObject_Create (  )
{
    return CEpetra::storeNewCompObject(new Epetra_CompObject());
}

CT_Epetra_CompObject_ID_t Epetra_CompObject_Duplicate ( 
  CT_Epetra_CompObject_ID_t SourceID )
{
    const Teuchos::RCP<const Epetra_CompObject> Source = 
        CEpetra::getConstCompObject(SourceID);
    return CEpetra::storeNewCompObject(new Epetra_CompObject(*Source));
}

void Epetra_CompObject_Destroy ( CT_Epetra_CompObject_ID_t * selfID )
{
    CEpetra::removeCompObject(selfID);
}

void Epetra_CompObject_SetFlopCounter ( 
  CT_Epetra_CompObject_ID_t selfID, 
  CT_Epetra_Flops_ID_t FlopCounter_inID )
{
    const Teuchos::RCP<const Epetra_Flops> FlopCounter_in = 
        CEpetra::getConstFlops(FlopCounter_inID);
    CEpetra::getCompObject(selfID)->SetFlopCounter(*FlopCounter_in);
}

void Epetra_CompObject_SetFlopCounter_Matching ( 
  CT_Epetra_CompObject_ID_t selfID, 
  CT_Epetra_CompObject_ID_t CompObjectID )
{
    const Teuchos::RCP<const Epetra_CompObject> CompObject = 
        CEpetra::getConstCompObject(CompObjectID);
    CEpetra::getCompObject(selfID)->SetFlopCounter(*CompObject);
}

void Epetra_CompObject_UnsetFlopCounter ( 
  CT_Epetra_CompObject_ID_t selfID )
{
    CEpetra::getCompObject(selfID)->UnsetFlopCounter();
}

CT_Epetra_Flops_ID_t Epetra_CompObject_GetFlopCounter ( 
  CT_Epetra_CompObject_ID_t selfID )
{
    return CEpetra::storeFlops(CEpetra::getConstCompObject(
        selfID)->GetFlopCounter());
}

void Epetra_CompObject_ResetFlops ( 
  CT_Epetra_CompObject_ID_t selfID )
{
    CEpetra::getConstCompObject(selfID)->ResetFlops();
}

double Epetra_CompObject_Flops ( CT_Epetra_CompObject_ID_t selfID )
{
    return CEpetra::getConstCompObject(selfID)->Flops();
}

void Epetra_CompObject_UpdateFlops_Int ( 
  CT_Epetra_CompObject_ID_t selfID, int Flops_in )
{
    CEpetra::getConstCompObject(selfID)->UpdateFlops(Flops_in);
}

void Epetra_CompObject_UpdateFlops_Long ( 
  CT_Epetra_CompObject_ID_t selfID, long int Flops_in )
{
    CEpetra::getConstCompObject(selfID)->UpdateFlops(Flops_in);
}

void Epetra_CompObject_UpdateFlops_Double ( 
  CT_Epetra_CompObject_ID_t selfID, double Flops_in )
{
    CEpetra::getConstCompObject(selfID)->UpdateFlops(Flops_in);
}

void Epetra_CompObject_UpdateFlops_Float ( 
  CT_Epetra_CompObject_ID_t selfID, float Flops_in )
{
    CEpetra::getConstCompObject(selfID)->UpdateFlops(Flops_in);
}

void Epetra_CompObject_Assign ( 
  CT_Epetra_CompObject_ID_t selfID, CT_Epetra_CompObject_ID_t srcID )
{
    Epetra_CompObject& self = *( CEpetra::getCompObject(selfID) );

    const Teuchos::RCP<const Epetra_CompObject> src = 
        CEpetra::getConstCompObject(srcID);
    self = *src;
}


} // extern "C"




