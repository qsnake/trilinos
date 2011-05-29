
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
#include "CEpetra_Flops.h"
#include "CEpetra_Flops_Cpp.hpp"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"


//
// Definitions from CEpetra_Flops.h
//


extern "C" {


CT_Epetra_Flops_ID_t Epetra_Flops_Degeneralize ( 
  CTrilinos_Universal_ID_t id )
{
    return CTrilinos::concreteType<CT_Epetra_Flops_ID_t>(id);
}

CTrilinos_Universal_ID_t Epetra_Flops_Generalize ( 
  CT_Epetra_Flops_ID_t id )
{
    return CTrilinos::abstractType<CT_Epetra_Flops_ID_t>(id);
}

CT_Epetra_Flops_ID_t Epetra_Flops_Create (  )
{
    return CEpetra::storeNewFlops(new Epetra_Flops());
}

CT_Epetra_Flops_ID_t Epetra_Flops_Duplicate ( 
  CT_Epetra_Flops_ID_t Flops_inID )
{
    const Teuchos::RCP<const Epetra_Flops> Flops_in = CEpetra::getConstFlops(
        Flops_inID);
    return CEpetra::storeNewFlops(new Epetra_Flops(*Flops_in));
}

void Epetra_Flops_Destroy ( CT_Epetra_Flops_ID_t * selfID )
{
    CEpetra::removeFlops(selfID);
}

double Epetra_Flops_Flops ( CT_Epetra_Flops_ID_t selfID )
{
    return CEpetra::getConstFlops(selfID)->Flops();
}

void Epetra_Flops_ResetFlops ( CT_Epetra_Flops_ID_t selfID )
{
    CEpetra::getFlops(selfID)->ResetFlops();
}

void Epetra_Flops_Assign ( 
  CT_Epetra_Flops_ID_t selfID, CT_Epetra_Flops_ID_t srcID )
{
    Epetra_Flops& self = *( CEpetra::getFlops(selfID) );

    const Teuchos::RCP<const Epetra_Flops> src = CEpetra::getConstFlops(srcID);
    self = *src;
}


} // extern "C"




