
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
#include "CEpetra_Map.h"
#include "CEpetra_Map_Cpp.hpp"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CEpetra_Comm_Cpp.hpp"


//
// Definitions from CEpetra_Map.h
//


extern "C" {


CT_Epetra_Map_ID_t Epetra_Map_Degeneralize ( 
  CTrilinos_Universal_ID_t id )
{
    return CTrilinos::concreteType<CT_Epetra_Map_ID_t>(id);
}

CTrilinos_Universal_ID_t Epetra_Map_Generalize ( 
  CT_Epetra_Map_ID_t id )
{
    return CTrilinos::abstractType<CT_Epetra_Map_ID_t>(id);
}

CT_Epetra_Map_ID_t Epetra_Map_Create ( 
  int NumGlobalElements, int IndexBase, CT_Epetra_Comm_ID_t CommID )
{
    const Teuchos::RCP<const Epetra_Comm> Comm = CEpetra::getConstComm(CommID);
    return CEpetra::storeNewMap(new Epetra_Map(NumGlobalElements, IndexBase, 
        *Comm));
}

CT_Epetra_Map_ID_t Epetra_Map_Create_Linear ( 
  int NumGlobalElements, int NumMyElements, int IndexBase, 
  CT_Epetra_Comm_ID_t CommID )
{
    const Teuchos::RCP<const Epetra_Comm> Comm = CEpetra::getConstComm(CommID);
    return CEpetra::storeNewMap(new Epetra_Map(NumGlobalElements, 
        NumMyElements, IndexBase, *Comm));
}

CT_Epetra_Map_ID_t Epetra_Map_Create_Arbitrary ( 
  int NumGlobalElements, int NumMyElements, 
  const int * MyGlobalElements, int IndexBase, 
  CT_Epetra_Comm_ID_t CommID )
{
    const Teuchos::RCP<const Epetra_Comm> Comm = CEpetra::getConstComm(CommID);
    return CEpetra::storeNewMap(new Epetra_Map(NumGlobalElements, 
        NumMyElements, MyGlobalElements, IndexBase, *Comm));
}

CT_Epetra_Map_ID_t Epetra_Map_Duplicate ( CT_Epetra_Map_ID_t mapID )
{
    const Teuchos::RCP<const Epetra_Map> map = CEpetra::getConstMap(mapID);
    return CEpetra::storeNewMap(new Epetra_Map(*map));
}

void Epetra_Map_Destroy ( CT_Epetra_Map_ID_t * selfID )
{
    CEpetra::removeMap(selfID);
}

void Epetra_Map_Assign ( 
  CT_Epetra_Map_ID_t selfID, CT_Epetra_Map_ID_t mapID )
{
    Epetra_Map& self = *( CEpetra::getMap(selfID) );

    const Teuchos::RCP<const Epetra_Map> map = CEpetra::getConstMap(mapID);
    self = *map;
}


} // extern "C"




