
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
#include "CEpetra_Import.h"
#include "CEpetra_Import_Cpp.hpp"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CEpetra_BlockMap_Cpp.hpp"
#include "CEpetra_Distributor_Cpp.hpp"


//
// Definitions from CEpetra_Import.h
//


extern "C" {


CT_Epetra_Import_ID_t Epetra_Import_Degeneralize ( 
  CTrilinos_Universal_ID_t id )
{
    return CTrilinos::concreteType<CT_Epetra_Import_ID_t>(id);
}

CTrilinos_Universal_ID_t Epetra_Import_Generalize ( 
  CT_Epetra_Import_ID_t id )
{
    return CTrilinos::abstractType<CT_Epetra_Import_ID_t>(id);
}

CT_Epetra_Import_ID_t Epetra_Import_Create ( 
  CT_Epetra_BlockMap_ID_t TargetMapID, 
  CT_Epetra_BlockMap_ID_t SourceMapID )
{
    const Teuchos::RCP<const Epetra_BlockMap> TargetMap = 
        CEpetra::getConstBlockMap(TargetMapID);
    const Teuchos::RCP<const Epetra_BlockMap> SourceMap = 
        CEpetra::getConstBlockMap(SourceMapID);
    return CEpetra::storeNewImport(new Epetra_Import(*TargetMap, *SourceMap));
}

CT_Epetra_Import_ID_t Epetra_Import_Duplicate ( 
  CT_Epetra_Import_ID_t ImporterID )
{
    const Teuchos::RCP<const Epetra_Import> Importer = CEpetra::getConstImport(
        ImporterID);
    return CEpetra::storeNewImport(new Epetra_Import(*Importer));
}

void Epetra_Import_Destroy ( CT_Epetra_Import_ID_t * selfID )
{
    CEpetra::removeImport(selfID);
}

int Epetra_Import_NumSameIDs ( CT_Epetra_Import_ID_t selfID )
{
    return CEpetra::getConstImport(selfID)->NumSameIDs();
}

int Epetra_Import_NumPermuteIDs ( CT_Epetra_Import_ID_t selfID )
{
    return CEpetra::getConstImport(selfID)->NumPermuteIDs();
}

int * Epetra_Import_PermuteFromLIDs ( CT_Epetra_Import_ID_t selfID )
{
    return CEpetra::getConstImport(selfID)->PermuteFromLIDs();
}

int * Epetra_Import_PermuteToLIDs ( CT_Epetra_Import_ID_t selfID )
{
    return CEpetra::getConstImport(selfID)->PermuteToLIDs();
}

int Epetra_Import_NumRemoteIDs ( CT_Epetra_Import_ID_t selfID )
{
    return CEpetra::getConstImport(selfID)->NumRemoteIDs();
}

int * Epetra_Import_RemoteLIDs ( CT_Epetra_Import_ID_t selfID )
{
    return CEpetra::getConstImport(selfID)->RemoteLIDs();
}

int Epetra_Import_NumExportIDs ( CT_Epetra_Import_ID_t selfID )
{
    return CEpetra::getConstImport(selfID)->NumExportIDs();
}

int * Epetra_Import_ExportLIDs ( CT_Epetra_Import_ID_t selfID )
{
    return CEpetra::getConstImport(selfID)->ExportLIDs();
}

int * Epetra_Import_ExportPIDs ( CT_Epetra_Import_ID_t selfID )
{
    return CEpetra::getConstImport(selfID)->ExportPIDs();
}

int Epetra_Import_NumSend ( CT_Epetra_Import_ID_t selfID )
{
    return CEpetra::getConstImport(selfID)->NumSend();
}

int Epetra_Import_NumRecv ( CT_Epetra_Import_ID_t selfID )
{
    return CEpetra::getConstImport(selfID)->NumRecv();
}

CT_Epetra_BlockMap_ID_t Epetra_Import_SourceMap ( 
  CT_Epetra_Import_ID_t selfID )
{
    return CEpetra::storeConstBlockMap(&( CEpetra::getConstImport(
        selfID)->SourceMap() ));
}

CT_Epetra_BlockMap_ID_t Epetra_Import_TargetMap ( 
  CT_Epetra_Import_ID_t selfID )
{
    return CEpetra::storeConstBlockMap(&( CEpetra::getConstImport(
        selfID)->TargetMap() ));
}

CT_Epetra_Distributor_ID_t Epetra_Import_Distributor ( 
  CT_Epetra_Import_ID_t selfID )
{
    return CEpetra::storeDistributor(&( CEpetra::getConstImport(
        selfID)->Distributor() ));
}


} // extern "C"




