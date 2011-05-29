
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
#include "CEpetra_DistObject.h"
#include "CEpetra_DistObject_Cpp.hpp"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CEpetra_SrcDistObject_Cpp.hpp"
#include "CEpetra_Import_Cpp.hpp"
#include "CEpetra_OffsetIndex_Cpp.hpp"
#include "CEpetra_Export_Cpp.hpp"
#include "CEpetra_BlockMap_Cpp.hpp"
#include "CEpetra_Comm_Cpp.hpp"


//
// Definitions from CEpetra_DistObject.h
//


extern "C" {


CT_Epetra_DistObject_ID_t Epetra_DistObject_Degeneralize ( 
  CTrilinos_Universal_ID_t id )
{
    return CTrilinos::concreteType<CT_Epetra_DistObject_ID_t>(id);
}

CTrilinos_Universal_ID_t Epetra_DistObject_Generalize ( 
  CT_Epetra_DistObject_ID_t id )
{
    return CTrilinos::abstractType<CT_Epetra_DistObject_ID_t>(id);
}

void Epetra_DistObject_Destroy ( CT_Epetra_DistObject_ID_t * selfID )
{
    CEpetra::removeDistObject(selfID);
}

int Epetra_DistObject_Import ( 
  CT_Epetra_DistObject_ID_t selfID, 
  CT_Epetra_SrcDistObject_ID_t AID, 
  CT_Epetra_Import_ID_t ImporterID, 
  CT_Epetra_CombineMode_E_t CombineMode, 
  CT_Epetra_OffsetIndex_ID_t IndexorID )
{
    const Teuchos::RCP<const Epetra_SrcDistObject> A = 
        CEpetra::getConstSrcDistObject(AID);
    const Teuchos::RCP<const Epetra_Import> Importer = CEpetra::getConstImport(
        ImporterID);
    const Teuchos::RCP<const Epetra_OffsetIndex> Indexor = 
        CEpetra::getConstOffsetIndex(IndexorID);
    return CEpetra::getDistObject(selfID)->Import(*A, *Importer, 
        (Epetra_CombineMode) CombineMode, Indexor.getRawPtr());
}

int Epetra_DistObject_Import_UsingExporter ( 
  CT_Epetra_DistObject_ID_t selfID, 
  CT_Epetra_SrcDistObject_ID_t AID, 
  CT_Epetra_Export_ID_t ExporterID, 
  CT_Epetra_CombineMode_E_t CombineMode, 
  CT_Epetra_OffsetIndex_ID_t IndexorID )
{
    const Teuchos::RCP<const Epetra_SrcDistObject> A = 
        CEpetra::getConstSrcDistObject(AID);
    const Teuchos::RCP<const Epetra_Export> Exporter = CEpetra::getConstExport(
        ExporterID);
    const Teuchos::RCP<const Epetra_OffsetIndex> Indexor = 
        CEpetra::getConstOffsetIndex(IndexorID);
    return CEpetra::getDistObject(selfID)->Import(*A, *Exporter, 
        (Epetra_CombineMode) CombineMode, Indexor.getRawPtr());
}

int Epetra_DistObject_Export_UsingImporter ( 
  CT_Epetra_DistObject_ID_t selfID, 
  CT_Epetra_SrcDistObject_ID_t AID, 
  CT_Epetra_Import_ID_t ImporterID, 
  CT_Epetra_CombineMode_E_t CombineMode, 
  CT_Epetra_OffsetIndex_ID_t IndexorID )
{
    const Teuchos::RCP<const Epetra_SrcDistObject> A = 
        CEpetra::getConstSrcDistObject(AID);
    const Teuchos::RCP<const Epetra_Import> Importer = CEpetra::getConstImport(
        ImporterID);
    const Teuchos::RCP<const Epetra_OffsetIndex> Indexor = 
        CEpetra::getConstOffsetIndex(IndexorID);
    return CEpetra::getDistObject(selfID)->Export(*A, *Importer, 
        (Epetra_CombineMode) CombineMode, Indexor.getRawPtr());
}

int Epetra_DistObject_Export ( 
  CT_Epetra_DistObject_ID_t selfID, 
  CT_Epetra_SrcDistObject_ID_t AID, 
  CT_Epetra_Export_ID_t ExporterID, 
  CT_Epetra_CombineMode_E_t CombineMode, 
  CT_Epetra_OffsetIndex_ID_t IndexorID )
{
    const Teuchos::RCP<const Epetra_SrcDistObject> A = 
        CEpetra::getConstSrcDistObject(AID);
    const Teuchos::RCP<const Epetra_Export> Exporter = CEpetra::getConstExport(
        ExporterID);
    const Teuchos::RCP<const Epetra_OffsetIndex> Indexor = 
        CEpetra::getConstOffsetIndex(IndexorID);
    return CEpetra::getDistObject(selfID)->Export(*A, *Exporter, 
        (Epetra_CombineMode) CombineMode, Indexor.getRawPtr());
}

CT_Epetra_BlockMap_ID_t Epetra_DistObject_Map ( 
  CT_Epetra_DistObject_ID_t selfID )
{
    return CEpetra::storeConstBlockMap(&( CEpetra::getConstDistObject(
        selfID)->Map() ));
}

CT_Epetra_Comm_ID_t Epetra_DistObject_Comm ( 
  CT_Epetra_DistObject_ID_t selfID )
{
    return CEpetra::storeConstComm(&( CEpetra::getConstDistObject(
        selfID)->Comm() ));
}

boolean Epetra_DistObject_DistributedGlobal ( 
  CT_Epetra_DistObject_ID_t selfID )
{
    return ((CEpetra::getConstDistObject(
        selfID)->DistributedGlobal()) ? TRUE : FALSE);
}


} // extern "C"




