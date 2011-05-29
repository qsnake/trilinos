
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
#include "CEpetra_OffsetIndex.h"
#include "CEpetra_OffsetIndex_Cpp.hpp"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CEpetra_CrsGraph_Cpp.hpp"
#include "CEpetra_Import_Cpp.hpp"
#include "CEpetra_Export_Cpp.hpp"


//
// Definitions from CEpetra_OffsetIndex.h
//


extern "C" {


CT_Epetra_OffsetIndex_ID_t Epetra_OffsetIndex_Degeneralize ( 
  CTrilinos_Universal_ID_t id )
{
    return CTrilinos::concreteType<CT_Epetra_OffsetIndex_ID_t>(id);
}

CTrilinos_Universal_ID_t Epetra_OffsetIndex_Generalize ( 
  CT_Epetra_OffsetIndex_ID_t id )
{
    return CTrilinos::abstractType<CT_Epetra_OffsetIndex_ID_t>(id);
}

CT_Epetra_OffsetIndex_ID_t Epetra_OffsetIndex_Create_FromImporter ( 
  CT_Epetra_CrsGraph_ID_t SourceGraphID, 
  CT_Epetra_CrsGraph_ID_t TargetGraphID, 
  CT_Epetra_Import_ID_t ImporterID )
{
    const Teuchos::RCP<const Epetra_CrsGraph> SourceGraph = 
        CEpetra::getConstCrsGraph(SourceGraphID);
    const Teuchos::RCP<const Epetra_CrsGraph> TargetGraph = 
        CEpetra::getConstCrsGraph(TargetGraphID);
    const Teuchos::RCP<Epetra_Import> Importer = CEpetra::getImport(
        ImporterID);
    return CEpetra::storeNewOffsetIndex(new Epetra_OffsetIndex(*SourceGraph, 
        *TargetGraph, *Importer));
}

CT_Epetra_OffsetIndex_ID_t Epetra_OffsetIndex_Create_FromExporter ( 
  CT_Epetra_CrsGraph_ID_t SourceGraphID, 
  CT_Epetra_CrsGraph_ID_t TargetGraphID, 
  CT_Epetra_Export_ID_t ExporterID )
{
    const Teuchos::RCP<const Epetra_CrsGraph> SourceGraph = 
        CEpetra::getConstCrsGraph(SourceGraphID);
    const Teuchos::RCP<const Epetra_CrsGraph> TargetGraph = 
        CEpetra::getConstCrsGraph(TargetGraphID);
    const Teuchos::RCP<Epetra_Export> Exporter = CEpetra::getExport(
        ExporterID);
    return CEpetra::storeNewOffsetIndex(new Epetra_OffsetIndex(*SourceGraph, 
        *TargetGraph, *Exporter));
}

CT_Epetra_OffsetIndex_ID_t Epetra_OffsetIndex_Duplicate ( 
  CT_Epetra_OffsetIndex_ID_t IndexorID )
{
    const Teuchos::RCP<const Epetra_OffsetIndex> Indexor = 
        CEpetra::getConstOffsetIndex(IndexorID);
    return CEpetra::storeNewOffsetIndex(new Epetra_OffsetIndex(*Indexor));
}

void Epetra_OffsetIndex_Destroy ( 
  CT_Epetra_OffsetIndex_ID_t * selfID )
{
    CEpetra::removeOffsetIndex(selfID);
}

int ** Epetra_OffsetIndex_SameOffsets ( 
  CT_Epetra_OffsetIndex_ID_t selfID )
{
    return CEpetra::getConstOffsetIndex(selfID)->SameOffsets();
}

int ** Epetra_OffsetIndex_PermuteOffsets ( 
  CT_Epetra_OffsetIndex_ID_t selfID )
{
    return CEpetra::getConstOffsetIndex(selfID)->PermuteOffsets();
}

int ** Epetra_OffsetIndex_RemoteOffsets ( 
  CT_Epetra_OffsetIndex_ID_t selfID )
{
    return CEpetra::getConstOffsetIndex(selfID)->RemoteOffsets();
}


} // extern "C"




