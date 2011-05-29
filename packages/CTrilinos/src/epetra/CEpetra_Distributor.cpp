
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
#include "CEpetra_Distributor.h"
#include "CEpetra_Distributor_Cpp.hpp"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"


//
// Definitions from CEpetra_Distributor.h
//


extern "C" {


CT_Epetra_Distributor_ID_t Epetra_Distributor_Degeneralize ( 
  CTrilinos_Universal_ID_t id )
{
    return CTrilinos::concreteType<CT_Epetra_Distributor_ID_t>(id);
}

CTrilinos_Universal_ID_t Epetra_Distributor_Generalize ( 
  CT_Epetra_Distributor_ID_t id )
{
    return CTrilinos::abstractType<CT_Epetra_Distributor_ID_t>(id);
}

void Epetra_Distributor_Destroy ( 
  CT_Epetra_Distributor_ID_t * selfID )
{
    CEpetra::removeDistributor(selfID);
}

CT_Epetra_Distributor_ID_t Epetra_Distributor_Clone ( 
  CT_Epetra_Distributor_ID_t selfID )
{
    return CEpetra::storeDistributor(CEpetra::getDistributor(selfID)->Clone());
}

int Epetra_Distributor_CreateFromSends ( 
  CT_Epetra_Distributor_ID_t selfID, int NumExportIDs, 
  const int * ExportPIDs, boolean Deterministic, 
  int * NumRemoteIDs )
{
    return CEpetra::getDistributor(selfID)->CreateFromSends(NumExportIDs, 
        ExportPIDs, ((Deterministic) != FALSE ? true : false), *NumRemoteIDs);
}

int Epetra_Distributor_CreateFromRecvs ( 
  CT_Epetra_Distributor_ID_t selfID, int NumRemoteIDs, 
  const int * RemoteGIDs, const int * RemotePIDs, 
  boolean Deterministic, int * NumExportIDs, int ** ExportGIDs, 
  int ** ExportPIDs )
{
    return CEpetra::getDistributor(selfID)->CreateFromRecvs(NumRemoteIDs, 
        RemoteGIDs, RemotePIDs, ((Deterministic) != FALSE ? true : false), 
        *NumExportIDs, *ExportGIDs, *ExportPIDs);
}

int Epetra_Distributor_Do ( 
  CT_Epetra_Distributor_ID_t selfID, char * export_objs, 
  int obj_size, int * len_import_objs, char ** import_objs )
{
    return CEpetra::getDistributor(selfID)->Do(export_objs, obj_size, 
        *len_import_objs, *import_objs);
}

int Epetra_Distributor_DoReverse ( 
  CT_Epetra_Distributor_ID_t selfID, char * export_objs, 
  int obj_size, int * len_import_objs, char ** import_objs )
{
    return CEpetra::getDistributor(selfID)->DoReverse(export_objs, obj_size, 
        *len_import_objs, *import_objs);
}

int Epetra_Distributor_DoPosts ( 
  CT_Epetra_Distributor_ID_t selfID, char * export_objs, 
  int obj_size, int * len_import_objs, char ** import_objs )
{
    return CEpetra::getDistributor(selfID)->DoPosts(export_objs, obj_size, 
        *len_import_objs, *import_objs);
}

int Epetra_Distributor_DoWaits ( CT_Epetra_Distributor_ID_t selfID )
{
    return CEpetra::getDistributor(selfID)->DoWaits();
}

int Epetra_Distributor_DoReversePosts ( 
  CT_Epetra_Distributor_ID_t selfID, char * export_objs, 
  int obj_size, int * len_import_objs, char ** import_objs )
{
    return CEpetra::getDistributor(selfID)->DoReversePosts(export_objs, 
        obj_size, *len_import_objs, *import_objs);
}

int Epetra_Distributor_DoReverseWaits ( 
  CT_Epetra_Distributor_ID_t selfID )
{
    return CEpetra::getDistributor(selfID)->DoReverseWaits();
}

int Epetra_Distributor_Do_VarLen ( 
  CT_Epetra_Distributor_ID_t selfID, char * export_objs, 
  int obj_size, int ** sizes, int * len_import_objs, 
  char ** import_objs )
{
    return CEpetra::getDistributor(selfID)->Do(export_objs, obj_size, *sizes, 
        *len_import_objs, *import_objs);
}

int Epetra_Distributor_DoReverse_VarLen ( 
  CT_Epetra_Distributor_ID_t selfID, char * export_objs, 
  int obj_size, int ** sizes, int * len_import_objs, 
  char ** import_objs )
{
    return CEpetra::getDistributor(selfID)->DoReverse(export_objs, obj_size, 
        *sizes, *len_import_objs, *import_objs);
}

int Epetra_Distributor_DoPosts_VarLen ( 
  CT_Epetra_Distributor_ID_t selfID, char * export_objs, 
  int obj_size, int ** sizes, int * len_import_objs, 
  char ** import_objs )
{
    return CEpetra::getDistributor(selfID)->DoPosts(export_objs, obj_size, 
        *sizes, *len_import_objs, *import_objs);
}

int Epetra_Distributor_DoReversePosts_VarLen ( 
  CT_Epetra_Distributor_ID_t selfID, char * export_objs, 
  int obj_size, int ** sizes, int * len_import_objs, 
  char ** import_objs )
{
    return CEpetra::getDistributor(selfID)->DoReversePosts(export_objs, 
        obj_size, *sizes, *len_import_objs, *import_objs);
}


} // extern "C"




