
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
#include "CEpetra_CrsGraph.h"
#include "CEpetra_CrsGraph_Cpp.hpp"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CEpetra_BlockMap_Cpp.hpp"
#include "CEpetra_Import_Cpp.hpp"
#include "CEpetra_Export_Cpp.hpp"
#include "CEpetra_Comm_Cpp.hpp"


//
// Definitions from CEpetra_CrsGraph.h
//


extern "C" {


CT_Epetra_CrsGraph_ID_t Epetra_CrsGraph_Degeneralize ( 
  CTrilinos_Universal_ID_t id )
{
    return CTrilinos::concreteType<CT_Epetra_CrsGraph_ID_t>(id);
}

CTrilinos_Universal_ID_t Epetra_CrsGraph_Generalize ( 
  CT_Epetra_CrsGraph_ID_t id )
{
    return CTrilinos::abstractType<CT_Epetra_CrsGraph_ID_t>(id);
}

CT_Epetra_CrsGraph_ID_t Epetra_CrsGraph_Create_VarPerRow ( 
  CT_Epetra_DataAccess_E_t CV, CT_Epetra_BlockMap_ID_t RowMapID, 
  const int * NumIndicesPerRow, boolean StaticProfile )
{
    const Teuchos::RCP<const Epetra_BlockMap> RowMap = 
        CEpetra::getConstBlockMap(RowMapID);
    return CEpetra::storeNewCrsGraph(new Epetra_CrsGraph(
        (Epetra_DataAccess) CV, *RowMap, NumIndicesPerRow, ((StaticProfile) != 
        FALSE ? true : false)));
}

CT_Epetra_CrsGraph_ID_t Epetra_CrsGraph_Create ( 
  CT_Epetra_DataAccess_E_t CV, CT_Epetra_BlockMap_ID_t RowMapID, 
  int NumIndicesPerRow, boolean StaticProfile )
{
    const Teuchos::RCP<const Epetra_BlockMap> RowMap = 
        CEpetra::getConstBlockMap(RowMapID);
    return CEpetra::storeNewCrsGraph(new Epetra_CrsGraph(
        (Epetra_DataAccess) CV, *RowMap, NumIndicesPerRow, ((StaticProfile) != 
        FALSE ? true : false)));
}

CT_Epetra_CrsGraph_ID_t Epetra_CrsGraph_Create_VarPerRow_WithColMap ( 
  CT_Epetra_DataAccess_E_t CV, CT_Epetra_BlockMap_ID_t RowMapID, 
  CT_Epetra_BlockMap_ID_t ColMapID, const int * NumIndicesPerRow, 
  boolean StaticProfile )
{
    const Teuchos::RCP<const Epetra_BlockMap> RowMap = 
        CEpetra::getConstBlockMap(RowMapID);
    const Teuchos::RCP<const Epetra_BlockMap> ColMap = 
        CEpetra::getConstBlockMap(ColMapID);
    return CEpetra::storeNewCrsGraph(new Epetra_CrsGraph(
        (Epetra_DataAccess) CV, *RowMap, *ColMap, NumIndicesPerRow, ((
        StaticProfile) != FALSE ? true : false)));
}

CT_Epetra_CrsGraph_ID_t Epetra_CrsGraph_Create_With_ColMap ( 
  CT_Epetra_DataAccess_E_t CV, CT_Epetra_BlockMap_ID_t RowMapID, 
  CT_Epetra_BlockMap_ID_t ColMapID, int NumIndicesPerRow, 
  boolean StaticProfile )
{
    const Teuchos::RCP<const Epetra_BlockMap> RowMap = 
        CEpetra::getConstBlockMap(RowMapID);
    const Teuchos::RCP<const Epetra_BlockMap> ColMap = 
        CEpetra::getConstBlockMap(ColMapID);
    return CEpetra::storeNewCrsGraph(new Epetra_CrsGraph(
        (Epetra_DataAccess) CV, *RowMap, *ColMap, NumIndicesPerRow, ((
        StaticProfile) != FALSE ? true : false)));
}

CT_Epetra_CrsGraph_ID_t Epetra_CrsGraph_Duplicate ( 
  CT_Epetra_CrsGraph_ID_t GraphID )
{
    const Teuchos::RCP<const Epetra_CrsGraph> Graph = CEpetra::getConstCrsGraph(
        GraphID);
    return CEpetra::storeNewCrsGraph(new Epetra_CrsGraph(*Graph));
}

void Epetra_CrsGraph_Destroy ( CT_Epetra_CrsGraph_ID_t * selfID )
{
    CEpetra::removeCrsGraph(selfID);
}

int Epetra_CrsGraph_InsertGlobalIndices ( 
  CT_Epetra_CrsGraph_ID_t selfID, int GlobalRow, int NumIndices, 
  int * Indices )
{
    return CEpetra::getCrsGraph(selfID)->InsertGlobalIndices(GlobalRow, 
        NumIndices, Indices);
}

int Epetra_CrsGraph_RemoveGlobalIndices ( 
  CT_Epetra_CrsGraph_ID_t selfID, int GlobalRow, int NumIndices, 
  int * Indices )
{
    return CEpetra::getCrsGraph(selfID)->RemoveGlobalIndices(GlobalRow, 
        NumIndices, Indices);
}

int Epetra_CrsGraph_RemoveGlobalIndices_LocalRow ( 
  CT_Epetra_CrsGraph_ID_t selfID, int Row )
{
    return CEpetra::getCrsGraph(selfID)->RemoveGlobalIndices(Row);
}

int Epetra_CrsGraph_InsertMyIndices ( 
  CT_Epetra_CrsGraph_ID_t selfID, int LocalRow, int NumIndices, 
  int * Indices )
{
    return CEpetra::getCrsGraph(selfID)->InsertMyIndices(LocalRow, NumIndices, 
        Indices);
}

int Epetra_CrsGraph_RemoveMyIndices ( 
  CT_Epetra_CrsGraph_ID_t selfID, int LocalRow, int NumIndices, 
  int * Indices )
{
    return CEpetra::getCrsGraph(selfID)->RemoveMyIndices(LocalRow, NumIndices, 
        Indices);
}

int Epetra_CrsGraph_RemoveMyIndices_LocalRow ( 
  CT_Epetra_CrsGraph_ID_t selfID, int Row )
{
    return CEpetra::getCrsGraph(selfID)->RemoveMyIndices(Row);
}

int Epetra_CrsGraph_FillComplete ( CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getCrsGraph(selfID)->FillComplete();
}

int Epetra_CrsGraph_FillComplete_UsingMaps ( 
  CT_Epetra_CrsGraph_ID_t selfID, 
  CT_Epetra_BlockMap_ID_t DomainMapID, 
  CT_Epetra_BlockMap_ID_t RangeMapID )
{
    const Teuchos::RCP<const Epetra_BlockMap> DomainMap = 
        CEpetra::getConstBlockMap(DomainMapID);
    const Teuchos::RCP<const Epetra_BlockMap> RangeMap = 
        CEpetra::getConstBlockMap(RangeMapID);
    return CEpetra::getCrsGraph(selfID)->FillComplete(*DomainMap, *RangeMap);
}

int Epetra_CrsGraph_OptimizeStorage ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getCrsGraph(selfID)->OptimizeStorage();
}

int Epetra_CrsGraph_ExtractGlobalRowCopy ( 
  CT_Epetra_CrsGraph_ID_t selfID, int GlobalRow, int LenOfIndices, 
  int * NumIndices, int * Indices )
{
    return CEpetra::getConstCrsGraph(selfID)->ExtractGlobalRowCopy(GlobalRow, 
        LenOfIndices, *NumIndices, Indices);
}

int Epetra_CrsGraph_ExtractMyRowCopy ( 
  CT_Epetra_CrsGraph_ID_t selfID, int LocalRow, int LenOfIndices, 
  int * NumIndices, int * Indices )
{
    return CEpetra::getConstCrsGraph(selfID)->ExtractMyRowCopy(LocalRow, 
        LenOfIndices, *NumIndices, Indices);
}

int Epetra_CrsGraph_ExtractGlobalRowView ( 
  CT_Epetra_CrsGraph_ID_t selfID, int GlobalRow, int * NumIndices, 
  int ** Indices )
{
    return CEpetra::getConstCrsGraph(selfID)->ExtractGlobalRowView(GlobalRow, 
        *NumIndices, *Indices);
}

int Epetra_CrsGraph_ExtractMyRowView ( 
  CT_Epetra_CrsGraph_ID_t selfID, int LocalRow, int * NumIndices, 
  int ** Indices )
{
    return CEpetra::getConstCrsGraph(selfID)->ExtractMyRowView(LocalRow, 
        *NumIndices, *Indices);
}

boolean Epetra_CrsGraph_Filled ( CT_Epetra_CrsGraph_ID_t selfID )
{
    return ((CEpetra::getConstCrsGraph(selfID)->Filled()) ? TRUE : FALSE);
}

boolean Epetra_CrsGraph_StorageOptimized ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return ((CEpetra::getConstCrsGraph(
        selfID)->StorageOptimized()) ? TRUE : FALSE);
}

boolean Epetra_CrsGraph_IndicesAreGlobal ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return ((CEpetra::getConstCrsGraph(
        selfID)->IndicesAreGlobal()) ? TRUE : FALSE);
}

boolean Epetra_CrsGraph_IndicesAreLocal ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return ((CEpetra::getConstCrsGraph(
        selfID)->IndicesAreLocal()) ? TRUE : FALSE);
}

boolean Epetra_CrsGraph_LowerTriangular ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return ((CEpetra::getConstCrsGraph(
        selfID)->LowerTriangular()) ? TRUE : FALSE);
}

boolean Epetra_CrsGraph_UpperTriangular ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return ((CEpetra::getConstCrsGraph(
        selfID)->UpperTriangular()) ? TRUE : FALSE);
}

boolean Epetra_CrsGraph_NoDiagonal ( CT_Epetra_CrsGraph_ID_t selfID )
{
    return ((CEpetra::getConstCrsGraph(selfID)->NoDiagonal()) ? TRUE : FALSE);
}

boolean Epetra_CrsGraph_MyGlobalRow ( 
  CT_Epetra_CrsGraph_ID_t selfID, int GID )
{
    return ((CEpetra::getConstCrsGraph(selfID)->MyGlobalRow(
        GID)) ? TRUE : FALSE);
}

boolean Epetra_CrsGraph_HaveColMap ( CT_Epetra_CrsGraph_ID_t selfID )
{
    return ((CEpetra::getConstCrsGraph(selfID)->HaveColMap()) ? TRUE : FALSE);
}

int Epetra_CrsGraph_NumMyRows ( CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getConstCrsGraph(selfID)->NumMyRows();
}

int Epetra_CrsGraph_NumGlobalRows ( CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getConstCrsGraph(selfID)->NumGlobalRows();
}

int Epetra_CrsGraph_NumMyCols ( CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getConstCrsGraph(selfID)->NumMyCols();
}

int Epetra_CrsGraph_NumGlobalCols ( CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getConstCrsGraph(selfID)->NumGlobalCols();
}

int Epetra_CrsGraph_NumGlobalNonzeros ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getConstCrsGraph(selfID)->NumGlobalNonzeros();
}

int Epetra_CrsGraph_NumGlobalDiagonals ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getConstCrsGraph(selfID)->NumGlobalDiagonals();
}

int Epetra_CrsGraph_NumMyDiagonals ( CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getConstCrsGraph(selfID)->NumMyDiagonals();
}

int Epetra_CrsGraph_NumMyBlockRows ( CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getConstCrsGraph(selfID)->NumMyBlockRows();
}

int Epetra_CrsGraph_NumGlobalBlockRows ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getConstCrsGraph(selfID)->NumGlobalBlockRows();
}

int Epetra_CrsGraph_NumMyBlockCols ( CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getConstCrsGraph(selfID)->NumMyBlockCols();
}

int Epetra_CrsGraph_NumGlobalBlockCols ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getConstCrsGraph(selfID)->NumGlobalBlockCols();
}

int Epetra_CrsGraph_NumMyBlockDiagonals ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getConstCrsGraph(selfID)->NumMyBlockDiagonals();
}

int Epetra_CrsGraph_NumGlobalBlockDiagonals ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getConstCrsGraph(selfID)->NumGlobalBlockDiagonals();
}

int Epetra_CrsGraph_NumGlobalEntries ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getConstCrsGraph(selfID)->NumGlobalEntries();
}

int Epetra_CrsGraph_NumMyEntries ( CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getConstCrsGraph(selfID)->NumMyEntries();
}

int Epetra_CrsGraph_MaxRowDim ( CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getConstCrsGraph(selfID)->MaxRowDim();
}

int Epetra_CrsGraph_GlobalMaxRowDim ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getConstCrsGraph(selfID)->GlobalMaxRowDim();
}

int Epetra_CrsGraph_MaxColDim ( CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getConstCrsGraph(selfID)->MaxColDim();
}

int Epetra_CrsGraph_GlobalMaxColDim ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getConstCrsGraph(selfID)->GlobalMaxColDim();
}

int Epetra_CrsGraph_NumMyNonzeros ( CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getConstCrsGraph(selfID)->NumMyNonzeros();
}

int Epetra_CrsGraph_NumGlobalIndices ( 
  CT_Epetra_CrsGraph_ID_t selfID, int Row )
{
    return CEpetra::getConstCrsGraph(selfID)->NumGlobalIndices(Row);
}

int Epetra_CrsGraph_NumAllocatedGlobalIndices ( 
  CT_Epetra_CrsGraph_ID_t selfID, int Row )
{
    return CEpetra::getConstCrsGraph(selfID)->NumAllocatedGlobalIndices(Row);
}

int Epetra_CrsGraph_MaxNumIndices ( CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getConstCrsGraph(selfID)->MaxNumIndices();
}

int Epetra_CrsGraph_GlobalMaxNumIndices ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getConstCrsGraph(selfID)->GlobalMaxNumIndices();
}

int Epetra_CrsGraph_MaxNumNonzeros ( CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getConstCrsGraph(selfID)->MaxNumNonzeros();
}

int Epetra_CrsGraph_GlobalMaxNumNonzeros ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getConstCrsGraph(selfID)->GlobalMaxNumNonzeros();
}

int Epetra_CrsGraph_NumMyIndices ( 
  CT_Epetra_CrsGraph_ID_t selfID, int Row )
{
    return CEpetra::getConstCrsGraph(selfID)->NumMyIndices(Row);
}

int Epetra_CrsGraph_NumAllocatedMyIndices ( 
  CT_Epetra_CrsGraph_ID_t selfID, int Row )
{
    return CEpetra::getConstCrsGraph(selfID)->NumAllocatedMyIndices(Row);
}

int Epetra_CrsGraph_IndexBase ( CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getConstCrsGraph(selfID)->IndexBase();
}

CT_Epetra_BlockMap_ID_t Epetra_CrsGraph_RowMap ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::storeConstBlockMap(&( CEpetra::getConstCrsGraph(
        selfID)->RowMap() ));
}

int Epetra_CrsGraph_ReplaceRowMap ( 
  CT_Epetra_CrsGraph_ID_t selfID, CT_Epetra_BlockMap_ID_t newmapID )
{
    const Teuchos::RCP<const Epetra_BlockMap> newmap = 
        CEpetra::getConstBlockMap(newmapID);
    return CEpetra::getCrsGraph(selfID)->ReplaceRowMap(*newmap);
}

int Epetra_CrsGraph_ReplaceColMap ( 
  CT_Epetra_CrsGraph_ID_t selfID, CT_Epetra_BlockMap_ID_t newmapID )
{
    const Teuchos::RCP<const Epetra_BlockMap> newmap = 
        CEpetra::getConstBlockMap(newmapID);
    return CEpetra::getCrsGraph(selfID)->ReplaceColMap(*newmap);
}

CT_Epetra_BlockMap_ID_t Epetra_CrsGraph_ColMap ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::storeConstBlockMap(&( CEpetra::getConstCrsGraph(
        selfID)->ColMap() ));
}

CT_Epetra_BlockMap_ID_t Epetra_CrsGraph_DomainMap ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::storeConstBlockMap(&( CEpetra::getConstCrsGraph(
        selfID)->DomainMap() ));
}

CT_Epetra_BlockMap_ID_t Epetra_CrsGraph_RangeMap ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::storeConstBlockMap(&( CEpetra::getConstCrsGraph(
        selfID)->RangeMap() ));
}

CT_Epetra_Import_ID_t Epetra_CrsGraph_Importer ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::storeConstImport(CEpetra::getConstCrsGraph(
        selfID)->Importer());
}

CT_Epetra_Export_ID_t Epetra_CrsGraph_Exporter ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::storeConstExport(CEpetra::getConstCrsGraph(
        selfID)->Exporter());
}

CT_Epetra_Comm_ID_t Epetra_CrsGraph_Comm ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::storeConstComm(&( CEpetra::getConstCrsGraph(
        selfID)->Comm() ));
}

int Epetra_CrsGraph_LRID ( 
  CT_Epetra_CrsGraph_ID_t selfID, int GRID_in )
{
    return CEpetra::getConstCrsGraph(selfID)->LRID(GRID_in);
}

int Epetra_CrsGraph_GRID ( 
  CT_Epetra_CrsGraph_ID_t selfID, int LRID_in )
{
    return CEpetra::getConstCrsGraph(selfID)->GRID(LRID_in);
}

int Epetra_CrsGraph_LCID ( 
  CT_Epetra_CrsGraph_ID_t selfID, int GCID_in )
{
    return CEpetra::getConstCrsGraph(selfID)->LCID(GCID_in);
}

int Epetra_CrsGraph_GCID ( 
  CT_Epetra_CrsGraph_ID_t selfID, int LCID_in )
{
    return CEpetra::getConstCrsGraph(selfID)->GCID(LCID_in);
}

boolean Epetra_CrsGraph_MyGRID ( 
  CT_Epetra_CrsGraph_ID_t selfID, int GRID_in )
{
    return ((CEpetra::getConstCrsGraph(selfID)->MyGRID(
        GRID_in)) ? TRUE : FALSE);
}

boolean Epetra_CrsGraph_MyLRID ( 
  CT_Epetra_CrsGraph_ID_t selfID, int LRID_in )
{
    return ((CEpetra::getConstCrsGraph(selfID)->MyLRID(
        LRID_in)) ? TRUE : FALSE);
}

boolean Epetra_CrsGraph_MyGCID ( 
  CT_Epetra_CrsGraph_ID_t selfID, int GCID_in )
{
    return ((CEpetra::getConstCrsGraph(selfID)->MyGCID(
        GCID_in)) ? TRUE : FALSE);
}

boolean Epetra_CrsGraph_MyLCID ( 
  CT_Epetra_CrsGraph_ID_t selfID, int LCID_in )
{
    return ((CEpetra::getConstCrsGraph(selfID)->MyLCID(
        LCID_in)) ? TRUE : FALSE);
}

CT_Epetra_BlockMap_ID_t Epetra_CrsGraph_ImportMap ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::storeConstBlockMap(&( CEpetra::getConstCrsGraph(
        selfID)->ImportMap() ));
}

int Epetra_CrsGraph_TransformToLocal ( 
  CT_Epetra_CrsGraph_ID_t selfID )
{
    return CEpetra::getCrsGraph(selfID)->TransformToLocal();
}

int Epetra_CrsGraph_TransformToLocal_UsingMaps ( 
  CT_Epetra_CrsGraph_ID_t selfID, 
  CT_Epetra_BlockMap_ID_t DomainMapID, 
  CT_Epetra_BlockMap_ID_t RangeMapID )
{
    const Teuchos::RCP<const Epetra_BlockMap> DomainMap = 
        CEpetra::getConstBlockMap(DomainMapID);
    const Teuchos::RCP<const Epetra_BlockMap> RangeMap = 
        CEpetra::getConstBlockMap(RangeMapID);
    return CEpetra::getCrsGraph(selfID)->TransformToLocal(
        DomainMap.getRawPtr(), RangeMap.getRawPtr());
}

int * Epetra_CrsGraph_getRow ( 
  CT_Epetra_CrsGraph_ID_t selfID, int Loc )
{
    const Epetra_CrsGraph& self = *( CEpetra::getConstCrsGraph(selfID) );

    return self[Loc];
}

void Epetra_CrsGraph_Assign ( 
  CT_Epetra_CrsGraph_ID_t selfID, CT_Epetra_CrsGraph_ID_t SourceID )
{
    Epetra_CrsGraph& self = *( CEpetra::getCrsGraph(selfID) );

    const Teuchos::RCP<const Epetra_CrsGraph> Source = 
        CEpetra::getConstCrsGraph(SourceID);
    self = *Source;
}


} // extern "C"




