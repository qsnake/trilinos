
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


/*! @file CEpetra_Import.h
 * @brief Wrappers for Epetra_Import */

/* True C header file! */


#ifndef CEPETRA_IMPORT_H
#define CEPETRA_IMPORT_H


#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif



/*! @name ID struct conversion functions */
/*@{*/

/*! @brief Changes the ID struct from the universal
   (generalized) struct type to the class-specific one.
*/
CT_Epetra_Import_ID_t Epetra_Import_Degeneralize ( 
  CTrilinos_Universal_ID_t id );

/*! @brief Changes the ID struct from the class-specific
   struct type to the universal (generalized) one.
*/
CTrilinos_Universal_ID_t Epetra_Import_Generalize ( 
  CT_Epetra_Import_ID_t id );

/*@}*/

/*! @name Epetra_Import constructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   Epetra_Import::Epetra_Import( const Epetra_BlockMap & TargetMap, const Epetra_BlockMap & SourceMap )
*/
CT_Epetra_Import_ID_t Epetra_Import_Create ( 
  CT_Epetra_BlockMap_ID_t TargetMapID, 
  CT_Epetra_BlockMap_ID_t SourceMapID );

/*! @brief Wrapper for 
   Epetra_Import::Epetra_Import(const Epetra_Import& Importer)
*/
CT_Epetra_Import_ID_t Epetra_Import_Duplicate ( 
  CT_Epetra_Import_ID_t ImporterID );

/*@}*/

/*! @name Epetra_Import destructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   virtual Epetra_Import::~Epetra_Import(void)
*/
void Epetra_Import_Destroy ( CT_Epetra_Import_ID_t * selfID );

/*@}*/

/*! @name Epetra_Import member wrappers */
/*@{*/

/*! @brief Wrapper for 
   int Epetra_Import::NumSameIDs() const
*/
int Epetra_Import_NumSameIDs ( CT_Epetra_Import_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_Import::NumPermuteIDs() const
*/
int Epetra_Import_NumPermuteIDs ( CT_Epetra_Import_ID_t selfID );

/*! @brief Wrapper for 
   int * Epetra_Import::PermuteFromLIDs() const
*/
int * Epetra_Import_PermuteFromLIDs ( CT_Epetra_Import_ID_t selfID );

/*! @brief Wrapper for 
   int * Epetra_Import::PermuteToLIDs() const
*/
int * Epetra_Import_PermuteToLIDs ( CT_Epetra_Import_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_Import::NumRemoteIDs() const
*/
int Epetra_Import_NumRemoteIDs ( CT_Epetra_Import_ID_t selfID );

/*! @brief Wrapper for 
   int * Epetra_Import::RemoteLIDs() const
*/
int * Epetra_Import_RemoteLIDs ( CT_Epetra_Import_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_Import::NumExportIDs() const
*/
int Epetra_Import_NumExportIDs ( CT_Epetra_Import_ID_t selfID );

/*! @brief Wrapper for 
   int * Epetra_Import::ExportLIDs() const
*/
int * Epetra_Import_ExportLIDs ( CT_Epetra_Import_ID_t selfID );

/*! @brief Wrapper for 
   int * Epetra_Import::ExportPIDs() const
*/
int * Epetra_Import_ExportPIDs ( CT_Epetra_Import_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_Import::NumSend() const
*/
int Epetra_Import_NumSend ( CT_Epetra_Import_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_Import::NumRecv() const
*/
int Epetra_Import_NumRecv ( CT_Epetra_Import_ID_t selfID );

/*! @brief Wrapper for 
   const Epetra_BlockMap & Epetra_Import::SourceMap() const
*/
CT_Epetra_BlockMap_ID_t Epetra_Import_SourceMap ( 
  CT_Epetra_Import_ID_t selfID );

/*! @brief Wrapper for 
   const Epetra_BlockMap & Epetra_Import::TargetMap() const
*/
CT_Epetra_BlockMap_ID_t Epetra_Import_TargetMap ( 
  CT_Epetra_Import_ID_t selfID );

/*! @brief Wrapper for 
   Epetra_Distributor & Epetra_Import::Distributor() const
*/
CT_Epetra_Distributor_ID_t Epetra_Import_Distributor ( 
  CT_Epetra_Import_ID_t selfID );

/*@}*/


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CEPETRA_IMPORT_H */

