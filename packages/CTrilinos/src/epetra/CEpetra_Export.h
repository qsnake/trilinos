
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


/*! @file CEpetra_Export.h
 * @brief Wrappers for Epetra_Export */

/* True C header file! */


#ifndef CEPETRA_EXPORT_H
#define CEPETRA_EXPORT_H


#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif



/*! @name ID struct conversion functions */
/*@{*/

/*! @brief Changes the ID struct from the universal
   (generalized) struct type to the class-specific one.
*/
CT_Epetra_Export_ID_t Epetra_Export_Degeneralize ( 
  CTrilinos_Universal_ID_t id );

/*! @brief Changes the ID struct from the class-specific
   struct type to the universal (generalized) one.
*/
CTrilinos_Universal_ID_t Epetra_Export_Generalize ( 
  CT_Epetra_Export_ID_t id );

/*@}*/

/*! @name Epetra_Export constructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   Epetra_Export::Epetra_Export( const Epetra_BlockMap & SourceMap, const Epetra_BlockMap & TargetMap )
*/
CT_Epetra_Export_ID_t Epetra_Export_Create ( 
  CT_Epetra_BlockMap_ID_t SourceMapID, 
  CT_Epetra_BlockMap_ID_t TargetMapID );

/*! @brief Wrapper for 
   Epetra_Export::Epetra_Export(const Epetra_Export& Exporter)
*/
CT_Epetra_Export_ID_t Epetra_Export_Duplicate ( 
  CT_Epetra_Export_ID_t ExporterID );

/*@}*/

/*! @name Epetra_Export destructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   virtual Epetra_Export::~Epetra_Export(void)
*/
void Epetra_Export_Destroy ( CT_Epetra_Export_ID_t * selfID );

/*@}*/

/*! @name Epetra_Export member wrappers */
/*@{*/

/*! @brief Wrapper for 
   int Epetra_Export::NumSameIDs() const
*/
int Epetra_Export_NumSameIDs ( CT_Epetra_Export_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_Export::NumPermuteIDs() const
*/
int Epetra_Export_NumPermuteIDs ( CT_Epetra_Export_ID_t selfID );

/*! @brief Wrapper for 
   int * Epetra_Export::PermuteFromLIDs() const
*/
int * Epetra_Export_PermuteFromLIDs ( CT_Epetra_Export_ID_t selfID );

/*! @brief Wrapper for 
   int * Epetra_Export::PermuteToLIDs() const
*/
int * Epetra_Export_PermuteToLIDs ( CT_Epetra_Export_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_Export::NumRemoteIDs() const
*/
int Epetra_Export_NumRemoteIDs ( CT_Epetra_Export_ID_t selfID );

/*! @brief Wrapper for 
   int * Epetra_Export::RemoteLIDs() const
*/
int * Epetra_Export_RemoteLIDs ( CT_Epetra_Export_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_Export::NumExportIDs() const
*/
int Epetra_Export_NumExportIDs ( CT_Epetra_Export_ID_t selfID );

/*! @brief Wrapper for 
   int * Epetra_Export::ExportLIDs() const
*/
int * Epetra_Export_ExportLIDs ( CT_Epetra_Export_ID_t selfID );

/*! @brief Wrapper for 
   int * Epetra_Export::ExportPIDs() const
*/
int * Epetra_Export_ExportPIDs ( CT_Epetra_Export_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_Export::NumSend() const
*/
int Epetra_Export_NumSend ( CT_Epetra_Export_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_Export::NumRecv() const
*/
int Epetra_Export_NumRecv ( CT_Epetra_Export_ID_t selfID );

/*! @brief Wrapper for 
   const Epetra_BlockMap & Epetra_Export::SourceMap() const
*/
CT_Epetra_BlockMap_ID_t Epetra_Export_SourceMap ( 
  CT_Epetra_Export_ID_t selfID );

/*! @brief Wrapper for 
   const Epetra_BlockMap & Epetra_Export::TargetMap() const
*/
CT_Epetra_BlockMap_ID_t Epetra_Export_TargetMap ( 
  CT_Epetra_Export_ID_t selfID );

/*! @brief Wrapper for 
   Epetra_Distributor & Epetra_Export::Distributor() const
*/
CT_Epetra_Distributor_ID_t Epetra_Export_Distributor ( 
  CT_Epetra_Export_ID_t selfID );

/*@}*/


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CEPETRA_EXPORT_H */

