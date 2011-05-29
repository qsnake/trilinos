
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


/*! @file CEpetra_OffsetIndex.h
 * @brief Wrappers for Epetra_OffsetIndex */

/* True C header file! */


#ifndef CEPETRA_OFFSETINDEX_H
#define CEPETRA_OFFSETINDEX_H


#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif



/*! @name ID struct conversion functions */
/*@{*/

/*! @brief Changes the ID struct from the universal
   (generalized) struct type to the class-specific one.
*/
CT_Epetra_OffsetIndex_ID_t Epetra_OffsetIndex_Degeneralize ( 
  CTrilinos_Universal_ID_t id );

/*! @brief Changes the ID struct from the class-specific
   struct type to the universal (generalized) one.
*/
CTrilinos_Universal_ID_t Epetra_OffsetIndex_Generalize ( 
  CT_Epetra_OffsetIndex_ID_t id );

/*@}*/

/*! @name Epetra_OffsetIndex constructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   Epetra_OffsetIndex::Epetra_OffsetIndex( const Epetra_CrsGraph & SourceGraph, const Epetra_CrsGraph & TargetGraph, Epetra_Import & Importer )
*/
CT_Epetra_OffsetIndex_ID_t Epetra_OffsetIndex_Create_FromImporter ( 
  CT_Epetra_CrsGraph_ID_t SourceGraphID, 
  CT_Epetra_CrsGraph_ID_t TargetGraphID, 
  CT_Epetra_Import_ID_t ImporterID );

/*! @brief Wrapper for 
   Epetra_OffsetIndex::Epetra_OffsetIndex( const Epetra_CrsGraph & SourceGraph, const Epetra_CrsGraph & TargetGraph, Epetra_Export & Exporter )
*/
CT_Epetra_OffsetIndex_ID_t Epetra_OffsetIndex_Create_FromExporter ( 
  CT_Epetra_CrsGraph_ID_t SourceGraphID, 
  CT_Epetra_CrsGraph_ID_t TargetGraphID, 
  CT_Epetra_Export_ID_t ExporterID );

/*! @brief Wrapper for 
   Epetra_OffsetIndex::Epetra_OffsetIndex(const Epetra_OffsetIndex & Indexor)
*/
CT_Epetra_OffsetIndex_ID_t Epetra_OffsetIndex_Duplicate ( 
  CT_Epetra_OffsetIndex_ID_t IndexorID );

/*@}*/

/*! @name Epetra_OffsetIndex destructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   virtual Epetra_OffsetIndex::~Epetra_OffsetIndex(void)
*/
void Epetra_OffsetIndex_Destroy ( 
  CT_Epetra_OffsetIndex_ID_t * selfID );

/*@}*/

/*! @name Epetra_OffsetIndex member wrappers */
/*@{*/

/*! @brief Wrapper for 
   int ** Epetra_OffsetIndex::SameOffsets() const
*/
int ** Epetra_OffsetIndex_SameOffsets ( 
  CT_Epetra_OffsetIndex_ID_t selfID );

/*! @brief Wrapper for 
   int ** Epetra_OffsetIndex::PermuteOffsets() const
*/
int ** Epetra_OffsetIndex_PermuteOffsets ( 
  CT_Epetra_OffsetIndex_ID_t selfID );

/*! @brief Wrapper for 
   int ** Epetra_OffsetIndex::RemoteOffsets() const
*/
int ** Epetra_OffsetIndex_RemoteOffsets ( 
  CT_Epetra_OffsetIndex_ID_t selfID );

/*@}*/


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CEPETRA_OFFSETINDEX_H */

