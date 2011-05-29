
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


/*! @file CEpetra_DistObject.h
 * @brief Wrappers for Epetra_DistObject */

/* True C header file! */


#ifndef CEPETRA_DISTOBJECT_H
#define CEPETRA_DISTOBJECT_H


#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif



/*! @name ID struct conversion functions */
/*@{*/

/*! @brief Changes the ID struct from the universal
   (generalized) struct type to the class-specific one.
*/
CT_Epetra_DistObject_ID_t Epetra_DistObject_Degeneralize ( 
  CTrilinos_Universal_ID_t id );

/*! @brief Changes the ID struct from the class-specific
   struct type to the universal (generalized) one.
*/
CTrilinos_Universal_ID_t Epetra_DistObject_Generalize ( 
  CT_Epetra_DistObject_ID_t id );

/*@}*/

/*! @name Epetra_DistObject destructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   virtual Epetra_DistObject::~Epetra_DistObject()
*/
void Epetra_DistObject_Destroy ( CT_Epetra_DistObject_ID_t * selfID );

/*@}*/

/*! @name Epetra_DistObject member wrappers */
/*@{*/

/*! @brief Wrapper for 
   int Epetra_DistObject::Import(const Epetra_SrcDistObject& A, const Epetra_Import& Importer, Epetra_CombineMode CombineMode, const Epetra_OffsetIndex * Indexor = 0)
*/
int Epetra_DistObject_Import ( 
  CT_Epetra_DistObject_ID_t selfID, 
  CT_Epetra_SrcDistObject_ID_t AID, 
  CT_Epetra_Import_ID_t ImporterID, 
  CT_Epetra_CombineMode_E_t CombineMode, 
  CT_Epetra_OffsetIndex_ID_t IndexorID );

/*! @brief Wrapper for 
   int Epetra_DistObject::Import(const Epetra_SrcDistObject& A, const Epetra_Export& Exporter, Epetra_CombineMode CombineMode, const Epetra_OffsetIndex * Indexor = 0)
*/
int Epetra_DistObject_Import_UsingExporter ( 
  CT_Epetra_DistObject_ID_t selfID, 
  CT_Epetra_SrcDistObject_ID_t AID, 
  CT_Epetra_Export_ID_t ExporterID, 
  CT_Epetra_CombineMode_E_t CombineMode, 
  CT_Epetra_OffsetIndex_ID_t IndexorID );

/*! @brief Wrapper for 
   int Epetra_DistObject::Export(const Epetra_SrcDistObject& A, const Epetra_Import & Importer, Epetra_CombineMode CombineMode, const Epetra_OffsetIndex * Indexor = 0)
*/
int Epetra_DistObject_Export_UsingImporter ( 
  CT_Epetra_DistObject_ID_t selfID, 
  CT_Epetra_SrcDistObject_ID_t AID, 
  CT_Epetra_Import_ID_t ImporterID, 
  CT_Epetra_CombineMode_E_t CombineMode, 
  CT_Epetra_OffsetIndex_ID_t IndexorID );

/*! @brief Wrapper for 
   int Epetra_DistObject::Export(const Epetra_SrcDistObject& A, const Epetra_Export& Exporter, Epetra_CombineMode CombineMode, const Epetra_OffsetIndex * Indexor = 0)
*/
int Epetra_DistObject_Export ( 
  CT_Epetra_DistObject_ID_t selfID, 
  CT_Epetra_SrcDistObject_ID_t AID, 
  CT_Epetra_Export_ID_t ExporterID, 
  CT_Epetra_CombineMode_E_t CombineMode, 
  CT_Epetra_OffsetIndex_ID_t IndexorID );

/*! @brief Wrapper for 
   const Epetra_BlockMap& Epetra_DistObject::Map() const
*/
CT_Epetra_BlockMap_ID_t Epetra_DistObject_Map ( 
  CT_Epetra_DistObject_ID_t selfID );

/*! @brief Wrapper for 
   const Epetra_Comm& Epetra_DistObject::Comm() const
*/
CT_Epetra_Comm_ID_t Epetra_DistObject_Comm ( 
  CT_Epetra_DistObject_ID_t selfID );

/*! @brief Wrapper for 
   bool Epetra_DistObject::DistributedGlobal() const
*/
boolean Epetra_DistObject_DistributedGlobal ( 
  CT_Epetra_DistObject_ID_t selfID );

/*@}*/


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CEPETRA_DISTOBJECT_H */

