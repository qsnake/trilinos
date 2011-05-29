
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


/*! @file CEpetra_Distributor.h
 * @brief Wrappers for Epetra_Distributor */

/* True C header file! */


#ifndef CEPETRA_DISTRIBUTOR_H
#define CEPETRA_DISTRIBUTOR_H


#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif



/*! @name ID struct conversion functions */
/*@{*/

/*! @brief Changes the ID struct from the universal
   (generalized) struct type to the class-specific one.
*/
CT_Epetra_Distributor_ID_t Epetra_Distributor_Degeneralize ( 
  CTrilinos_Universal_ID_t id );

/*! @brief Changes the ID struct from the class-specific
   struct type to the universal (generalized) one.
*/
CTrilinos_Universal_ID_t Epetra_Distributor_Generalize ( 
  CT_Epetra_Distributor_ID_t id );

/*@}*/

/*! @name Epetra_Distributor destructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   virtual Epetra_Distributor::~Epetra_Distributor()
*/
void Epetra_Distributor_Destroy ( 
  CT_Epetra_Distributor_ID_t * selfID );

/*@}*/

/*! @name Epetra_Distributor member wrappers */
/*@{*/

/*! @brief Wrapper for 
   virtual Epetra_Distributor * Epetra_Distributor::Clone() = 0
*/
CT_Epetra_Distributor_ID_t Epetra_Distributor_Clone ( 
  CT_Epetra_Distributor_ID_t selfID );

/*! @brief Wrapper for 
   virtual int Epetra_Distributor::CreateFromSends( const int & NumExportIDs, const int * ExportPIDs, bool Deterministic, int & NumRemoteIDs ) = 0
*/
int Epetra_Distributor_CreateFromSends ( 
  CT_Epetra_Distributor_ID_t selfID, int NumExportIDs, 
  const int * ExportPIDs, boolean Deterministic, 
  int * NumRemoteIDs );

/*! @brief Wrapper for 
   virtual int Epetra_Distributor::CreateFromRecvs( const int & NumRemoteIDs, const int * RemoteGIDs, const int * RemotePIDs, bool Deterministic, int & NumExportIDs, int *& ExportGIDs, int *& ExportPIDs) = 0
*/
int Epetra_Distributor_CreateFromRecvs ( 
  CT_Epetra_Distributor_ID_t selfID, int NumRemoteIDs, 
  const int * RemoteGIDs, const int * RemotePIDs, 
  boolean Deterministic, int * NumExportIDs, int ** ExportGIDs, 
  int ** ExportPIDs );

/*! @brief Wrapper for 
   virtual int Epetra_Distributor::Do( char * export_objs, int obj_size, int & len_import_objs, char *& import_objs) = 0
*/
int Epetra_Distributor_Do ( 
  CT_Epetra_Distributor_ID_t selfID, char * export_objs, 
  int obj_size, int * len_import_objs, char ** import_objs );

/*! @brief Wrapper for 
   virtual int Epetra_Distributor::DoReverse( char * export_objs, int obj_size, int & len_import_objs, char *& import_objs ) = 0
*/
int Epetra_Distributor_DoReverse ( 
  CT_Epetra_Distributor_ID_t selfID, char * export_objs, 
  int obj_size, int * len_import_objs, char ** import_objs );

/*! @brief Wrapper for 
   virtual int Epetra_Distributor::DoPosts( char * export_objs, int obj_size, int & len_import_objs, char *& import_objs ) = 0
*/
int Epetra_Distributor_DoPosts ( 
  CT_Epetra_Distributor_ID_t selfID, char * export_objs, 
  int obj_size, int * len_import_objs, char ** import_objs );

/*! @brief Wrapper for 
   virtual int Epetra_Distributor::DoWaits() = 0
*/
int Epetra_Distributor_DoWaits ( CT_Epetra_Distributor_ID_t selfID );

/*! @brief Wrapper for 
   virtual int Epetra_Distributor::DoReversePosts( char * export_objs, int obj_size, int & len_import_objs, char *& import_objs) = 0
*/
int Epetra_Distributor_DoReversePosts ( 
  CT_Epetra_Distributor_ID_t selfID, char * export_objs, 
  int obj_size, int * len_import_objs, char ** import_objs );

/*! @brief Wrapper for 
   virtual int Epetra_Distributor::DoReverseWaits() = 0
*/
int Epetra_Distributor_DoReverseWaits ( 
  CT_Epetra_Distributor_ID_t selfID );

/*! @brief Wrapper for 
   virtual int Epetra_Distributor::Do( char * export_objs, int obj_size, int *& sizes, int & len_import_objs, char *& import_objs) = 0
*/
int Epetra_Distributor_Do_VarLen ( 
  CT_Epetra_Distributor_ID_t selfID, char * export_objs, 
  int obj_size, int ** sizes, int * len_import_objs, 
  char ** import_objs );

/*! @brief Wrapper for 
   virtual int Epetra_Distributor::DoReverse( char * export_objs, int obj_size, int *& sizes, int & len_import_objs, char *& import_objs) = 0
*/
int Epetra_Distributor_DoReverse_VarLen ( 
  CT_Epetra_Distributor_ID_t selfID, char * export_objs, 
  int obj_size, int ** sizes, int * len_import_objs, 
  char ** import_objs );

/*! @brief Wrapper for 
   virtual int Epetra_Distributor::DoPosts( char * export_objs, int obj_size, int *& sizes, int & len_import_objs, char *& import_objs) = 0
*/
int Epetra_Distributor_DoPosts_VarLen ( 
  CT_Epetra_Distributor_ID_t selfID, char * export_objs, 
  int obj_size, int ** sizes, int * len_import_objs, 
  char ** import_objs );

/*! @brief Wrapper for 
   virtual int Epetra_Distributor::DoReversePosts( char * export_objs, int obj_size, int *& sizes, int & len_import_objs, char *& import_objs) = 0
*/
int Epetra_Distributor_DoReversePosts_VarLen ( 
  CT_Epetra_Distributor_ID_t selfID, char * export_objs, 
  int obj_size, int ** sizes, int * len_import_objs, 
  char ** import_objs );

/*@}*/


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CEPETRA_DISTRIBUTOR_H */

