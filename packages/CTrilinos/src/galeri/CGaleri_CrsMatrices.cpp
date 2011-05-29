
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


#ifdef HAVE_CTRILINOS_GALERI


#include "CTrilinos_enums.h"
#include "CGaleri_CrsMatrices.h"
#include "CGaleri_CrsMatrices_Cpp.hpp"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CEpetra_CrsMatrix_Cpp.hpp"
#include "CEpetra_Map_Cpp.hpp"
#include "CTeuchos_ParameterList_Cpp.hpp"


//
// Definitions from CGaleri_CrsMatrices.h
//


extern "C" {


CT_Epetra_CrsMatrix_ID_t Galeri_CrsMatrices_CreateCrsMatrix ( 
  char MatrixType[], CT_Epetra_Map_ID_t MapID, 
  CT_Teuchos_ParameterList_ID_t ListID )
{
    const Teuchos::RCP<const Epetra_Map> Map = CEpetra::getConstMap(MapID);
    const Teuchos::RCP<Teuchos::ParameterList> List = 
        CTeuchos::getParameterList(ListID);
    return CEpetra::storeCrsMatrix(Galeri::CreateCrsMatrix(std::string(
        MatrixType), Map.getRawPtr(), *List));
}


} // extern "C"


//
// Definitions from CGaleri_CrsMatrices_Cpp.hpp
//




#endif /* HAVE_CTRILINOS_GALERI */


