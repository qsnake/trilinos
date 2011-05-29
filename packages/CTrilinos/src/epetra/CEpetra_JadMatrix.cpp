
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
#include "CEpetra_JadMatrix.h"
#include "CEpetra_JadMatrix_Cpp.hpp"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CEpetra_RowMatrix_Cpp.hpp"
#include "CEpetra_MultiVector_Cpp.hpp"


//
// Definitions from CEpetra_JadMatrix.h
//


extern "C" {


CT_Epetra_JadMatrix_ID_t Epetra_JadMatrix_Degeneralize ( 
  CTrilinos_Universal_ID_t id )
{
    return CTrilinos::concreteType<CT_Epetra_JadMatrix_ID_t>(id);
}

CTrilinos_Universal_ID_t Epetra_JadMatrix_Generalize ( 
  CT_Epetra_JadMatrix_ID_t id )
{
    return CTrilinos::abstractType<CT_Epetra_JadMatrix_ID_t>(id);
}

CT_Epetra_JadMatrix_ID_t Epetra_JadMatrix_Create ( 
  CT_Epetra_RowMatrix_ID_t MatrixID )
{
    const Teuchos::RCP<const Epetra_RowMatrix> Matrix = 
        CEpetra::getConstRowMatrix(MatrixID);
    return CEpetra::storeNewJadMatrix(new Epetra_JadMatrix(*Matrix));
}

void Epetra_JadMatrix_Destroy ( CT_Epetra_JadMatrix_ID_t * selfID )
{
    CEpetra::removeJadMatrix(selfID);
}

int Epetra_JadMatrix_UpdateValues ( 
  CT_Epetra_JadMatrix_ID_t selfID, 
  CT_Epetra_RowMatrix_ID_t MatrixID, boolean CheckStructure )
{
    const Teuchos::RCP<const Epetra_RowMatrix> Matrix = 
        CEpetra::getConstRowMatrix(MatrixID);
    return CEpetra::getJadMatrix(selfID)->UpdateValues(*Matrix, ((
        CheckStructure) != FALSE ? true : false));
}

int Epetra_JadMatrix_ExtractMyRowCopy ( 
  CT_Epetra_JadMatrix_ID_t selfID, int MyRow, int Length, 
  int * NumEntries, double * Values, int * Indices )
{
    return CEpetra::getConstJadMatrix(selfID)->ExtractMyRowCopy(MyRow, Length, 
        *NumEntries, Values, Indices);
}

int Epetra_JadMatrix_ExtractMyEntryView ( 
  CT_Epetra_JadMatrix_ID_t selfID, int CurEntry, double * * Value, 
  int * RowIndex, int * ColIndex )
{
    return CEpetra::getJadMatrix(selfID)->ExtractMyEntryView(CurEntry, *Value, 
        *RowIndex, *ColIndex);
}

int Epetra_JadMatrix_ExtractMyEntryView_Const ( 
  CT_Epetra_JadMatrix_ID_t selfID, int CurEntry, 
  double const ** Value, int * RowIndex, int * ColIndex )
{
    return CEpetra::getConstJadMatrix(selfID)->ExtractMyEntryView(CurEntry, 
        *Value, *RowIndex, *ColIndex);
}

int Epetra_JadMatrix_NumMyRowEntries ( 
  CT_Epetra_JadMatrix_ID_t selfID, int MyRow, int * NumEntries )
{
    return CEpetra::getConstJadMatrix(selfID)->NumMyRowEntries(MyRow, 
        *NumEntries);
}

int Epetra_JadMatrix_Multiply ( 
  CT_Epetra_JadMatrix_ID_t selfID, boolean TransA, 
  CT_Epetra_MultiVector_ID_t XID, CT_Epetra_MultiVector_ID_t YID )
{
    const Teuchos::RCP<const Epetra_MultiVector> X = 
        CEpetra::getConstMultiVector(XID);
    const Teuchos::RCP<Epetra_MultiVector> Y = CEpetra::getMultiVector(YID);
    return CEpetra::getConstJadMatrix(selfID)->Multiply(((TransA) != 
        FALSE ? true : false), *X, *Y);
}

int Epetra_JadMatrix_Solve ( 
  CT_Epetra_JadMatrix_ID_t selfID, boolean Upper, boolean Trans, 
  boolean UnitDiagonal, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t YID )
{
    const Teuchos::RCP<const Epetra_MultiVector> X = 
        CEpetra::getConstMultiVector(XID);
    const Teuchos::RCP<Epetra_MultiVector> Y = CEpetra::getMultiVector(YID);
    return CEpetra::getConstJadMatrix(selfID)->Solve(((Upper) != 
        FALSE ? true : false), ((Trans) != FALSE ? true : false), ((
        UnitDiagonal) != FALSE ? true : false), *X, *Y);
}


} // extern "C"




