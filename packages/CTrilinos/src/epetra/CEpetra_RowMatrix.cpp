
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
#include "CEpetra_RowMatrix.h"
#include "CEpetra_RowMatrix_Cpp.hpp"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CEpetra_Vector_Cpp.hpp"
#include "CEpetra_MultiVector_Cpp.hpp"
#include "CEpetra_Map_Cpp.hpp"
#include "CEpetra_Import_Cpp.hpp"


//
// Definitions from CEpetra_RowMatrix.h
//


extern "C" {


CT_Epetra_RowMatrix_ID_t Epetra_RowMatrix_Degeneralize ( 
  CTrilinos_Universal_ID_t id )
{
    return CTrilinos::concreteType<CT_Epetra_RowMatrix_ID_t>(id);
}

CTrilinos_Universal_ID_t Epetra_RowMatrix_Generalize ( 
  CT_Epetra_RowMatrix_ID_t id )
{
    return CTrilinos::abstractType<CT_Epetra_RowMatrix_ID_t>(id);
}

void Epetra_RowMatrix_Destroy ( CT_Epetra_RowMatrix_ID_t * selfID )
{
    CEpetra::removeRowMatrix(selfID);
}

int Epetra_RowMatrix_NumMyRowEntries ( 
  CT_Epetra_RowMatrix_ID_t selfID, int MyRow, int * NumEntries )
{
    return CEpetra::getConstRowMatrix(selfID)->NumMyRowEntries(MyRow, 
        *NumEntries);
}

int Epetra_RowMatrix_MaxNumEntries ( 
  CT_Epetra_RowMatrix_ID_t selfID )
{
    return CEpetra::getConstRowMatrix(selfID)->MaxNumEntries();
}

int Epetra_RowMatrix_ExtractMyRowCopy ( 
  CT_Epetra_RowMatrix_ID_t selfID, int MyRow, int Length, 
  int * NumEntries, double * Values, int * Indices )
{
    return CEpetra::getConstRowMatrix(selfID)->ExtractMyRowCopy(MyRow, Length, 
        *NumEntries, Values, Indices);
}

int Epetra_RowMatrix_ExtractDiagonalCopy ( 
  CT_Epetra_RowMatrix_ID_t selfID, CT_Epetra_Vector_ID_t DiagonalID )
{
    const Teuchos::RCP<Epetra_Vector> Diagonal = CEpetra::getVector(
        DiagonalID);
    return CEpetra::getConstRowMatrix(selfID)->ExtractDiagonalCopy(*Diagonal);
}

int Epetra_RowMatrix_Multiply ( 
  CT_Epetra_RowMatrix_ID_t selfID, boolean TransA, 
  CT_Epetra_MultiVector_ID_t XID, CT_Epetra_MultiVector_ID_t YID )
{
    const Teuchos::RCP<const Epetra_MultiVector> X = 
        CEpetra::getConstMultiVector(XID);
    const Teuchos::RCP<Epetra_MultiVector> Y = CEpetra::getMultiVector(YID);
    return CEpetra::getConstRowMatrix(selfID)->Multiply(((TransA) != 
        FALSE ? true : false), *X, *Y);
}

int Epetra_RowMatrix_Solve ( 
  CT_Epetra_RowMatrix_ID_t selfID, boolean Upper, boolean Trans, 
  boolean UnitDiagonal, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t YID )
{
    const Teuchos::RCP<const Epetra_MultiVector> X = 
        CEpetra::getConstMultiVector(XID);
    const Teuchos::RCP<Epetra_MultiVector> Y = CEpetra::getMultiVector(YID);
    return CEpetra::getConstRowMatrix(selfID)->Solve(((Upper) != 
        FALSE ? true : false), ((Trans) != FALSE ? true : false), ((
        UnitDiagonal) != FALSE ? true : false), *X, *Y);
}

int Epetra_RowMatrix_InvRowSums ( 
  CT_Epetra_RowMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID )
{
    const Teuchos::RCP<Epetra_Vector> x = CEpetra::getVector(xID);
    return CEpetra::getConstRowMatrix(selfID)->InvRowSums(*x);
}

int Epetra_RowMatrix_LeftScale ( 
  CT_Epetra_RowMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID )
{
    const Teuchos::RCP<const Epetra_Vector> x = CEpetra::getConstVector(xID);
    return CEpetra::getRowMatrix(selfID)->LeftScale(*x);
}

int Epetra_RowMatrix_InvColSums ( 
  CT_Epetra_RowMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID )
{
    const Teuchos::RCP<Epetra_Vector> x = CEpetra::getVector(xID);
    return CEpetra::getConstRowMatrix(selfID)->InvColSums(*x);
}

int Epetra_RowMatrix_RightScale ( 
  CT_Epetra_RowMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID )
{
    const Teuchos::RCP<const Epetra_Vector> x = CEpetra::getConstVector(xID);
    return CEpetra::getRowMatrix(selfID)->RightScale(*x);
}

boolean Epetra_RowMatrix_Filled ( CT_Epetra_RowMatrix_ID_t selfID )
{
    return ((CEpetra::getConstRowMatrix(selfID)->Filled()) ? TRUE : FALSE);
}

double Epetra_RowMatrix_NormInf ( CT_Epetra_RowMatrix_ID_t selfID )
{
    return CEpetra::getConstRowMatrix(selfID)->NormInf();
}

double Epetra_RowMatrix_NormOne ( CT_Epetra_RowMatrix_ID_t selfID )
{
    return CEpetra::getConstRowMatrix(selfID)->NormOne();
}

int Epetra_RowMatrix_NumGlobalNonzeros ( 
  CT_Epetra_RowMatrix_ID_t selfID )
{
    return CEpetra::getConstRowMatrix(selfID)->NumGlobalNonzeros();
}

int Epetra_RowMatrix_NumGlobalRows ( 
  CT_Epetra_RowMatrix_ID_t selfID )
{
    return CEpetra::getConstRowMatrix(selfID)->NumGlobalRows();
}

int Epetra_RowMatrix_NumGlobalCols ( 
  CT_Epetra_RowMatrix_ID_t selfID )
{
    return CEpetra::getConstRowMatrix(selfID)->NumGlobalCols();
}

int Epetra_RowMatrix_NumGlobalDiagonals ( 
  CT_Epetra_RowMatrix_ID_t selfID )
{
    return CEpetra::getConstRowMatrix(selfID)->NumGlobalDiagonals();
}

int Epetra_RowMatrix_NumMyNonzeros ( 
  CT_Epetra_RowMatrix_ID_t selfID )
{
    return CEpetra::getConstRowMatrix(selfID)->NumMyNonzeros();
}

int Epetra_RowMatrix_NumMyRows ( CT_Epetra_RowMatrix_ID_t selfID )
{
    return CEpetra::getConstRowMatrix(selfID)->NumMyRows();
}

int Epetra_RowMatrix_NumMyCols ( CT_Epetra_RowMatrix_ID_t selfID )
{
    return CEpetra::getConstRowMatrix(selfID)->NumMyCols();
}

int Epetra_RowMatrix_NumMyDiagonals ( 
  CT_Epetra_RowMatrix_ID_t selfID )
{
    return CEpetra::getConstRowMatrix(selfID)->NumMyDiagonals();
}

boolean Epetra_RowMatrix_LowerTriangular ( 
  CT_Epetra_RowMatrix_ID_t selfID )
{
    return ((CEpetra::getConstRowMatrix(
        selfID)->LowerTriangular()) ? TRUE : FALSE);
}

boolean Epetra_RowMatrix_UpperTriangular ( 
  CT_Epetra_RowMatrix_ID_t selfID )
{
    return ((CEpetra::getConstRowMatrix(
        selfID)->UpperTriangular()) ? TRUE : FALSE);
}

CT_Epetra_Map_ID_t Epetra_RowMatrix_RowMatrixRowMap ( 
  CT_Epetra_RowMatrix_ID_t selfID )
{
    return CEpetra::storeConstMap(&( CEpetra::getConstRowMatrix(
        selfID)->RowMatrixRowMap() ));
}

CT_Epetra_Map_ID_t Epetra_RowMatrix_RowMatrixColMap ( 
  CT_Epetra_RowMatrix_ID_t selfID )
{
    return CEpetra::storeConstMap(&( CEpetra::getConstRowMatrix(
        selfID)->RowMatrixColMap() ));
}

CT_Epetra_Import_ID_t Epetra_RowMatrix_RowMatrixImporter ( 
  CT_Epetra_RowMatrix_ID_t selfID )
{
    return CEpetra::storeConstImport(CEpetra::getConstRowMatrix(
        selfID)->RowMatrixImporter());
}


} // extern "C"




