
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
#include "CEpetra_SerialDenseMatrix.h"
#include "CEpetra_SerialDenseMatrix_Cpp.hpp"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"


//
// Definitions from CEpetra_SerialDenseMatrix.h
//


extern "C" {


CT_Epetra_SerialDenseMatrix_ID_t Epetra_SerialDenseMatrix_Degeneralize ( 
  CTrilinos_Universal_ID_t id )
{
    return CTrilinos::concreteType<CT_Epetra_SerialDenseMatrix_ID_t>(id);
}

CTrilinos_Universal_ID_t Epetra_SerialDenseMatrix_Generalize ( 
  CT_Epetra_SerialDenseMatrix_ID_t id )
{
    return CTrilinos::abstractType<CT_Epetra_SerialDenseMatrix_ID_t>(id);
}

CT_Epetra_SerialDenseMatrix_ID_t Epetra_SerialDenseMatrix_Create_Empty ( 
  boolean set_object_label )
{
    return CEpetra::storeNewSerialDenseMatrix(new Epetra_SerialDenseMatrix(
        ((set_object_label) != FALSE ? true : false)));
}

CT_Epetra_SerialDenseMatrix_ID_t Epetra_SerialDenseMatrix_Create ( 
  int NumRows, int NumCols, boolean set_object_label )
{
    return CEpetra::storeNewSerialDenseMatrix(new Epetra_SerialDenseMatrix(
        NumRows, NumCols, ((set_object_label) != FALSE ? true : false)));
}

CT_Epetra_SerialDenseMatrix_ID_t Epetra_SerialDenseMatrix_Create_FromArray ( 
  CT_Epetra_DataAccess_E_t CV, double * A_in, int LDA_in, 
  int NumRows, int NumCols, boolean set_object_label )
{
    return CEpetra::storeNewSerialDenseMatrix(new Epetra_SerialDenseMatrix(
        (Epetra_DataAccess) CV, A_in, LDA_in, NumRows, NumCols, ((
        set_object_label) != FALSE ? true : false)));
}

CT_Epetra_SerialDenseMatrix_ID_t Epetra_SerialDenseMatrix_Duplicate ( 
  CT_Epetra_SerialDenseMatrix_ID_t SourceID )
{
    const Teuchos::RCP<const Epetra_SerialDenseMatrix> Source = 
        CEpetra::getConstSerialDenseMatrix(SourceID);
    return CEpetra::storeNewSerialDenseMatrix(new Epetra_SerialDenseMatrix(
        *Source));
}

void Epetra_SerialDenseMatrix_Destroy ( 
  CT_Epetra_SerialDenseMatrix_ID_t * selfID )
{
    CEpetra::removeSerialDenseMatrix(selfID);
}

int Epetra_SerialDenseMatrix_Shape ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID, int NumRows, int NumCols )
{
    return CEpetra::getSerialDenseMatrix(selfID)->Shape(NumRows, NumCols);
}

int Epetra_SerialDenseMatrix_Reshape ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID, int NumRows, int NumCols )
{
    return CEpetra::getSerialDenseMatrix(selfID)->Reshape(NumRows, NumCols);
}

int Epetra_SerialDenseMatrix_Multiply_Matrix ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID, char TransA, char TransB, 
  double ScalarAB, CT_Epetra_SerialDenseMatrix_ID_t AID, 
  CT_Epetra_SerialDenseMatrix_ID_t BID, double ScalarThis )
{
    const Teuchos::RCP<const Epetra_SerialDenseMatrix> A = 
        CEpetra::getConstSerialDenseMatrix(AID);
    const Teuchos::RCP<const Epetra_SerialDenseMatrix> B = 
        CEpetra::getConstSerialDenseMatrix(BID);
    return CEpetra::getSerialDenseMatrix(selfID)->Multiply(TransA, TransB, 
        ScalarAB, *A, *B, ScalarThis);
}

int Epetra_SerialDenseMatrix_Multiply_Vector ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID, boolean transA, 
  CT_Epetra_SerialDenseMatrix_ID_t xID, 
  CT_Epetra_SerialDenseMatrix_ID_t yID )
{
    const Teuchos::RCP<const Epetra_SerialDenseMatrix> x = 
        CEpetra::getConstSerialDenseMatrix(xID);
    const Teuchos::RCP<Epetra_SerialDenseMatrix> y = 
        CEpetra::getSerialDenseMatrix(yID);
    return CEpetra::getSerialDenseMatrix(selfID)->Multiply(((transA) != 
        FALSE ? true : false), *x, *y);
}

int Epetra_SerialDenseMatrix_Scale ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID, double ScalarA )
{
    return CEpetra::getSerialDenseMatrix(selfID)->Scale(ScalarA);
}

double Epetra_SerialDenseMatrix_NormOne ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID )
{
    return CEpetra::getConstSerialDenseMatrix(selfID)->NormOne();
}

double Epetra_SerialDenseMatrix_NormInf ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID )
{
    return CEpetra::getConstSerialDenseMatrix(selfID)->NormInf();
}

int Epetra_SerialDenseMatrix_Random ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID )
{
    return CEpetra::getSerialDenseMatrix(selfID)->Random();
}

int Epetra_SerialDenseMatrix_M ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID )
{
    return CEpetra::getConstSerialDenseMatrix(selfID)->M();
}

int Epetra_SerialDenseMatrix_N ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID )
{
    return CEpetra::getConstSerialDenseMatrix(selfID)->N();
}

double * Epetra_SerialDenseMatrix_A_Const ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID )
{
    return CEpetra::getConstSerialDenseMatrix(selfID)->A();
}

double * Epetra_SerialDenseMatrix_A ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID )
{
    return CEpetra::getSerialDenseMatrix(selfID)->A();
}

int Epetra_SerialDenseMatrix_LDA ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID )
{
    return CEpetra::getConstSerialDenseMatrix(selfID)->LDA();
}

CT_Epetra_DataAccess_E_t Epetra_SerialDenseMatrix_CV ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID )
{
    return (CT_Epetra_DataAccess_E_t)( CEpetra::getConstSerialDenseMatrix(
        selfID)->CV() );
}

double Epetra_SerialDenseMatrix_OneNorm ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID )
{
    return CEpetra::getConstSerialDenseMatrix(selfID)->OneNorm();
}

double Epetra_SerialDenseMatrix_InfNorm ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID )
{
    return CEpetra::getConstSerialDenseMatrix(selfID)->InfNorm();
}

int Epetra_SerialDenseMatrix_SetUseTranspose ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID, boolean UseTranspose_in )
{
    return CEpetra::getSerialDenseMatrix(selfID)->SetUseTranspose(
        ((UseTranspose_in) != FALSE ? true : false));
}

int Epetra_SerialDenseMatrix_Apply ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID, 
  CT_Epetra_SerialDenseMatrix_ID_t XID, 
  CT_Epetra_SerialDenseMatrix_ID_t YID )
{
    const Teuchos::RCP<const Epetra_SerialDenseMatrix> X = 
        CEpetra::getConstSerialDenseMatrix(XID);
    const Teuchos::RCP<Epetra_SerialDenseMatrix> Y = 
        CEpetra::getSerialDenseMatrix(YID);
    return CEpetra::getSerialDenseMatrix(selfID)->Apply(*X, *Y);
}

int Epetra_SerialDenseMatrix_ApplyInverse ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID, 
  CT_Epetra_SerialDenseMatrix_ID_t XID, 
  CT_Epetra_SerialDenseMatrix_ID_t YID )
{
    const Teuchos::RCP<const Epetra_SerialDenseMatrix> X = 
        CEpetra::getConstSerialDenseMatrix(XID);
    const Teuchos::RCP<Epetra_SerialDenseMatrix> Y = 
        CEpetra::getSerialDenseMatrix(YID);
    return CEpetra::getSerialDenseMatrix(selfID)->ApplyInverse(*X, *Y);
}

const char * Epetra_SerialDenseMatrix_Label ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID )
{
    return CEpetra::getConstSerialDenseMatrix(selfID)->Label();
}

boolean Epetra_SerialDenseMatrix_UseTranspose ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID )
{
    return ((CEpetra::getConstSerialDenseMatrix(
        selfID)->UseTranspose()) ? TRUE : FALSE);
}

boolean Epetra_SerialDenseMatrix_HasNormInf ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID )
{
    return ((CEpetra::getConstSerialDenseMatrix(
        selfID)->HasNormInf()) ? TRUE : FALSE);
}

int Epetra_SerialDenseMatrix_RowDim ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID )
{
    return CEpetra::getConstSerialDenseMatrix(selfID)->RowDim();
}

int Epetra_SerialDenseMatrix_ColDim ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID )
{
    return CEpetra::getConstSerialDenseMatrix(selfID)->ColDim();
}

void Epetra_SerialDenseMatrix_Assign ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID, 
  CT_Epetra_SerialDenseMatrix_ID_t SourceID )
{
    Epetra_SerialDenseMatrix& self = *( CEpetra::getSerialDenseMatrix(
        selfID) );

    const Teuchos::RCP<const Epetra_SerialDenseMatrix> Source = 
        CEpetra::getConstSerialDenseMatrix(SourceID);
    self = *Source;
}

boolean Epetra_SerialDenseMatrix_IsEqual ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID, 
  CT_Epetra_SerialDenseMatrix_ID_t rhsID )
{
    const Epetra_SerialDenseMatrix& self = *( 
        CEpetra::getConstSerialDenseMatrix(selfID) );

    const Teuchos::RCP<const Epetra_SerialDenseMatrix> rhs = 
        CEpetra::getConstSerialDenseMatrix(rhsID);
    return ((self == *rhs) ? TRUE : FALSE);
}

boolean Epetra_SerialDenseMatrix_NotEqual ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID, 
  CT_Epetra_SerialDenseMatrix_ID_t rhsID )
{
    const Epetra_SerialDenseMatrix& self = *( 
        CEpetra::getConstSerialDenseMatrix(selfID) );

    const Teuchos::RCP<const Epetra_SerialDenseMatrix> rhs = 
        CEpetra::getConstSerialDenseMatrix(rhsID);
    return ((self != *rhs) ? TRUE : FALSE);
}

void Epetra_SerialDenseMatrix_AddTo ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID, 
  CT_Epetra_SerialDenseMatrix_ID_t SourceID )
{
    Epetra_SerialDenseMatrix& self = *( CEpetra::getSerialDenseMatrix(
        selfID) );

    const Teuchos::RCP<const Epetra_SerialDenseMatrix> Source = 
        CEpetra::getConstSerialDenseMatrix(SourceID);
    self += *Source;
}

void Epetra_SerialDenseMatrix_setElement ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID, int RowIndex, 
  int ColIndex, double * value )
{
    Epetra_SerialDenseMatrix& self = *( CEpetra::getSerialDenseMatrix(
        selfID) );

    self(RowIndex, ColIndex) = *value;
}

double Epetra_SerialDenseMatrix_getElement ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID, int RowIndex, 
  int ColIndex )
{
    const Epetra_SerialDenseMatrix& self = *( 
        CEpetra::getConstSerialDenseMatrix(selfID) );

    return self(RowIndex, ColIndex);
}

const double * Epetra_SerialDenseMatrix_getColumn ( 
  CT_Epetra_SerialDenseMatrix_ID_t selfID, int ColIndex )
{
    const Epetra_SerialDenseMatrix& self = *( 
        CEpetra::getConstSerialDenseMatrix(selfID) );

    return self[ColIndex];
}


} // extern "C"




