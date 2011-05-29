// @HEADER
// ***********************************************************************
// 
//                      Didasko Tutorial Package
//                 Copyright (2005) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
//
// Questions about Didasko? Contact Marzio Sala (marzio.sala _AT_ gmail.com)
// 
// ***********************************************************************
// @HEADER

#include "Tpetra_ConfigDefs.hpp"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Tpetra_MpiPlatform.hpp"
#include "Tpetra_MpiComm.hpp"
#else
#include "Tpetra_SerialPlatform.hpp"
#include "Tpetra_SerialComm.hpp"
#endif
#include "Tpetra_ElementSpace.hpp"
#include "Tpetra_VectorSpace.hpp"
#include "Tpetra_CisMatrix.hpp"
#include "Teuchos_ScalarTraits.hpp"

// \author Marzio Sala, ETHZ/COLAB
//
// \date Last updated on 28-Nov-05

typedef int OrdinalType;
typedef double ScalarType;

int main(int argc, char *argv[]) 
{
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Tpetra::MpiComm<OrdinalType, ScalarType> Comm(MPI_COMM_WORLD);
#else
  Tpetra::SerialComm<OrdinalType, ScalarType> Comm;
#endif

  // Get zero and one for the OrdinalType
  
  OrdinalType const OrdinalZero = Teuchos::ScalarTraits<OrdinalType>::zero();
  OrdinalType const OrdinalOne  = Teuchos::ScalarTraits<OrdinalType>::one();

  // Get zero and one for the ScalarType
  
  ScalarType const ScalarZero = Teuchos::ScalarTraits<ScalarType>::zero();
  ScalarType const ScalarOne  = Teuchos::ScalarTraits<ScalarType>::one();

  // Creates a vector of size `length', then set the elements values.
  
  OrdinalType length    = OrdinalOne * 10;
  OrdinalType indexBase = OrdinalZero;

  // 1) Creation of a platform
  
#ifdef HAVE_MPI
  const Tpetra::MpiPlatform <OrdinalType, OrdinalType> platformE(MPI_COMM_WORLD);
  const Tpetra::MpiPlatform <OrdinalType, ScalarType> platformV(MPI_COMM_WORLD);
#else
  const Tpetra::SerialPlatform <OrdinalType, OrdinalType> platformE;
  const Tpetra::SerialPlatform <OrdinalType, ScalarType> platformV;
#endif

  // 2) We can now create a space:

  Tpetra::ElementSpace<OrdinalType> elementSpace(length, indexBase, platformE);
  Tpetra::VectorSpace<OrdinalType, ScalarType> 
    vectorSpace(elementSpace, platformV);

  // 3) We now setup the matrix, which is diagonal.
  //    To that aim, we need to extract the list of locally owned
  //    ID's from elementSpace.
  
  OrdinalType NumMyElements = elementSpace.getNumMyElements();
  vector<OrdinalType> MyGlobalElements = elementSpace.getMyGlobalElements();
  
  Tpetra::CisMatrix<OrdinalType,ScalarType> matrix(vectorSpace);

  for (OrdinalType LID = OrdinalZero ; LID < NumMyElements ; ++LID)
  {
    OrdinalType GID = MyGlobalElements[LID];
    // add the diagonal element of value `GID'
    matrix.submitEntry(Tpetra::Insert, GID, (ScalarType)GID, GID);
  }

  matrix.fillComplete();

  cout << matrix;

  OrdinalType NumGlobalNonzeros   = matrix.getNumGlobalNonzeros();
  OrdinalType NumGlobalRows       = matrix.getNumGlobalRows();
  OrdinalType NumGlobalCols       = matrix.getNumGlobalCols();
  OrdinalType NumGlobalDiagonals  = matrix.getNumGlobalDiagonals();
  OrdinalType GlobalMaxNumEntries = matrix.getGlobalMaxNumEntries();
  OrdinalType isRowOriented       = matrix.isRowOriented();

  if (Comm.getMyImageID() == 0)
  {
    cout << endl << "*) Global information:" << endl << endl;
    cout << "NumGlobalNonzeros   = " << NumGlobalNonzeros << endl;
    cout << "NumGlobalRows       = " << NumGlobalRows << endl;
    cout << "NumGlobalCols       = " << NumGlobalCols << endl;
    cout << "NumGlobalDiagonals  = " << NumGlobalDiagonals << endl;
    cout << "GlobalMaxNumEntries = " << GlobalMaxNumEntries << endl;
    cout << "isRowOriented       = " << matrix.isRowOriented() << endl;
    cout << "isFillCompleted     = " << matrix.isFillCompleted() << endl;

    cout << endl;
    cout << endl << "*) Local information:" << endl << endl;
    cout << "NumMyNonzeros   = " << matrix.getNumMyNonzeros() << endl;
    cout << "NumMyRows       = " << matrix.getNumMyRows() << endl;
    cout << "NumMyCols       = " << matrix.getNumMyCols() << endl;
    cout << "NumMyDiagonals  = " << matrix.getNumMyDiagonals() << endl;
    cout << "MyMaxNumEntries = " << matrix.getMyMaxNumEntries() << endl;
  }

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

  return(EXIT_SUCCESS);
}
