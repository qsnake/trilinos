//@HEADER
// ***********************************************************************
// 
//                      Didasko Tutorial Package
//                 Copyright (2009) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
//@HEADER

#include "Ifpack.h"
#include "Ifpack_Euclid.h"
#include "AztecOO.h"
#include "Galeri_Maps.h"
#include "Galeri_CrsMatrices.h"
#include "Galeri_Utils.h"
#include "Epetra_MultiVector.h"

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_RCP.hpp"
#include "Epetra_ConfigDefs.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "euclid_Helpers.hpp"
#include "Teuchos_Array.hpp"
#include <string>
#include <stdio.h>
#include <map>

using Teuchos::RCP;
using Teuchos::rcp;

const double tol = 1E-7;
const int numVec = 1;

TEUCHOS_UNIT_TEST( Ifpack_Hypre, Euclid){
  RCP<Epetra_CrsMatrix> Matrix = rcp(newCrsMatrix(27));
  //cout << endl << *Matrix << endl;
  TEST_EQUALITY(Matrix->RowMap().LinearMap(), true);
  Ifpack_Euclid preconditioner(Matrix.get());
  TEST_EQUALITY(preconditioner.Initialize(),0);
  TEST_EQUALITY(preconditioner.SetParameter("setmem", 0),0);
  TEST_EQUALITY(preconditioner.SetParameter("setstats", 0), 0);
  TEST_EQUALITY(preconditioner.SetParameter("setlevel", 2), 0);
  TEST_EQUALITY(preconditioner.Compute(),0);
  cout << endl << preconditioner << endl;
  Epetra_MultiVector KnownX(Matrix->DomainMap(), numVec);
  KnownX.Random();

  Epetra_MultiVector B(Matrix->RangeMap(), numVec);
  TEST_EQUALITY(Matrix->Apply(KnownX, B), 0);

  Epetra_MultiVector X(Matrix->DomainMap(), numVec);

  AztecOO Solver(Matrix.get(), &X, &B);
  TEST_EQUALITY(Solver.SetPrecOperator(&preconditioner),0);
  Solver.Iterate(1000, 1E-9);
  //TEST_EQUALITY(preconditioner.ApplyInverse(B,X),0);
  TEST_EQUALITY(EquivalentVectors(X, KnownX, tol), true);
  cout << endl << preconditioner << endl;
}

