// @HEADER
// ************************************************************************
//
//           Galeri: Finite Element and Matrix Generation Package
//                 Copyright (2006) ETHZ/Sandia Corporation
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
// Questions about Galeri? Contact Marzio Sala (marzio.sala _AT_ gmail.com)
//
// ************************************************************************
// @HEADER

#ifndef GALERI_STRETCHED2DMATRIX_H
#define GALERI_STRETCHED2DMATRIX_H

#include "Galeri_Exception.h"
#include "Galeri_Utils.h"
#include "Epetra_Comm.h"
#include "Epetra_BlockMap.h"
#include "Epetra_CrsMatrix.h"

namespace Galeri {
namespace Matrices {

  //def for espilon muyst be 1e-5
inline
Epetra_CrsMatrix* 
Stretched2D(const Epetra_Map* Map, const int nx, const int ny,
            const double epsilon)
{
  Epetra_CrsMatrix* Matrix = new Epetra_CrsMatrix(Copy, *Map,  9);

  int NumMyElements = Map->NumMyElements();
  int* MyGlobalElements = Map->MyGlobalElements();

  int left, right, lower, upper;
    
  double Values[9];
  int    Indices[9];

  for (int i = 0 ; i < NumMyElements ; ++i) 
  {
    int NumEntries = 0;
    GetNeighboursCartesian2d(MyGlobalElements[i], nx, ny, 
			     left, right, lower, upper);

    if (left != -1) 
    {
      Indices[NumEntries] = left;
      Values[NumEntries] = 2.0 - epsilon;
      ++NumEntries;
    }
    if (right != -1) 
    {
      Indices[NumEntries] = right;
      Values[NumEntries] = 2.0 - epsilon;
      ++NumEntries;
    }
    if (lower != -1) 
    {
      Indices[NumEntries] = lower;
      Values[NumEntries] = -4.0 + epsilon;
      ++NumEntries;
    }
    if (upper != -1) 
    {
      Indices[NumEntries] = upper;
      Values[NumEntries] = -4.0 + epsilon;
      ++NumEntries;
    }
    if (left != -1 && lower != -1) 
    {
      Indices[NumEntries] = lower - 1;
      Values[NumEntries] = -1.0;
      ++NumEntries;
    }
    if (right != -1 && lower != -1) 
    {
      Indices[NumEntries] = lower + 1;
      Values[NumEntries] = -1.0;
      ++NumEntries;
    }
    if (left != -1 && upper != -1) 
    {
      Indices[NumEntries] = upper - 1;
      Values[NumEntries] = -1.0;
      ++NumEntries;
    }
    if (right != -1 && upper != -1) 
    {
      Indices[NumEntries] = upper + 1;
      Values[NumEntries] = -1.0;
      ++NumEntries;
    }
    
    Indices[NumEntries] = MyGlobalElements[i];
    Values[NumEntries] = 8.0;
    ++NumEntries;

    Matrix->InsertGlobalValues(MyGlobalElements[i], NumEntries, 
                               Values, Indices);
  }

  Matrix->FillComplete();
  Matrix->OptimizeStorage();

  return(Matrix);
}

} // namespace Matrices
} // namespace Galeri
#endif
