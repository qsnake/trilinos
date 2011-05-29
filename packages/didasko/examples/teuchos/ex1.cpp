
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

#include "Didasko_ConfigDefs.h"
#if defined(HAVE_DIDASKO_TEUCHOS)
#ifdef HAVE_MPI
#include "mpi.h"
#endif

#include "Teuchos_SerialDenseMatrix.hpp"

int main(int argc, char* argv[])
{
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  // Creating a double-precision matrix can be done in several ways:
  // Create an empty matrix with no dimension
  Teuchos::SerialDenseMatrix<int,double> Empty_Matrix;
  // Create an empty 3x4 matrix
  Teuchos::SerialDenseMatrix<int,double> My_Matrix( 3, 4 );
  // Basic copy of My_Matrix
  Teuchos::SerialDenseMatrix<int,double> My_Copy1( My_Matrix ),
    // (Deep) Copy of principle 3x3 submatrix of My_Matrix
    My_Copy2( Teuchos::Copy, My_Matrix, 3, 3 ),
    // (Shallow) Copy of 2x3 submatrix of My_Matrix
    My_Copy3( Teuchos::View, My_Matrix, 2, 3, 1, 1 );

  // The matrix dimensions and strided storage information can be obtained:
  int rows, cols, stride;
  rows = My_Copy3.numRows();  // number of rows
  cols = My_Copy3.numCols();  // number of columns
  stride = My_Copy3.stride(); // storage stride

  // Matrices can change dimension:
  Empty_Matrix.shape( 3, 3 );      // size non-dimensional matrices
  My_Matrix.reshape( 3, 3 );       // resize matrices and save values

  // Filling matrices with numbers can be done in several ways:
  My_Matrix.random();             // random numbers
  My_Copy1.putScalar( 1.0 );      // every entry is 1.0
  My_Copy2(1,1) = 10.0;           // individual element access
  Empty_Matrix = My_Matrix;       // copy My_Matrix to Empty_Matrix 

  // Basic matrix arithmetic can be performed:
  Teuchos::SerialDenseMatrix<int,double> My_Prod( 3, 2 );
  // Matrix multiplication ( My_Prod = 1.0*My_Matrix*My_Copy^T )
  My_Prod.multiply( Teuchos::NO_TRANS, Teuchos::TRANS, 
		    1.0, My_Matrix, My_Copy3, 0.0 );
  My_Copy2 += My_Matrix;         // Matrix addition
  My_Copy2.scale( 0.5 );         // Matrix scaling
  
  // The pointer to the array of matrix values can be obtained:
  double *My_Array, *My_Column;
  My_Array = My_Matrix.values();   // pointer to matrix values
  My_Column = My_Matrix[2];        // pointer to third column values

  // The norm of a matrix can be computed:
  double norm_one, norm_inf, norm_fro;
  norm_one = My_Matrix.normOne();        // one norm
  norm_inf = My_Matrix.normInf();        // infinity norm
  norm_fro = My_Matrix.normFrobenius();  // frobenius norm

  // Matrices can be compared:
  // Check if the matrices are equal in dimension and values
  if (Empty_Matrix == My_Matrix) {
    cout<< "The matrices are the same!" <<endl;
  }
  // Check if the matrices are different in dimension or values
  if (My_Copy2 != My_Matrix) {
    cout<< "The matrices are different!" <<endl;
  }

  // A matrix can be sent to the output stream:
  cout<< My_Matrix << endl;

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  return 0;
}

#else

#include <stdlib.h>
#include <stdio.h>
#ifdef HAVE_MPI
#include "mpi.h"
#endif

int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  puts("Please configure Didasko with:\n"
       "--enable-teuchos");

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  return 0;
}
#endif
