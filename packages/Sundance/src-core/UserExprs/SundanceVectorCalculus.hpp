/* @HEADER@ */
// ************************************************************************
// 
//                              Sundance
//                 Copyright (2005) Sandia Corporation
// 
// Copyright (year first published) Sandia Corporation.  Under the terms 
// of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government 
// retains certain rights in this software.
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
// Questions? Contact Kevin Long (krlong@sandia.gov), 
// Sandia National Laboratories, Livermore, California, USA
// 
// ************************************************************************
/* @HEADER@ */

#ifndef SUNDANCE_VECTORCALCULUS_H
#define SUNDANCE_VECTORCALCULUS_H

#include "SundanceExpr.hpp"

namespace Sundance
{


/** \relates Expr */
Expr gradient(int dim);
  
/** \relates Expr */
Expr div(const Expr& f);
  
/** \relates Expr */
Expr cross(const Expr& a, const Expr& b);
  
/** \relates Expr */
Expr curl(const Expr& f);

/** \relates Expr 
 * Compute the colon (Frobenius) product of two matrices. The result
 * is a scalar. 
 **/
Expr colonProduct(const Expr& A, const Expr& B);

/** \relates Expr 
 * Compute the outer (Kronecker) product of two vectors. The result
 * is a dim(a) by dim(b) rectangular matrix.
 **/
Expr outerProduct(const Expr& a, const Expr& b);

/** \relates Expr
 * Indicate whether the given expression is a square matrix. If so,
 * return by reference argument the size of the matrix.
 *
 * A scalar is a 1 by 1 matrix.
 */
bool isSquareMatrix(const Expr& x, int& N);

/** \relates Expr
 * Indicate whether the given expression is a vector. If so,
 * return by reference argument the size of the vector.
 *
 * A scalar is a vector of dimension 1.
 * */
bool isVector(const Expr& x, int& N);

}

#endif
