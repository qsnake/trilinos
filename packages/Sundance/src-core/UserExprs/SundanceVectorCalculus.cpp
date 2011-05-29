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

#include "SundanceVectorCalculus.hpp"
#include "SundanceDerivative.hpp"
#include "SundanceListExpr.hpp"


namespace Sundance
{

Expr gradient(int dim)
{
  Array<Expr> rtn(dim);
  for (int i=0; i<dim; i++)
  {
    rtn[i] = new Derivative(i);
  }
  return new ListExpr(rtn);
}

Expr div(const Expr& f)
{
  Expr rtn = 0.0;
  for (int i=0; i<f.size(); i++)
  {
    Expr d = new Derivative(i);
    rtn = rtn + d*f[i];
  }
  return rtn;
}


Expr cross(const Expr& a, const Expr& b)
{
  TEST_FOR_EXCEPTION(a.size() != b.size(), RuntimeError,
    "mismatched vector sizes in cross(a,b): a.size()=" << a.size()
    << ", b.size()=" << b.size());

  TEST_FOR_EXCEPTION(a.size() < 2 || a.size() > 3, RuntimeError,
    "cross(a,b) undefined for dim=" << a.size());

  if (a.size()==2)
  {
    return a[0]*b[1] - a[1]*b[0];
  }
  else
  {
    return List(
      cross(List(a[1],a[2]), List(b[1],b[2])),
      -cross(List(a[0],a[2]), List(b[0],b[2])),
      cross(List(a[0],a[1]), List(b[0],b[1]))
      );
  }
}

Expr curl(const Expr& f)
{
  Expr del = gradient(f.size());

  return cross(del, f);
}

Expr colonProduct(const Expr& A, const Expr& B)
{
  int nA = 0;
  int nB = 0;
  TEST_FOR_EXCEPTION(!isSquareMatrix(A, nA), RuntimeError,
    "Colon product expected argument A=" << A << " to be a square matrix");
  TEST_FOR_EXCEPTION(!isSquareMatrix(B, nB), RuntimeError,
    "Colon product expected argument B=" << B << " to be a square matrix");

  TEST_FOR_EXCEPTION(nA!=nB, RuntimeError,
    "Colon product expected operands A=" << A << " and B=" << B 
    << " to have identical sizes");

  Expr rtn = 0.0;
  for (int i=0; i<nA; i++)
  {
    for (int j=0; j<nA; j++)
    {
      rtn = rtn + A[i][j] * B[i][j];
    }
  }

  return rtn;
}

Expr outerProduct(const Expr& A, const Expr& B)
{
  int nA = 0;
  int nB = 0;
  TEST_FOR_EXCEPTION(!isVector(A, nA), RuntimeError,
    "Outer product expected argument A=" << A << " to be a vector");
  TEST_FOR_EXCEPTION(!isVector(B, nB), RuntimeError,
    "Outer product expected argument B=" << B << " to be a vector");

  TEST_FOR_EXCEPTION(nA!=nB, RuntimeError,
    "Colon product expected operands A=" << A << " and B=" << B 
    << " to have identical sizes");

  Array<Expr> rtn(nA);
  for (int i=0; i<nA; i++)
  {
    rtn[i] = A[i] * B;
  }

  return new ListExpr(rtn);
}

bool isVector(const Expr& x, int& N)
{
  N = 0;
  if (x.size() == x.totalSize()) 
  {
    N = x.size();
    return true;
  }
  else
  {
    return false;
  }
}


bool isSquareMatrix(const Expr& x, int& N)
{
  N = 0;
  /* do the simplest checks first */
  if (x.size()*x.size() != x.totalSize()) return false;
  N = x.size();
  for (int i=0; i<N; i++)
  {
    int M = 0;
    if (!isVector(x[i],M)) return false;
    if (M != N) return false;
  }
  return true;
}




}

