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

#ifndef SUNDANCE_GAUSS1D_H
#define SUNDANCE_GAUSS1D_H

#include "SundanceDefs.hpp"
#include "Teuchos_Array.hpp"

namespace Sundance
{
using namespace Teuchos;

/**
 * Get abscissas and weights for Gauss-Legendre quadrature 
 * on a line segement.
 */
class Gauss1D
{
public:
  /** create an n-point rule on the interval [-1, 1] */
  Gauss1D(int n);
  /** create an n-point rule on the interval [a, b] */
  Gauss1D(int n, double a, double b);

  /** return the number of points in the rule */
  int nPoints() const {return nodes_.length();}
  /** get the abscissas */
  const Array<double>& nodes() const {return nodes_;}
  /** get the weights */
  const Array<double>& weights() const {return weights_;}

  static bool unitTest();
private:
  void computeWeights(int n, double a, double b);

  Array<double> nodes_;
  Array<double> weights_;
};
}


#endif
