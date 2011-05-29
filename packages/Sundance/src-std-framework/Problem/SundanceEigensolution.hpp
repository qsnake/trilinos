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

#ifndef SUNDANCE_EIGENSOLUTION_H
#define SUNDANCE_EIGENSOLUTION_H

#include "SundanceDefs.hpp"
#include "SundanceExpr.hpp"
#include <complex>
#include "TSFVectorDecl.hpp"

namespace Sundance
{
  using namespace Sundance;
  using namespace Sundance;
  using namespace Sundance;
  using namespace Sundance;
  using namespace Sundance;
  using namespace Teuchos;

/** */
class Eigensolution
{
public:
  /** */
  Eigensolution() : eigenfunctions_(), eigenvalues_(){;}

  /** */
  Eigensolution(const Array<Expr>& eigenfunctions, 
    const Array<std::complex<double> >& eigenvalues)
    : eigenfunctions_(eigenfunctions),
      eigenvalues_(eigenvalues)
    {}

  /** */
  bool isValid() const {return eigenvalues_.size() == 0;}

  /** */
  int numEigenfunctions() const {return eigenvalues_.size();}

  /** */
  const Expr& eigenfunction(int i) const {return eigenfunctions_[i];}

  /** */
  const std::complex<double>& eigenvalue(int i) const 
    {return eigenvalues_[i];}

private:
  Array<Expr> eigenfunctions_;
  Array<std::complex<double> > eigenvalues_;
};


}


#endif
