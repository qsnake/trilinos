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

#ifndef SUNDANCE_FIATQUADRATURE_H
#define SUNDANCE_FIATQUADRATURE_H

#include "SundanceDefs.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "SundanceQuadratureFamilyBase.hpp"
#include "Python.h"

namespace Sundance
{
  using namespace Teuchos;

  /** 
   * Family of optimal Gaussian integration rules, e.g., Gauss-Legendre on 
   * lines, Dunavant on triangles. 
   */
  class FIATQuadratureAdapter : public QuadratureFamilyBase
  {
  public:
    /** */
    FIATQuadratureAdapter( PyObject *py_quad_factory , int order);

    /** */
    virtual ~FIATQuadratureAdapter(){;}

    void getPoints( const CellType & cellType ,
			    Array<Point>& quadPoints,
			    Array<double>& quadWeights ) const;

    /** */
    virtual XMLObject toXML() const ;

    /** Describable interface */
    virtual std::string description() const 
    {return "FIATQuadrature[order=" + Teuchos::toString(order()) 
         +  "]";}

    /* handleable boilerplate */
    GET_RCP(QuadratureFamilyStub);

  protected:
    Array<Array<Point> > pts_;
    Array<Array<double> > wts_;

  };
}


#endif
