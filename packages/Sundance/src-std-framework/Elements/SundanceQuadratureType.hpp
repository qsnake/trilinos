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

#ifndef SUNDANCE_QUADRATURETYPE_H
#define SUNDANCE_QUADRATURETYPE_H

#include "SundanceDefs.hpp"
#include "SundanceQuadratureTypeBase.hpp"
#include "SundanceHandle.hpp"

namespace Sundance
{

/** 
 * QuadratureFamily is a geometry-independent specification of
 * a method by which quadrature is to be carried out. For example,
 * a GaussianQuadrature family will generate Gaussian
 * quadrature points on any cell type.
 */
class QuadratureType : public Sundance::Handle<QuadratureTypeBase>
{
public:
  /* */
  HANDLE_CTORS(QuadratureType, QuadratureTypeBase);

  /** */
  XMLObject toXML() const {return ptr()->toXML();}

  /** Indicate whether the given cell type is supported at any order */
  bool supportsCellType(const CellType& cellType) const 
    {return ptr()->supportsCellType(cellType);}
    
  /** Indicate whether the given cell type is supported at the
   * specified order */
  bool supports(const CellType& cellType, int order) const 
    {return ptr()->supports(cellType, order);}
    
  /** Return the max quadrature order available on the given cell type */
  int maxOrder(const CellType& cellType) const 
    {return ptr()->maxOrder(cellType);}

  /** Indicate whether there is a maximum order for quadrature rules
   * available on the given cell type. */
  bool hasLimitedOrder(const CellType& cellType) const 
    {return ptr()->hasLimitedOrder(cellType);}

  /** Create a quadrature family of the specified order */
  QuadratureFamily createQuadFamily(int order) const 
    {return ptr()->createQuadFamily(order);}

  /** */
  int findValidOrder(const CellType& cellType, int minReqOrder) const ;
private:
};
}

#endif
