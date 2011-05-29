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

#ifndef SUNDANCE_QUADRATURETYPEBASE_H
#define SUNDANCE_QUADRATURETYPEBASE_H

#include "SundanceDefs.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_Array.hpp"
#include "SundanceCellType.hpp"
#include "SundanceQuadratureFamily.hpp"
#include "SundanceHandleable.hpp"
#include "Teuchos_XMLObject.hpp"


namespace Sundance
{
using namespace Teuchos;

/** 
 * QuadratureTypeBase 
 */
class QuadratureTypeBase 
  : public Sundance::Handleable<QuadratureTypeBase>
{
public:
  /** */
  QuadratureTypeBase() {;}

  /** */
  virtual ~QuadratureTypeBase(){;}

  /** */
  virtual XMLObject toXML() const = 0 ;

  /** Indicate whether the given cell type is supported at any order */
  virtual bool supportsCellType(const CellType& cellType) const = 0 ;

  /** Indicate whether the given cell type is supported at the
   * specified order */
  virtual bool supports(const CellType& cellType, int order) const = 0 ;

  /** Indicate whether there is a maximum order for quadrature rules
   * available on the given cell type. */
  virtual bool hasLimitedOrder(const CellType& cellType) const = 0 ;

  /** Return the max quadrature order available on the given cell type */
  virtual int maxOrder(const CellType& cellType) const {return -1;}

  /** Create a quadrature family of the specified order */
  virtual QuadratureFamily createQuadFamily(int order) const = 0 ;
      
private:
};
}

                  

#endif

