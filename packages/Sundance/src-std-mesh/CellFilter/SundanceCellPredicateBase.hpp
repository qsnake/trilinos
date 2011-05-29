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


#ifndef SUNDANCE_CELLPREDICATEBASE_H
#define SUNDANCE_CELLPREDICATEBASE_H


#include "SundanceDefs.hpp"
#include "SundanceMesh.hpp"
#include "SundanceHandleable.hpp"
#include "Teuchos_XMLObject.hpp"
#include "SundancePrintable.hpp"
#include "Teuchos_Describable.hpp"
#include "SundanceObjectWithVerbosity.hpp"
#include <typeinfo>

namespace Sundance
{
using namespace Teuchos;
  
/** 
 * CellPredicateBase is the base class for predicate objects
 * that test cells
 * against some condition. A simulation developer needing
 * some specialized method for identifying cells might implement
 * a custom cell predicate by extending this function. However,
 * the most common cases, selection by cell label or cell position,
 * have already been implemented 
 * in LabelCellPredicate and PositionalCellPredicate.
 */
class CellPredicateBase 
  : public Sundance::Handleable<CellPredicateBase>,
    public Noncopyable,
    public ObjectWithClassVerbosity<CellPredicateBase>
{
public:
  /** Empty ctor */
  CellPredicateBase();

  /** virtual dtor */
  virtual ~CellPredicateBase();
      
  /** Test the predicate on a batch of cells */
  virtual void testBatch(const Array<int>& cellLID,
    Array<int>& results) const = 0 ;

      
  /** Set the current mesh and dimension on which cells are to be tested */
  virtual void setMesh(const Mesh& mesh, int cellDim) const 
    {mesh_ = mesh; cellDim_ = cellDim;}

  /** Write to XML */
  virtual XMLObject toXML() const = 0 ;

  /** */
  virtual bool lessThan(const CellPredicateBase* other) const = 0 ;

  /** */
  virtual std::string description() const = 0 ;


  /** */
  virtual std::string typeName() const {return typeid(*this).name();}
protected:

  /** */
  const Mesh& mesh() const {return mesh_;}


  /** */
  int cellDim() const {return cellDim_;}

private:
  mutable Mesh mesh_;

  mutable int cellDim_;



};

}

#endif
