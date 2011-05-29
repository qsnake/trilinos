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


#ifndef SUNDANCE_GIDCELLPREDICATE_H
#define SUNDANCE_GIDCELLPREDICATE_H

#include "SundanceDefs.hpp"
#include "SundanceCellPredicateBase.hpp"

namespace Sundance
{
using namespace Teuchos;
  
/** 
 * GIDCellPredicate tests whether a cell's GID is contained
 * in a specified set of GIDs. 
 */
class GIDCellPredicate : public CellPredicateBase 
{
public:
  /** Construct with a GID set */
  GIDCellPredicate(int cellDim, const Set<int>& gids) 
    : CellPredicateBase(), cellDim_(cellDim), gids_(gids){;}

  /** virtual dtor */
  virtual ~GIDCellPredicate(){;}
      
  /** Test the predicate on a batch of cells */
  virtual void testBatch(const Array<int>& cellLID,
    Array<int>& results) const ;

  /** Write to XML */
  virtual XMLObject toXML() const ;

  /** comparison */
  virtual bool lessThan(const CellPredicateBase* other) const ;

  /** */
  virtual std::string description() const 
    {return "GIDCellPredicate(" + gids_.toString() + ")";}

  /* */
  GET_RCP(CellPredicateBase);

private:
  int cellDim_;
  Set<int> gids_;

};

}

#endif
