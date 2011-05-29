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

#ifndef SUNDANCE_BLOCK_H
#define SUNDANCE_BLOCK_H

#include "SundanceDefs.hpp"
#include "SundanceExpr.hpp"
#include "TSFVectorType.hpp"


namespace Sundance
{
using namespace Teuchos;
using namespace TSFExtended;
  
/** 
 * 
 */
class Block 
{
public:
  /** */
  Block()
    : expr_(), vecType_() {;}

  /** */
  Block(const Expr& expr, const VectorType<double>& vecType)
    : expr_(expr), vecType_(vecType) {;}

  /** */
  const Expr& expr() const {return expr_;}

  /** */
  const VectorType<double>& vecType() const {return vecType_;}

private:
  Expr expr_;

  VectorType<double> vecType_;
};

/** */
class BlockArray : public Array<Block>
{
public:
  /** */
  BlockArray(const Array<Block>& a) : Array<Block>(a) {;}
  /** explicit conversion needed for python wrappers */
  BlockArray(int n) : Array<Block>(n) {;}

};

/** \relates Block */
inline std::ostream& operator<<(std::ostream& os, const Block& block)
{
  os << "Block[" << block.expr() << ", " << block.vecType() << "]";
  return os;
}

}


#endif
