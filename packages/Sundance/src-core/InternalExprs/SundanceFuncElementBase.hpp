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

#ifndef SUNDANCE_FUNCELEMENTBASE_H
#define SUNDANCE_FUNCELEMENTBASE_H


#include "SundanceDefs.hpp"
#include "SundanceEvaluatableExpr.hpp"
#include "SundanceEvalContext.hpp"
#include "SundanceMultiSet.hpp"
#include "SundanceSet.hpp"
#include "SundanceMultipleDeriv.hpp"
#include "SundanceFunctionWithID.hpp"


namespace Sundance
{
using namespace Sundance;
using namespace Teuchos;





/** 
 * FuncElementBase defines the interface for scalar-valued elements
 * of Sundance functions. At the user level, Sundance functions can be
 * list (e.g, vector or tensor) valued; internally, however, compound
 * expressions use only scalar functions deriving from the 
 * FuncElementBase class. 
 */
class FuncElementBase : public virtual ScalarExpr,
  public FunctionWithID
{
public:
  /** */
  FuncElementBase(const std::string& rootName,
    const std::string& suffix,
    const FunctionIdentifier& fid);
  /** */
  FuncElementBase(const std::string& rootName);

  /** virtual destructor */
  virtual ~FuncElementBase() {;}

  /** Return the name of this function */
  const std::string& name() const {return name_;}

  /** Return the root name of this function */
  const std::string& rootName() const {return rootName_;}

  /** Return the root name of this function */
  const std::string& suffix() const {return suffix_;}

  /** Write self in text form */
  virtual std::ostream& toText(std::ostream& os, bool paren) const ;

  /** Ordering operator for use in transforming exprs to standard form */
  virtual bool lessThan(const ScalarExpr* other) const ;


protected:
private:

  std::string name_;

  std::string rootName_;

  std::string suffix_;
};

}

#endif
