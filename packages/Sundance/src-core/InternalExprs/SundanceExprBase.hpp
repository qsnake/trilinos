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

#ifndef SUNDANCE_EXPRBASE_H
#define SUNDANCE_EXPRBASE_H


#include "SundanceDefs.hpp"
#include "SundanceExceptions.hpp"
#include "SundanceSet.hpp"
#include "Teuchos_XMLObject.hpp"
#include "Teuchos_RefCountPtrDecl.hpp"
#include "SundanceHandleable.hpp"



namespace Sundance
{
using namespace Sundance;
using namespace Teuchos;





/** */
class ExprBase : public Sundance::Handleable<ExprBase>
{
public:
  /** empty ctor */
  ExprBase();

  /** virtual destructor */
  virtual ~ExprBase() {;}

  /** Write a simple text description suitable 
   * for output to a terminal */
  virtual std::ostream& toText(std::ostream& os, bool paren) const = 0 ;

  /** Append to the set of func IDs present in this expression.
   * Base class does nothing */
  virtual void accumulateFuncSet(Set<int>& funcIDs, 
    const Set<int>& activeSet) const {;}

  /** Indicate whether this expression contains any test 
   * functions. Default is to return false. This will be
   * overridden by TestFuncElement and ExprWithChildren. */
  virtual bool hasTestFunctions() const {return false;}
  /** 
   * Indicate whether the expression contains unknown functions */
  virtual bool hasUnkFunctions() const {return false;}

  /** */
  std::string toString() const ;

  /** Write in XML */
  virtual XMLObject toXML() const = 0 ;

  /** Return a descriptive name for the expression subtype */
  virtual std::string typeName() const ;

protected:
};



}

#endif
