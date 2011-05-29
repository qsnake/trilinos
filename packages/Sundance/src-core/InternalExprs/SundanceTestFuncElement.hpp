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

#ifndef SUNDANCE_TESTFUNCELEMENT_H
#define SUNDANCE_TESTFUNCELEMENT_H


#include "SundanceDefs.hpp"
#include "SundanceSymbolicFuncElement.hpp"
#include "SundanceTestFuncDataStub.hpp"

namespace Sundance
{
using namespace Sundance;
using namespace Teuchos;




/** 
 * TestFuncElement represents a scalar-valued element of a (possibly)
 * list-valued TestFunction
 */
class TestFuncElement : public SymbolicFuncElement
{
public:
  /** */
  TestFuncElement(const RCP<const TestFuncDataStub>& commonData,
    const std::string& name,
    const std::string& suffix, const FunctionIdentifier& fid);

  /** virtual destructor */
  virtual ~TestFuncElement() {;}

  /** Test whether all terms have test functions. 
   * I'm a test function, so return true */
  virtual bool everyTermHasTestFunctions() const {return true;}

  /** Test whether this expr contains a test function. 
   * I'm a test function, so return true. */
  virtual bool hasTestFunctions() const {return true;}

  /** */
  virtual bool isTestFunction() const {return true;}

  /** 
   * Indicate whether the expression is linear 
   * with respect to test functions */
  virtual bool isLinearInTests() const {return true;}
      
  /** Ordering operator for use in transforming exprs to standard form */
  virtual bool lessThan(const ScalarExpr* other) const ;


  /** */
  virtual XMLObject toXML() const ;

  /** */
  virtual RCP<ExprBase> getRcp() {return rcp(this);}
      
private:
};
}

#endif
