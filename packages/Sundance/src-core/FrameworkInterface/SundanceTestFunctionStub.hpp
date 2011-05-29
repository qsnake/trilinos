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

#ifndef SUNDANCE_TESTFUNCTIONSTUB_H
#define SUNDANCE_TESTFUNCTIONSTUB_H


#include "SundanceDefs.hpp"
#include "SundanceSymbolicFunc.hpp"
#include "SundanceTestFuncDataStub.hpp"
#include "SundanceSpectralBasis.hpp"
#include "SundanceSpectralExpr.hpp"

namespace Sundance
{
class DiscreteFunctionStub;

using namespace Teuchos;

/** 
 * TestFunctionStub is the base class for unknown functions. 
 * Each framework will need to implement its own subclass of
 * TestFunctionStub. 
 *
 * The interface is left very minimal so as to not place
 * any constraints on how a framework might specify the basis.
 * When a framework needs any information about the
 * unknown function, it will have to get it by downcasting
 * to the appropriate framework-specific subclass.
 *
 * <h4> Writing a TestFunctionStub subclass </h4>
 *
 * For purposes of interaction with the Sundance core, no 
 * additional methods are required.
 * However, most frameworks will require extensions to 
 * TestFunctionStub that can supply the framework with information
 * on the basis used by the unknown func. See the
 * demo and standard frameworks for information on how to do this.
 */
class TestFunctionStub : public SymbolicFunc
{
public:
  /** */
  TestFunctionStub(const std::string& name, 
    int tensorOrder=0, int dim=1,
    const RCP<const TestFuncDataStub>& data=RCP<const TestFuncDataStub>());

  /** */
  TestFunctionStub(const std::string& name, const SpectralBasis& sbasis, 
    int tensorOrder=0, int dim=1,
    const RCP<const TestFuncDataStub>& data=RCP<const TestFuncDataStub>());

  /** virtual destructor */
  virtual ~TestFunctionStub() {;}

  /** */
  bool isTestFunction() const {return true;}

  /** */
  virtual RCP<ExprBase> getRcp() {return rcp(this);}

protected:

  /** */
  const RCP<const TestFuncDataStub>& dataStub() const {return data_;}

private:
  RCP<const TestFuncDataStub> data_;


};

}
                  

#endif
