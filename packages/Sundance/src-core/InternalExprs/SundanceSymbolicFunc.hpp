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

#ifndef SUNDANCE_SYMBOLICFUNC_H
#define SUNDANCE_SYMBOLICFUNC_H


#include "SundanceDefs.hpp"
#include "SundanceListExpr.hpp"
#include "SundanceSymbolicFuncDescriptor.hpp"
#include "SundanceFunctionWithID.hpp"
#include "SundanceDiscreteFuncElement.hpp"
#include "SundanceDiscreteFunctionStub.hpp"

namespace Sundance
{
using namespace Sundance;
using namespace Teuchos;
    

/** 
 * SymbolicFunc is a base class for functions such as test and unknown
 * functions that are "variables" in a weak form. Symbolic functions
 * cannot be evaluated directly; before evaluating a weak form,
 a value must be substituted for 
 * each symbolic func using either the substituteZero() or
 * substituteFunction() method. 
 *
 * A symbolic function will construct itself as a list of
 * SymbolicFuncElement objects that point back to the SymbolicFunction.
 */
class SymbolicFunc : public ListExpr, public SymbolicFuncDescriptor,
  public FunctionWithID
{
public:
  /** Ctor records common data but initializes list to empty */
  SymbolicFunc(const FunctionWithID& fid, 
    const RCP<const CommonFuncDataStub>& data);

  /** virtual destructor */
  virtual ~SymbolicFunc() {;}

  /** Specify that expressions involving this function are to be evaluated
   * with this function set to zero. This is appropriate for computing
   * the functional derivatives that arise in a linear problem */
  void substituteZero() const ;

  /** Specify that expressions involving this function are to be evaluated
   * with this function set to the discrete function \f$u_0\f$. 
   * This is appropriate for computing
   * the functional derivatives that arise in a nonlinear expression
   * being linearized about \f$u_0\f$. 
   */
  void substituteFunction(const RCP<DiscreteFunctionStub>& u0) const ;

  /** */
  const RCP<const CommonFuncDataStub>& commonData() const 
    {return commonData_;}
  

private:
  RCP<const CommonFuncDataStub> commonData_;

};
}


#endif
