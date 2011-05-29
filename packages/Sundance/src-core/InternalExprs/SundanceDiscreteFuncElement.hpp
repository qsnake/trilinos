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

#ifndef SUNDANCE_DISCRETEFUNCELEMENT_H
#define SUNDANCE_DISCRETEFUNCELEMENT_H


#include "SundanceDefs.hpp"
#include "SundanceFuncElementBase.hpp"
#include "SundanceDiscreteFuncEvaluator.hpp"
#include "SundanceDiscreteFuncDataStub.hpp"

namespace Sundance
{
using namespace Sundance;


using namespace Teuchos;




/** 
 * DiscreteFuncElement represents a scalar-valued element
 * of a (possibly) vector-valued discrete function. 
 *
 * DiscreteFuncElement is framework-independent. Any framework-specific
 * information should go in a subclass of DiscreteFuncDataStub.
 * The DiscreteFuncDataStub object can be accessed through the
 * <tt>master()</tt> method of this class.
 */
class DiscreteFuncElement : public virtual EvaluatableExpr,
                            public FuncElementBase,
                            public virtual GenericEvaluatorFactory<DiscreteFuncElement, DiscreteFuncElementEvaluator>
{
public:
  /** */
  DiscreteFuncElement(const RCP<DiscreteFuncDataStub>& data,
    const std::string& name,
    const std::string& suffix,
    const FunctionIdentifier& fid,
    int myIndexIntoVector);

  /** virtual destructor */
  virtual ~DiscreteFuncElement() {;}


  /** Get the data associated with the vector-valued function 
   * that contains this function element. */
  RCP<const DiscreteFuncDataStub> commonData() const {return commonData_;}

  /** Get the data associated with the vector-valued function 
   * that contains this function element. */
  DiscreteFuncDataStub* commonData() {return commonData_.get();}

  /** Get my index into the master's list of elements */
  int myIndex() const {return myIndex_;}

  /** Inform this function that it will need to be evaluated using the specified
   * multiIndex*/
  void addMultiIndex(const MultiIndex& newMi) const ;

  /**
   * Find the maximum differentiation order acting on discrete
   * functions in this expression. 
   */
  int maxDiffOrderOnDiscreteFunctions() const {return 0;}
      
  /**
   * Indicate whether this expression contains discrete functions.
   * This object is a discrete function, so return true.
   */
  virtual bool hasDiscreteFunctions() const {return true;}
      
  /**
   * Indicate whether this expression contains test functions.
   * This object is a discrete function, so return false.
   */
  virtual bool hasTestFunctions() const {return false;}

  /** */
  virtual Set<MultipleDeriv> 
  internalFindW(int order, const EvalContext& context) const ;
  /** */
  virtual Set<MultipleDeriv> 
  internalFindV(int order, const EvalContext& context) const ;
  /** */
  virtual Set<MultipleDeriv> 
  internalFindC(int order, const EvalContext& context) const ;

  /** */
  virtual RCP<Array<Set<MultipleDeriv> > > 
  internalDetermineR(const EvalContext& context,
    const Array<Set<MultipleDeriv> >& RInput) const ;

  /** */
  virtual XMLObject toXML() const ;

  /** */
  const Set<MultiIndex>& multiIndexSet() const {return miSet_;}

  /** */
  virtual RCP<ExprBase> getRcp() {return rcp(this);}

  /** */
  bool lessThan(const ScalarExpr* other) const ;
      
private:

  RCP<DiscreteFuncDataStub> commonData_;

  mutable Set<MultiIndex> miSet_;

  int myIndex_;
      

};
}

#endif
