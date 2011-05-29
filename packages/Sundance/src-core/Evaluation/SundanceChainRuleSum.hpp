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

#ifndef SUNDANCE_CHAINRULESUM_H
#define SUNDANCE_CHAINRULESUM_H

#include "SundanceDefs.hpp"
#include "SundanceEvalManager.hpp"
#include "SundanceEvaluator.hpp"
#include "Teuchos_Array.hpp"

namespace Sundance 
{

using namespace Teuchos;


/** */
class IndexPair 
{
public:
  /** */
  IndexPair(int argIndex, int valueIndex)
    : argIndex_(argIndex), valueIndex_(valueIndex) {;}
      
  /** */
  int argIndex() const {return argIndex_;}
      
  /** */
  int valueIndex() const {return valueIndex_;}
      
private:
  int argIndex_;
  int valueIndex_;
};

/** */
class DerivProduct
{
public:
  /** */
  DerivProduct() : coeff_(1.0), constants_(), variables_() {}
  /** */
  DerivProduct(const double& coeff) : coeff_(coeff), constants_(), variables_() {}

  /** */
  void addConstantFactor(const IndexPair& p) {constants_.append(p);}

  /** */
  void addVariableFactor(const IndexPair& p) {variables_.append(p);}

  /** */
  bool isConstant() const {return numVariables()==0;}

  /** */
  int numConstants() const {return constants_.size();}

  /** */
  int numVariables() const {return variables_.size();}

  /** */
  const double& coeff() const {return coeff_;}

  /** */
  const IndexPair& constant(int i) const {return constants_[i];}

  /** */
  const IndexPair& variable(int i) const {return variables_[i];}
        
private:

  double coeff_;

  Array<IndexPair> constants_;

  Array<IndexPair> variables_;
};


/** */
class ChainRuleSum : public ObjectWithClassVerbosity<Evaluator>
{
public:
  /** */
  ChainRuleSum(const MultipleDeriv& md, 
    int resultIndex,
    bool resultIsConstant);

  /** */
  void addTerm(int argDerivIndex, 
    bool argDerivIsConstant,
    const Array<DerivProduct>& sum);


  /** */
  void evalConstant(const EvalManager& mgr,
    const Array<RCP<Array<double> > >& constantArgResults,
    const Array<double>& constantArgDerivs,
    double& constResult) const ;

  /** */
  void evalVar(const EvalManager& mgr,
    const Array<RCP<Array<double> > >& constantArgResults,
    const Array<RCP<Array<RCP<EvalVector> > > > & vArgResults,
    const Array<double>& constantArgDerivs,
    const Array<RCP<EvalVector> >& varArgDerivs,
    RCP<EvalVector>& varResult) const ;

  /** */
  int resultIndex() const {return resultIndex_;}

  /** */
  bool resultIsConstant() const {return resultIsConstant_;}

  /** */
  int numTerms() const {return terms_.size();}

  /** */
  bool argDerivIsConstant(int i) const {return argDerivIsConstant_[i];}

  /** */
  int argDerivIndex(int i) const {return argDerivIndex_[i];}

  /** */
  const Array<DerivProduct>& terms(int i) const {return terms_[i];}

  /** */
  const MultipleDeriv& deriv() const {return md_;}

private:
  MultipleDeriv md_;
  int resultIndex_;
  bool resultIsConstant_;

  Array<int> argDerivIndex_;
  Array<int> argDerivIsConstant_;
  Array<Array<DerivProduct> > terms_;
};

}
               
#endif
