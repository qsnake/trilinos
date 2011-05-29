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

#ifndef SUNDANCE_SPARSITYSUPERSET_H
#define SUNDANCE_SPARSITYSUPERSET_H



#include "SundanceDefs.hpp"
#include "SundanceDefs.hpp"
#include "SundanceDerivSet.hpp"
#include "SundanceDeriv.hpp"
#include "SundanceMultipleDeriv.hpp"
#include "SundanceMap.hpp"
#include "SundanceEvalVector.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_Array.hpp"


namespace Sundance
{
using namespace Sundance;
using namespace Teuchos;

/**
 * DerivState can be used to classify the known state of
 * each functional derivative at any node in the expression
 * tree.  ZeroDeriv means the derivative is structurally
 * zero at that node -- note that derivatives that are
 * nonzero at a root node can be structurally zero at some
 * child node. ConstantDeriv means that the derivative is
 * known to be a constant (in space) at that
 * node. VectorDeriv means that the derivative is non-zero,
 * non-constant, i.e., a vector of values.
 */
enum DerivState {ZeroDeriv, ConstantDeriv, VectorDeriv};

/**
 *
 */
class SparsitySuperset 
  : public ObjectWithClassVerbosity<SparsitySuperset>
{
public:
          
  /** Create a sparsity set */
  SparsitySuperset(const Set<MultipleDeriv>& C,
    const Set<MultipleDeriv>& V);

  /** \name Access to information about individual derivatives */
  //@{
  /** Detect whether a given derivative exists in this set */
  bool containsDeriv(const MultipleDeriv& d) const ;

  /** Find the index at which the results for the
      given functional derivative are stored in the results array */
  int getIndex(const MultipleDeriv& d) const ;

  /** Return the results stored at index i */
  inline const MultipleDeriv& deriv(int i) const
    {return derivs_[i];}

  /** Return the constancy state of deriv i */
  inline const DerivState& state(int i) const
    {return states_[i];}

  /** */
  inline int numDerivs() const {return derivs_.size();}

  /** */
  int numConstantDerivs() const {return numConstantDerivs_;}

  /** */
  int numVectorDerivs() const {return numVectorDerivs_;}

          
          
  /** */
  int maxOrder() const {return maxOrder_;}


  /** Indicate whether the specified derivative is
   * spatially constant at this node */
  bool isConstant(int i) const {return states_[i]==ConstantDeriv;}

  /** Indicate whether the specified multiple derivative contains
   * at least one order of spatial differentiation */
  bool isSpatialDeriv(int i) const 
    {
      return multiIndex_[i].order() != 0;
    }

  /** Return the spatial multi index for the i-th derivative */
  const MultiIndex& multiIndex(int i) const 
    {return multiIndex_[i];}
  //@}
          
  /** */
  void print(std::ostream& os) const ;
  /** */
  void print(std::ostream& os, 
    const Array<RCP<EvalVector> >& vecResults,
    const Array<double>& constantResults) const ;

  /** */
  std::string toString() const ;

  /** */
  DerivSet derivSet() const ;

private:

  /** Add a derivative to the set. Called by the subset's addDeriv()
   * method. */
  void addDeriv(const MultipleDeriv& d, 
    const DerivState& state);

  /** Add a derivative to the set. Called by the subset's addDeriv()
   * method. */
  void addDeriv(const Deriv& d, 
    const DerivState& state);

  /** */
  int maxOrder_;

  /** Map from deriv to position of the derivative's
   * value in the results array */
  Sundance::Map<MultipleDeriv, int> derivToIndexMap_;

  /** The list of functional derivatives whose values are
   * stored in this results set */
  Array<MultipleDeriv> derivs_;

  /** The state of each derivative at this node in the expression */
  Array<DerivState> states_;

  /** Multiindices */
  Array<MultiIndex> multiIndex_;

  /** */
  int numConstantDerivs_;

  /** */
  int numVectorDerivs_;

          

};

}

namespace std
{
/** \relates SparsitySuperset */
inline std::ostream& operator<<(std::ostream& os,
  const Sundance::SparsitySuperset& s)
{
  s.print(os);
  return os;
}
}



#endif
