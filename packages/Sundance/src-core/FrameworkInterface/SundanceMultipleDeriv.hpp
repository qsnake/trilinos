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

#ifndef SUNDANCE_MULTIPLEDERIV_H
#define SUNDANCE_MULTIPLEDERIV_H

#include "SundanceDefs.hpp"
#include "SundanceDeriv.hpp"
#include "SundanceMultiSet.hpp"
#include "SundanceOrderedTuple.hpp"
#include "SundanceMap.hpp"
#include "Teuchos_Array.hpp"



namespace Sundance
{
using namespace Sundance;
class MultipleDeriv;
/** */
typedef OrderedPair<MultipleDeriv, MultipleDeriv> DerivPair; 

/** */
typedef Sundance::Map<DerivPair, int> ProductRulePerms;
/** Class MultipleDeriv is a multiple functional derivative operator
 * represented as a multiset of first-order derivatives.
 * The derivatives are each Deriv objects, so the multiple
 * derivative can be an arbitrary mixture of spatial and
 * functional derivatives.
 *
 * <h3> Product rule application </h3>
 *
 * Class MultipleDeriv contains a utility method for computing the
 * derivatives that appear in the general product rule, as
 * used in automatic differentiation of products and spatial
 * derivatives.
 * The arbitrary order, multivariable product rule gives
 * the derivative of a product \f$L\cdot R\f$
 * with respect to a multiset of variables
 * \f$\{u_1, u_2, \dots, u_N\}\f$ as the sum
 * \f[ D_{u_1} D_{u_2} \cdots D_{u_N} (L \cdot R)
 * =\sum_{i=1}^{2^N} \left[\prod_{j=1}^N D_{u_j}^{b_{i,j}}\right] L
 * \times \left[\prod_{j=1}^N D_{u_j}^{1-b_{i,j}}\right] R\f]
 * where \f$b_{i,j}\f$ is the \f$j\f$-th bit in the binary
 * representation of \f$i\f$. Method productRulePermutations()
 * enumerates all the permutations that apply to the left and right
 * operands in this sum, and returns the
 * arrays of left operators and right operators
 * through non-const reference arguments.
 * The left and right arguments return
 *
 * \f[ {\rm left = Array}_{i=1}^{2^N}
 * \{\prod_{j=1}^N D_{u_j}^{b_{i,j}}\}\f]
 *
 * \f[ {\rm right = Array}_{i=1}^{2^N}
 * \{\prod_{j=1}^N D_{u_j}^{1-b_{i,j}}\}\f]
 *
 * Each element in these arrays is a multiple derivative, and
 * is naturally represented with a MultipleDeriv object.
 */
class MultipleDeriv : public MultiSet<Deriv>
{
public:

  /** Construct an empty multiple derivative */
  MultipleDeriv();

  /** Construct a first-order multiple deriv */
  MultipleDeriv(const Deriv& d);

  /** Construct a second-order multiple deriv */
  MultipleDeriv(const Deriv& d1, const Deriv& d2);

  /** Return the order of this derivative. Since a
   * multiple derivative is a multiset of first-order
   * derivatives, the order of the multiple derivative is the
   * size of the set. */
  int order() const {return this->size();}

  /** Return the order of spatial differentiation in this
   * derivative */
  int spatialOrder() const ;

  /** Return a multiindex representing the spatial derivs in this
   * multiple derivative */
  MultiIndex spatialDeriv() const ;
          
          
  /** Return by reference argument the partitioning of
   * derivatives when this derivative is applied to a product.
   */
  void productRulePermutations(ProductRulePerms& perms) const ;

  /** return the product of two multiple derivatives, which
   * is the union of the two multisets */
  MultipleDeriv product(const MultipleDeriv& other) const ;


  /** 
   * Return a copy of this derivative, but with the specified single
   * derivative removed. For example, if <t>this</t> is
   * \f$ \{u,v,w\} \f$, calling <t>factorOutDeriv(v)</t> will
   * return \f$\{u,w\}\f$.
   */
  MultipleDeriv factorOutDeriv(const Deriv& x) const ;


  /** 
   * Return a copy of this derivative, but with all single derivs
   * in the the argument removed. For example, if <t>this</t> is
   * \f$ \{u,v,w\} \f$, calling <t>factorOutDeriv(u,v)</t> will
   * return \f$\{w\}\f$.
   */
  MultipleDeriv factorOutDeriv(const MultipleDeriv& x) const ;


  /** 
   * Return true is all single derivatives appearing in x also appear in this.
   */
  bool containsDeriv(const MultipleDeriv& x) const ;


  /** 
   *
   */
  bool isInRequiredSet(const Set<MultiSet<int> >& funcCombinations,
    const Set<MultiIndex>& multiIndices) const ;

  /** */
  MultiSet<FunctionIdentifier> funcIDs() const ;

  /** */
  MultiSet<int> dofIDs() const ;

  /** \name Utilities used in computing product rule permutations */
  //@
  /** Compute the n-th power of 2 */
  static int pow2(int n);

  /** Compute the n bits of an integer x < 2^n. The bits are ordered
   * least-to-most significant.
   */
  static Array<int> bitsOfAnInteger(int x, int n);
  //@}
private:
};


/** 
 * \brief Tranform the input set S as follows: for each multiple derivative
 * in S, apply the multiindex x to each of its component functional derivatives 
 * in turn. Note that the multiindex may have negative indices. Results with
 * negative multiindices are ignored.
 */
Set<MultipleDeriv> applyTx(const Set<MultipleDeriv>& S,
  const MultiIndex& x);
/** 
 * \brief Filter the input set W, allowing only coordinate derivatives in the 
 * direction x and functional derivatives whose associated functions
 * have nonzero evaluation points. This function is used 
 * in the spatial/functional chain rule to identify those terms resulting
 * from differentiation wrt x or functional derivatives. 
 */
Set<MultipleDeriv> applyZx(const Set<MultipleDeriv>& W,
  const MultiIndex& x);


/** */
Set<MultipleDeriv> Xx(const MultiIndex& x) ;

/** */
int factorial(const MultipleDeriv& ms);


/**
 *
 */
template <class T> class increasingOrder
: std::binary_function<T, T, bool>
{
public:
  bool operator()(const MultipleDeriv& a,
    const MultipleDeriv& b) const;
};
/**
 * Functor increasingOrder() is used as a comparison method
 * for MultipleDeriv objects. When used in a set, it will sort
 * derivatives in order of increasing order of differentiation.
 * Sorting by differentiation order is useful in evaluation of
 * product and chain rules: if we evaluate higher-order derivatives
 * first, we can evaluate in place without destroying lower-order
 * derivatives.
 */
template <> class increasingOrder<MultipleDeriv>
{
public:
  bool operator()(const MultipleDeriv& a,
    const MultipleDeriv& b) const
    {
      if (a.order() < b.order()) return true;
      if (a.order() > b.order()) return false;

      MultipleDeriv::const_iterator aIter;
      MultipleDeriv::const_iterator bIter;

      for (aIter=a.begin(), bIter=b.begin();
           aIter != a.end() && bIter != b.end();
           aIter++, bIter++)
      {
        if (*aIter < *bIter) return true;
        if (*bIter < *aIter) return false;
      }

      return false;
    }
};


/** */
MultipleDeriv makeMultiDeriv(const Deriv& d);

/** */
bool hasParameter(const MultipleDeriv& d) ;


}


#endif
