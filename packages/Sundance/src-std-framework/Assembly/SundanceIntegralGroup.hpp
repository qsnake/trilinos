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

#ifndef SUNDANCE_INTEGRALGROUP_H
#define SUNDANCE_INTEGRALGROUP_H

#include "SundanceDefs.hpp"
#include "SundanceElementIntegral.hpp"
#include "SundanceEvalVector.hpp"
#include "SundanceMultipleDeriv.hpp"
#include "SundanceIntegrationCellSpecifier.hpp"



namespace Sundance
{
using namespace Teuchos;
/** 
 * IntegralGroup collects together several integrals having common properties that
 * can be used either to eliminate redundant computations or to share storage
 * an eliminate repeated allocations. 
 *
 * The integrations are done using calls to methods of the ElementIntegral class. Which
 * subtype of ElementIntegral is used (e.g., RefIntegral, QuadratureIntegral) controls
 * how the integrations are actually done. Reference integration does the integral exactly
 * on a reference element, then transforms to a physical cell. This is only possible
 * with constant coefficients on an affine cell. Quadrature integration proceeds
 * by quadrature (duh!), and is needed for variable coefficients or non-affine cells. 
 */
class IntegralGroup 
{
public:
  /** */
  IntegralGroup(const Array<RCP<ElementIntegral> >& integrals,
    const Array<int>& resultIndices,
    int verb);

  /** */
  IntegralGroup(const Array<int>& testID,
    const Array<int>& testBlock,
    const Array<int>& mvIndices,
    const Array<RCP<ElementIntegral> >& integrals,
    const Array<int>& resultIndices,
    const Array<MultipleDeriv>& derivs,
    int verb);

  /** */
  IntegralGroup(const Array<int>& testID,
    const Array<int>& testBlock,
    const Array<int>& unkID,
    const Array<int>& unkBlock,
    const Array<RCP<ElementIntegral> >& integrals,
    const Array<int>& resultIndices,
    const Array<MultipleDeriv>& derivs,
    int verb);


  /** Indicate whether this is a group of two-forms */
  bool isTwoForm() const {return order_==2;}

  /** Indicate whether this is a group of one-forms */
  bool isOneForm() const {return order_==1;}

  /** Indicate whether this is a group of zero-forms */
  bool isZeroForm() const {return order_==0;}

  /** Return the number of rows in the local matrices or vectors
   * computed by this integral group */
  int nTestNodes() const {return nTestNodes_;}

  /** Return the number of columns in the local matrices 
   * computed by this integral group */
  int nUnkNodes() const {return nUnkNodes_;}

  /** Return the test functions using this integral group */
  const Array<int>& testID() const {return testID_;}

  /** Return the unknown functions using this integral group */
  const Array<int>& unkID() const {return unkID_;}

  /** Return the block numbers for the test functions */
  const Array<int>& testBlock() const {return testBlock_;}

  /** Return the block numbers for the unk functions */
  const Array<int>& unkBlock() const {return unkBlock_;}

  /** Return the multivector column indices for the one-forms */
  const Array<int>& mvIndices() const {return mvIndices_;}
  

  /** Whether the group requires transformations based on a maximal cofacet */
  IntegrationCellSpecifier usesMaximalCofacets() const 
    {return requiresMaximalCofacet_;}

  /** Array specifying which terms need maximal cofacets */
  const Array<int>& termUsesMaximalCofacets() const
    {return termUsesMaximalCofacets_;}

  /** \brief Extract basis used for test functions in this integral group */
  const BasisFamily & getTestBasis() const { return integrals_[0]->getTestBasis(); }
    
  /** \brief Extract basis used for unknown functions in this integral group */
  const BasisFamily & getUnknownBasis() const { return integrals_[0]->getUnknownBasis(); }

  /** Evaluate this integral group */
  bool evaluate(const CellJacobianBatch& JTrans,
    const CellJacobianBatch& JVol,
    const Array<int>& isLocalFlag,
    const Array<int>& facetNum, 
    const RCP<Array<int> >& cellLIDs,
    const Array<RCP<EvalVector> >& vectorCoeffs,
    const Array<double>& constantCoeffs,
    RCP<Array<double> >& A) const ;


  /** */
  int integrationVerb() const {return integrationVerb_;}
    
  /** */
  int transformVerb() const {return transformVerb_;}

  /** */
  const Array<MultipleDeriv>& derivs() const 
    {return derivs_;}



private:
  /** */
  int findIntegrationVerb(const Array<RCP<ElementIntegral> >& integrals) const ;
  /** */
  int findTransformVerb(const Array<RCP<ElementIntegral> >& integrals) const ;
  /** */
  int integrationVerb_;

  /** */
  int transformVerb_;
      
  /** */
  int order_;

  /** */
  int nTestNodes_;

  /** */
  int nUnkNodes_;

  /** */
  Array<int> testID_;

  /** */
  Array<int> unkID_;

  /** */
  Array<int> testBlock_;

  /** */
  Array<int> unkBlock_;

  /** */
  Array<int> mvIndices_;

  /** */
  Array<RCP<ElementIntegral> > integrals_;

  /** */
  Array<int> resultIndices_;

  /** */
  Array<int> termUsesMaximalCofacets_;

  /** */
  IntegrationCellSpecifier requiresMaximalCofacet_;

  /** */
  Array<MultipleDeriv> derivs_;
};


}


#endif
