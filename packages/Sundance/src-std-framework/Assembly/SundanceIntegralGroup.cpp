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

#include "SundanceIntegralGroup.hpp"
#include "SundanceRefIntegral.hpp"
#include "SundanceQuadratureIntegral.hpp"
#include "SundanceMaximalQuadratureIntegral.hpp"
#include "SundanceCurveQuadratureIntegral.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"
#include "Teuchos_TimeMonitor.hpp"

using namespace Sundance;
using namespace Teuchos;


static Time& integrationTimer() 
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("integration"); 
  return *rtn;
}


IntegralGroup
::IntegralGroup(const Array<RCP<ElementIntegral> >& integrals,
  const Array<int>& resultIndices,
  int verb)
  : integrationVerb_(findIntegrationVerb(integrals)),
    transformVerb_(findTransformVerb(integrals)),
    order_(0),
    nTestNodes_(0),
    nUnkNodes_(0),
    testID_(),
    unkID_(),
    testBlock_(),
    unkBlock_(),
    mvIndices_(),
    integrals_(integrals),
    resultIndices_(resultIndices),
    termUsesMaximalCofacets_(integrals_.size()),
    requiresMaximalCofacet_(SomeTermsNeedCofacets),
    derivs_()
{
  Tabs tab;
  SUNDANCE_MSG2(verb, tab << "forming 0-form integral group");
  bool allReqMaximalCofacets = true;
  bool noneReqMaximalCofacets = true;
  bool someReqMaximalCofacets = false;

  for (int i=0; i<integrals_.size(); i++)
  {
    Tabs tab1;
    SUNDANCE_MSG2(verb, tab1 << "integral #" << i << ", numFacetCases="
      << integrals[i]->nFacetCases());
    if (integrals[i]->nFacetCases() > 1) 
    {
      Tabs tab2;
      someReqMaximalCofacets = true;
      noneReqMaximalCofacets = false;
      termUsesMaximalCofacets_[i] = true;
      SUNDANCE_MSG2(verb, tab2 << "I need maximal cofacets");
    }
    else
    {
      Tabs tab2;
      allReqMaximalCofacets = false;
      termUsesMaximalCofacets_[i] = false;
      SUNDANCE_MSG2(verb, tab2 << "I do not need maximal cofacets");
    }
  }

  if (allReqMaximalCofacets) 
  {
    requiresMaximalCofacet_ = AllTermsNeedCofacets;
  }
  else if (noneReqMaximalCofacets) 
  {
    requiresMaximalCofacet_ = NoTermsNeedCofacets;
  }
  
  Tabs tab1;
  SUNDANCE_MSG2(verb, tab1 << "result=" << requiresMaximalCofacet_);
}




IntegralGroup
::IntegralGroup(const Array<int>& testID,
  const Array<int>& testBlock,
  const Array<int>& mvIndices,
  const Array<RCP<ElementIntegral> >& integrals,
  const Array<int>& resultIndices,
  const Array<MultipleDeriv>& derivs,
  int verb)
  : integrationVerb_(findIntegrationVerb(integrals)),
    transformVerb_(findTransformVerb(integrals)),
    order_(1),
    nTestNodes_(0),
    nUnkNodes_(0),
    testID_(testID),
    unkID_(),
    testBlock_(testBlock),
    unkBlock_(),
    mvIndices_(mvIndices),
    integrals_(integrals),
    resultIndices_(resultIndices),
    termUsesMaximalCofacets_(integrals_.size()),
    requiresMaximalCofacet_(SomeTermsNeedCofacets),
    derivs_(derivs)
{
  Tabs tab;
  SUNDANCE_MSG2(verb, tab << "forming 1-form integral group");
  bool allReqMaximalCofacets = true;
  bool noneReqMaximalCofacets = true;
  bool someReqMaximalCofacets = false;

  for (int i=0; i<integrals.size(); i++)
  {
    int nt = integrals[i]->nNodesTest();
    if (i > 0) 
    {
      TEST_FOR_EXCEPTION(nt != nTestNodes_, InternalError,
        "IntegralGroup ctor detected integrals with inconsistent numbers of test nodes");
    }
    Tabs tab1;
    SUNDANCE_MSG2(verb, tab1 << "integral #" << i << ", numFacetCases="
      << integrals[i]->nFacetCases());
    nTestNodes_ = nt;
    if (integrals[i]->nFacetCases() > 1) 
    {
      Tabs tab2;
      SUNDANCE_MSG2(verb, tab2 << "I need maximal cofacets");
      someReqMaximalCofacets = true;
      noneReqMaximalCofacets = false;
      termUsesMaximalCofacets_[i] = true;
    }
    else
    {
      Tabs tab2;
      SUNDANCE_MSG2(verb, tab2 << "I do not need maximal cofacets");
      allReqMaximalCofacets = false;
      termUsesMaximalCofacets_[i] = false;
    }
  }

  if (allReqMaximalCofacets) 
  {
    requiresMaximalCofacet_ = AllTermsNeedCofacets;
  }
  else if (noneReqMaximalCofacets) 
  {
    requiresMaximalCofacet_ = NoTermsNeedCofacets;
  }
  
  Tabs tab1;
  SUNDANCE_MSG2(verb, tab1 << "result=" << requiresMaximalCofacet_);
}

IntegralGroup
::IntegralGroup(const Array<int>& testID,
  const Array<int>& testBlock,
  const Array<int>& unkID,
  const Array<int>& unkBlock,
  const Array<RCP<ElementIntegral> >& integrals,
  const Array<int>& resultIndices,
  const Array<MultipleDeriv>& derivs,
  int verb)
  : integrationVerb_(findIntegrationVerb(integrals)),
    transformVerb_(findTransformVerb(integrals)),
    order_(2),
    nTestNodes_(0),
    nUnkNodes_(0),
    testID_(testID),
    unkID_(unkID),
    testBlock_(testBlock),
    unkBlock_(unkBlock),
    mvIndices_(),
    integrals_(integrals),
    resultIndices_(resultIndices),
    termUsesMaximalCofacets_(integrals_.size()),
    requiresMaximalCofacet_(SomeTermsNeedCofacets),
    derivs_(derivs)
{
  Tabs tab;
  SUNDANCE_MSG2(verb, tab << "forming 2-form integral group");
  bool allReqMaximalCofacets = true;
  bool noneReqMaximalCofacets = true;
  bool someReqMaximalCofacets = false;

  for (int i=0; i<integrals.size(); i++)
  {
    Tabs tab1;
    SUNDANCE_MSG2(verb, tab1 << "integral #" << i << ", numFacetCases="
      << integrals[i]->nFacetCases());
    int nt = integrals[i]->nNodesTest();
    int nu = integrals[i]->nNodesUnk();
    if (i > 0) 
    {
      TEST_FOR_EXCEPTION(nt != nTestNodes_, InternalError,
        "IntegralGroup ctor detected integrals with inconsistent numbers of test nodes");
      TEST_FOR_EXCEPTION(nu != nUnkNodes_, InternalError,
        "IntegralGroup ctor detected integrals with inconsistent numbers of unk nodes");
    }
    nTestNodes_ = nt;
    nUnkNodes_ = nu;

    if (integrals[i]->nFacetCases() > 1) 
    {
      someReqMaximalCofacets = true;
      noneReqMaximalCofacets = false;
      termUsesMaximalCofacets_[i] = true;
      Tabs tab2;
      SUNDANCE_MSG2(verb, tab2 << "I need maximal cofacets");
    }
    else
    {
      allReqMaximalCofacets = false;
      termUsesMaximalCofacets_[i] = false;
      Tabs tab2;
      SUNDANCE_MSG2(verb, tab2 << "I do not need maximal cofacets");
    }
  }

  if (allReqMaximalCofacets) 
  {
    requiresMaximalCofacet_ = AllTermsNeedCofacets;
  }
  else if (noneReqMaximalCofacets) 
  {
    requiresMaximalCofacet_ = NoTermsNeedCofacets;
  }
  
  Tabs tab1;
  SUNDANCE_MSG2(verb, tab1 << "result=" << requiresMaximalCofacet_);
}


bool IntegralGroup
::evaluate(const CellJacobianBatch& JTrans,
  const CellJacobianBatch& JVol,
  const Array<int>& isLocalFlag, 
  const Array<int>& facetIndex, 
  const RCP<Array<int> >& cellLIDs,
  const Array<RCP<EvalVector> >& vectorCoeffs,
  const Array<double>& constantCoeffs,
  RCP<Array<double> >& A) const
{
  TimeMonitor timer(integrationTimer());
  Tabs tab0(0);


  SUNDANCE_MSG1(integrationVerb(), tab0 << "evaluating integral group with "
    << integrals_.size() << " integrals");

  SUNDANCE_MSG3(integrationVerb(), 
    tab0 << "num integration cells = " << JVol.numCells());
  SUNDANCE_MSG3(integrationVerb(), 
    tab0 << "num nodes in output = " << integrals_[0]->nNodes());

  /* initialize the return vector */
  if (integrals_[0]->nNodes() == -1) A->resize(1);
  else A->resize(JVol.numCells() * integrals_[0]->nNodes());
  double* aPtr = &((*A)[0]);
  int n = A->size();
  for (int i=0; i<n; i++) aPtr[i] = 0.0;

  SUNDANCE_MSG5(integrationVerb(), tab0 << "begin A=");
  if (integrationVerb() >=5) writeTable(Out::os(), tab0, *A, 6);

  /* do the integrals */
  for (int i=0; i<integrals_.size(); i++)
  {
    Tabs tab1;
    SUNDANCE_MSG1(integrationVerb(), tab1 << "group member i=" << i 
      << " of " << integrals_.size());
    Tabs tab2;

    const RefIntegral* ref 
      = dynamic_cast<const RefIntegral*>(integrals_[i].get());
    const QuadratureIntegral* quad 
      = dynamic_cast<const QuadratureIntegral*>(integrals_[i].get());
    const MaximalQuadratureIntegral* maxQuad 
      = dynamic_cast<const MaximalQuadratureIntegral*>(integrals_[i].get());
    const CurveQuadratureIntegral* curveQuad
      = dynamic_cast<const CurveQuadratureIntegral*>(integrals_[i].get());

    if (ref!=0)
    {
      SUNDANCE_MSG2(integrationVerb(),
        tab2 << "Integrating term group " << i 
        << " by transformed reference integral");
      double f = constantCoeffs[resultIndices_[i]];
      SUNDANCE_MSG2(integrationVerb(),
        tab2 << "Coefficient is " << f);
      ref->transform(JTrans, JVol, isLocalFlag, facetIndex, cellLIDs , f, A);
    }
    else if (quad != 0)
    {
      SUNDANCE_MSG2(integrationVerb(),
        tab2 << "Integrating term group " << i 
        << " by quadrature");
          
      TEST_FOR_EXCEPTION(vectorCoeffs[resultIndices_[i]]->length()==0,
        InternalError,
        "zero-length coeff vector detected in "
        "quadrature integration branch of "
        "IntegralGroup::evaluate(). std::string value is ["
        << vectorCoeffs[resultIndices_[i]]->str()
        << "]");

      Tabs tab3;
      SUNDANCE_MSG3(integrationVerb(),
        tab3 << "coefficients are " <<  vectorCoeffs[resultIndices_[i]]->str());

      const double* const f = vectorCoeffs[resultIndices_[i]]->start();
      quad->transform(JTrans, JVol, isLocalFlag, facetIndex, cellLIDs , f, A);
    }
    else if (maxQuad != 0)
    {
      SUNDANCE_MSG2(integrationVerb(),
        tab2 << "Integrating term group " << i 
        << " by quadrature");
          
      TEST_FOR_EXCEPTION(vectorCoeffs[resultIndices_[i]]->length()==0,
        InternalError,
        "zero-length coeff vector detected in "
        "quadrature integration branch of "
        "IntegralGroup::evaluate(). std::string value is ["
        << vectorCoeffs[resultIndices_[i]]->str()
        << "]");

      Tabs tab3;
      SUNDANCE_MSG3(integrationVerb(),
        tab3 << "coefficients are " <<  vectorCoeffs[resultIndices_[i]]->str());

      const double* const f = vectorCoeffs[resultIndices_[i]]->start();
      maxQuad->transform(JTrans, JVol, isLocalFlag, facetIndex, cellLIDs , f, A);
    }
    else if (curveQuad != 0)
    {
        SUNDANCE_MSG2(integrationVerb(),
          tab2 << "Integrating term group " << i
          << " by curve integral (quadrature by default) , result index: " << resultIndices_[i]);

        double f_const = 0.0;
        if (constantCoeffs.size() > resultIndices_[i]){
        	f_const = constantCoeffs[resultIndices_[i]];
        }

        SUNDANCE_MSG2(integrationVerb(),
          tab2 << "Coefficient is " << f_const);

        // set this
        if (vectorCoeffs.size() > resultIndices_[i]){
            Tabs tab3;
            double* const f = vectorCoeffs[resultIndices_[i]]->start();
        	SUNDANCE_MSG3(integrationVerb(),
        			tab3 << "coefficients are " <<  vectorCoeffs[resultIndices_[i]]->str());
            curveQuad->transform(JTrans, JVol, isLocalFlag, facetIndex, cellLIDs , f_const , f , A);
        } else{
            const double* f_null = 0;
            curveQuad->transform(JTrans, JVol, isLocalFlag, facetIndex, cellLIDs , f_const , f_null , A);
        }

    }
    else
    {
      TEST_FOR_EXCEPT(1);
    }

    SUNDANCE_MSG4(integrationVerb(),
      tab1 << "i=" << i << " integral values=");
    if (integrationVerb() >=4) writeTable(Out::os(), tab1, *A, 6);
  }
  SUNDANCE_MSG1(integrationVerb(), tab0 << "done integral group");

  return true;
}


int IntegralGroup::findIntegrationVerb(const Array<RCP<ElementIntegral> >& integrals) const
{
  int rtn = 0;
  for (int i=0; i<integrals.size(); i++)
  {
    rtn = std::max(rtn, integrals[i]->integrationVerb());
  }
  return rtn;
}


int IntegralGroup::findTransformVerb(const Array<RCP<ElementIntegral> >& integrals) const
{
  int rtn = 0;
  for (int i=0; i<integrals.size(); i++)
  {
    rtn = std::max(rtn, integrals[i]->transformVerb());
  }
  return rtn;
}
