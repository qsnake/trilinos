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

#include "SundanceBasisFamily.hpp"
#include "SundanceUnknownFunction.hpp"
#include "SundanceTestFunction.hpp"
#include "SundanceDiscreteFunction.hpp"
#include "SundanceUnknownFunctionData.hpp"
#include "SundanceTestFunctionData.hpp"
#include "SundanceDiscreteFunctionData.hpp"
#include "SundanceLagrange.hpp"
#include "SundanceEdgeLocalizedBasis.hpp"


using namespace Sundance;
using namespace Teuchos;





int BasisFamily::order() const 
{
  return ptr()->order();
}

int BasisFamily::dim() const 
{
  return ptr()->dim();
}

bool BasisFamily::operator==(const BasisFamily& other) const
{
  return !(*this < other || other < *this);
}

int BasisFamily::size(const Array<BasisFamily>& b)
{
  int rtn = 0;
  for (int i=0; i<b.size(); i++) rtn += b[i].dim();
  return rtn;
}

int BasisFamily::nReferenceDOFsWithFacets(const CellType& maximalCellType,
					  const CellType& cellType) const 
{
  return ptr()->nReferenceDOFsWithFacets(maximalCellType, cellType);
}

int BasisFamily::nReferenceDOFsWithoutFacets(const CellType& maximalCellType,
					     const CellType& cellType) const 
{
  return ptr()->nReferenceDOFsWithoutFacets(maximalCellType, cellType);
}

BasisFamily BasisFamily::getBasis(const RCP<const CommonFuncDataStub>& funcData)
{
  /* First try to process assuming the input is an unknown function */
  const UnknownFunctionData* u 
    = dynamic_cast<const UnknownFunctionData*>(funcData.get());
  if (u != 0)
    {
      return u->basis()[0];
    }

  /* Next try to process assuming the input is a test function */
  const TestFunctionData* t 
    = dynamic_cast<const TestFunctionData*>(funcData.get());
  if (t != 0)
    {
      return t->basis()[0];
    }


  /* Next try to process assuming the input is a discrete function */
  const DiscreteFunctionData* d
    = dynamic_cast<const DiscreteFunctionData*>(funcData.get());
  if (d != 0)
    {
      return d->discreteSpace().basis()[0];
    }

  TEST_FOR_EXCEPTION(u==0 && t==0 && d==0, RuntimeError,
		     "BasisFamily::getBasis() argument is not a recognized "
		     "type of function data");
  return BasisFamily();
  
}



RCP<BasisDOFTopologyBase> BasisFamily::getBasisTopology(const RCP<const CommonFuncDataStub>& funcData)
{
  BasisFamily b = getBasis(funcData);
  TEST_FOR_EXCEPT(b.ptr().get()==0);

  return rcp_dynamic_cast<BasisDOFTopologyBase>(b.ptr());
}


void BasisFamily::getConstrainsForHNDoF( const int indexInParent,
					 const int maxCellDim,
					 const int maxNrChild,
					 const int facetDim,
					 const int facetIndex,
					 const int nodeIndex,
					 Array<int>& localDoFs,
					 Array<int>& parentFacetDim,
					 Array<int>& parentFacetIndex,
					 Array<int>& parentFacetNode,
					 Array<double>& coefs
					 ) const
{
  Tabs tab;
  int verb = 4;
  SUNDANCE_MSG3( verb ,tab << "BasisFamily::getConstrainsForHNDoF IN: indexInParent:" << indexInParent
		 << "  maxCellDim:" << maxCellDim << " maxNrChild:" << maxNrChild
		 << " facetDim:" << facetDim << "  facetIndex:" << facetIndex << " nodeIndex:" << nodeIndex);
  ptr()->getConstrainsForHNDoF( indexInParent, maxCellDim ,
				maxNrChild , facetDim, facetIndex, nodeIndex, localDoFs, parentFacetDim ,
				parentFacetIndex , parentFacetNode , coefs );
  SUNDANCE_MSG3( verb , tab << "BasisFamily::getConstrainsForHNDoF OUT: localDoFs:" << localDoFs
		 << " coefs:" << coefs);
  SUNDANCE_MSG3( verb , tab << "BasisFamily::getConstrainsForHNDoF OUT: parentFacetDim:" << parentFacetDim );
  SUNDANCE_MSG3( verb , tab << "BasisFamily::getConstrainsForHNDoF OUT: parentFacetIndex:" << parentFacetIndex );
  SUNDANCE_MSG3( verb , tab << "BasisFamily::getConstrainsForHNDoF OUT: parentFacetNode:" << parentFacetNode );
}


void BasisFamily::refEval(
			  const CellType& cellType,
			  const Array<Point>& pts,
			  const SpatialDerivSpecifier& deriv,
			  Array<Array<Array<double> > >& result,
			  int verbosity) const
{
  using std::setw;
  Tabs tab;
  SUNDANCE_MSG3(verbosity, tab << "evaluating basis " << *this 
		<< " with spatial derivative " << deriv);
  ptr()->refEval(cellType, pts, deriv, result, verbosity);
  std::string f = deriv.toString()+ "[phi_n]";

  if (verbosity >= 4)
    {
      Tabs tab1;
      for (int q=0; q<pts.size(); q++)
	{
	  Tabs tab2;
	  Out::os() << tab1 << "evaluation point = " << pts[q] << std::endl;
	  Out::os() << tab2 << setw(5) << "n";
	  int nComps = result.size();
	  if (nComps == 1)
	    {
	      Out::os() << setw(20) <<  f << std::endl;
	    }
	  else
	    {
	      for (int d=0; d<nComps; d++)
		{
		  std::string fd = f + "[" + Teuchos::toString(d) + "]";
		  Out::os() << setw(20) <<  fd;
		}
	      Out::os() << std::endl;
	    }
	  for (int n=0; n<result[0][q].size(); n++)
	    {
	      Out::os() << tab2 << setw(5) << n;
	      for (int d=0; d<nComps; d++)
		{
		  Out::os() << setw(20) <<  result[d][q][n];
		}
	      Out::os() << std::endl;
	    }
	}
    }
}


namespace Sundance
{

  Array<std::pair<int, int> > 
  vectorDimStructure(const Array<BasisFamily>& basis)
  {
    Array<std::pair<int, int> > rtn;
    for (int i=0; i<basis.size(); i++) 
      {
	rtn.append(std::pair<int, int>(basis[i].tensorOrder(), basis[i].dim()));
      }
    return rtn;
  }


  Array<std::pair<int, int> > vectorDimStructure(const BasisFamily& basis)
  {
    return vectorDimStructure(tuple(basis));
  }

  bool basisRestrictableToBoundary(const BasisFamily& b)
  {
    const Lagrange* lag = dynamic_cast<const Lagrange*>(b.ptr().get());
    const EdgeLocalizedBasis* elb = dynamic_cast<const EdgeLocalizedBasis*>(b.ptr().get());
    return lag != 0 || elb != 0;
  }

}
