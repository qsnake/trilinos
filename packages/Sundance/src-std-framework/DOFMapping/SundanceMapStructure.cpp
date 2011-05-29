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

#include "SundanceMapStructure.hpp"
#include "SundanceTabs.hpp"
#include "SundanceExceptions.hpp"
#include "SundanceBasisDOFTopologyBase.hpp"
#include "SundanceOut.hpp"


namespace Sundance
{
using namespace Teuchos;
using std::endl;
using std::setw;

MapStructure::MapStructure(int nTotalFuncs,
                           const Array<RCP<BasisDOFTopologyBase> >& bases,
                           const Array<Array<int> >& funcs)
{
  init(nTotalFuncs, bases, funcs);
}


MapStructure::MapStructure(int nTotalFuncs,
                           const RCP<BasisDOFTopologyBase>& bases,
                           const Array<Array<int> >& funcs)
{
  init(nTotalFuncs, replicate(bases, funcs.size()), funcs);
}

MapStructure::MapStructure(int nTotalFuncs,
                           const RCP<BasisDOFTopologyBase>& bases)
{
  Array<int> f(nTotalFuncs);
  for (int i=0; i<nTotalFuncs; i++) f[i] = i;
  init(nTotalFuncs, tuple(bases), tuple(f));
}



void MapStructure::init(int nTotalFuncs,
                        const Array<RCP<BasisDOFTopologyBase> >& bases,
                        const Array<Array<int> >& funcs)
{
  bases_ = bases;
  funcs_ = funcs;
  chunkForFuncID_.resize(nTotalFuncs);
  indexForFuncID_.resize(nTotalFuncs);

  TEST_FOR_EXCEPTION(bases.size() != funcs.size(), InternalError,
                     "mismatched number of basis chunks=" << bases.size()
                     << " and number of function chunks=" << funcs.size());

  for (int f=0; f<indexForFuncID_.size(); f++) 
    {
      indexForFuncID_[f] = -1;
      chunkForFuncID_[f] = -1;
    }

  for (int b=0; b<funcs_.size(); b++)
    {
      for (int f=0; f<funcs_[b].size(); f++)
        {
          int fid = funcs_[b][f];
          TEST_FOR_EXCEPTION(fid >= nTotalFuncs, InternalError,
                             "bad funcID=" << fid 
                             << ". nTotalFuncs=" << nTotalFuncs);
          indexForFuncID_[fid] = f;
          chunkForFuncID_[fid] = b;
        }
    }
}

int MapStructure::chunkForFuncID(int funcID) const 
{
  int rtn = chunkForFuncID_[funcID];
  TEST_FOR_EXCEPTION(rtn < 0, InternalError,
                     "funcID=" << funcID << " not defined in map chunk."
                     " The functions defined there are " 
                     << funcs_ << ". The most likely cause of this error is "
                     "that you are trying to access a discrete function on "
                     "subdomain for which it is not defined.");
  return rtn;
}

int MapStructure::indexForFuncID(int funcID) const 
{
  int rtn = indexForFuncID_[funcID];
  TEST_FOR_EXCEPTION(rtn < 0, InternalError,
                     "funcID=" << funcID << " not defined in map chunk."
                     " The functions defined there are " 
                     << funcs_ << ". The most likely cause of this error is "
                     "that you are trying to access a discrete function on "
                     "subdomain for which it is not defined.");

  return rtn;
}


std::ostream& MapStructure::print(std::ostream& os) const
{
  Tabs tab;
  os << tab << "Map structure: nChunks=" << numBasisChunks() << std::endl;
  for (int c=0; c<numBasisChunks(); c++) 
  {
    Tabs tab1;
    os << tab1 << "chunk " << c << " funcIDs=" << funcs(c) << std::endl;
  }
  os << tab << "## end map structure" << std::endl;
  return os;
}


Array<RCP<BasisDOFTopologyBase> > replicate(
  const RCP<BasisDOFTopologyBase>& model,
  int n)
{
  Array<RCP<BasisDOFTopologyBase> > rtn(n);
  for (int i=0; i<n; i++) rtn[i] = model;
  return rtn;
}

} // end namespace Sundance

