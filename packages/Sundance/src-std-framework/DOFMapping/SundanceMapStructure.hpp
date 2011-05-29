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

#ifndef SUNDANCE_MAPSTRUCTURE_H
#define SUNDANCE_MAPSTRUCTURE_H


#include "SundanceDefs.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_RefCountPtr.hpp"

namespace Sundance
{
using namespace Teuchos;

class BasisDOFTopologyBase;

/** 
 * 
 */
class MapStructure
{
public:
  /** */
  MapStructure(int nTotalFuncs,
    const Array<RCP<BasisDOFTopologyBase> >& bases,
    const Array<Array<int> >& funcs);
  /** */
  MapStructure(int nTotalFuncs,
    const RCP<BasisDOFTopologyBase>& basis,
    const Array<Array<int> >& funcs);
  /** */
  MapStructure(int nTotalFuncs,
    const RCP<BasisDOFTopologyBase>& basis);

  /** */
  int numBasisChunks() const {return bases_.size();}

  /** */
  const RCP<BasisDOFTopologyBase>& basis(int basisChunk) const
    {return bases_[basisChunk];}

  /** */
  int numFuncs(int basisChunk) const 
    {return funcs_[basisChunk].size();}

  /** */
  const Array<int>& funcs(int basisChunk) const 
    {return funcs_[basisChunk];}

  /** */
  int chunkForFuncID(int funcID) const ;

  /** */
  int indexForFuncID(int funcID) const ;

  /** */
  std::ostream& print(std::ostream& os) const ;

private:
  /** */
  void init(int nTotalFuncs,
    const Array<RCP<BasisDOFTopologyBase> >& bases,
    const Array<Array<int> >& funcs);

  Array<RCP<BasisDOFTopologyBase> > bases_;
  Array<Array<int> > funcs_;
  Array<int> chunkForFuncID_;
  Array<int> indexForFuncID_;
};

/** \relates BasisDOFTopologyBase */
Array<RCP<BasisDOFTopologyBase> > replicate(
  const RCP<BasisDOFTopologyBase>& model,
  int n);



/** \relates MapStructure */
inline std::ostream& operator<<(std::ostream& os,
  const MapStructure& m)
{
  return m.print(os);
}

}
#endif
