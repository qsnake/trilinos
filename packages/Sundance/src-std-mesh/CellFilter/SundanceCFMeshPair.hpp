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

#ifndef SUNDANCE_CFMESHPAIR_H
#define SUNDANCE_CFMESHPAIR_H

#include "SundanceDefs.hpp"
#include "SundanceCellFilter.hpp"
#include "SundanceCellSet.hpp"
#include "SundanceCellFilterBase.hpp"

namespace Sundance
{
using namespace Teuchos;
  

/** 
 * 
 */
class CFMeshPair 
{
public:
  /** */
  CFMeshPair();
  /** */
  CFMeshPair(const CellFilter& cf,
    const Mesh& mesh,
    const Set<int>& funcs);

  /** */
  bool operator<(const CFMeshPair& other) const ;

  /** */
  bool isEmpty() const ;

  /** */
  const CellFilter& filter() const {return filter_;}

  /** */
  const Mesh& mesh() const {return mesh_;}

  /** */
  const CellSet& cellSet() const {return cellSet_;}

  /** */
  const Set<int>& funcs() const {return funcs_;}

  /** */
  CFMeshPair setMinus(const CFMeshPair& other) const ;

  /** */
  CFMeshPair intersection(const CFMeshPair& other) const ;

private:
  CellFilter filter_;
  Mesh mesh_;
  CellSet cellSet_;
  Set<int> funcs_;
};

/** */
Array<CFMeshPair> resolvePair(const CFMeshPair& a, 
  const CFMeshPair& b);
/** */
Array<CFMeshPair> resolveSets(const Array<CFMeshPair>& s);

/** */
Array<CFMeshPair>
findDisjointFilters(const Array<CellFilter>& filters,
  const Array<Set<int> >& funcs,
  const Mesh& mesh);


inline std::ostream& operator<<(std::ostream& os, 
  const CFMeshPair& c)
{
  os << c.filter().getCells(c.mesh());
  return os;
}


}



#endif
