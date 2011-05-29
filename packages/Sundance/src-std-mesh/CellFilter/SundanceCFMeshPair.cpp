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

#include "SundanceCFMeshPair.hpp"
#include "SundanceTabs.hpp"
#include "SundanceExplicitCellSet.hpp"
#include "SundanceBinaryCellFilter.hpp"
#include "SundanceSubsetCellFilter.hpp"
#include "SundanceLabelCellPredicate.hpp"
#include "SundanceNullCellFilterStub.hpp"
#include "SundanceNullCellFilter.hpp"

#include "Teuchos_Time.hpp"
#include "Teuchos_TimeMonitor.hpp"

using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;


static Time& csPartitionTimer() 
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("cell set partitioning"); 
  return *rtn;
}


CFMeshPair::CFMeshPair()
  : filter_(),
    mesh_(),
    cellSet_(),
    funcs_()
{}

CFMeshPair::CFMeshPair(const CellFilter& cf,
                       const Mesh& mesh,
                       const Set<int>& funcs)
  : filter_(cf),
    mesh_(mesh),
    cellSet_(),
    funcs_(funcs)
{
  if (filter_.ptr().get() != 0) cellSet_ = filter_.getCells(mesh);
}

bool CFMeshPair::operator<(const CFMeshPair& other) const
{
  if (isEmpty()) return true;
  if (other.isEmpty()) return false;

  TEST_FOR_EXCEPTION(mesh_.id() != other.mesh_.id(),
                     RuntimeError,
                     "mismatched meshes!");

  return cellSet_ < other.cellSet_;
}

bool CFMeshPair::isEmpty() const
{
  return filter_.ptr().get() == 0 
    || cellSet_.ptr().get() == 0
    || cellSet_.begin()==cellSet_.end();
}

CFMeshPair CFMeshPair::setMinus(const CFMeshPair& other) const
{

  if (isEmpty()) return CFMeshPair();
  if (other.isEmpty()) return *this;

  TEST_FOR_EXCEPTION(mesh().id() != other.mesh().id(),
                     RuntimeError,
                     "mismatched meshes!");

  CellFilter diff = filter() - other.filter();
  return CFMeshPair(diff, mesh(), funcs_);
}

CFMeshPair CFMeshPair::intersection(const CFMeshPair& other) const
{
  if (isEmpty() || other.isEmpty()) return CFMeshPair();

  TEST_FOR_EXCEPTION(mesh().id() != other.mesh().id(),
                     RuntimeError,
                     "mismatched meshes!");

  CellFilter inter = filter().intersection(other.filter());

  return CFMeshPair(inter, mesh(), funcs_.setUnion(other.funcs_));
}


namespace Sundance
{
  Array<CFMeshPair> resolvePair(const CFMeshPair& a, 
                                const CFMeshPair& b)
  {
    Tabs tab0;

    CFMeshPair inter = a.intersection(b);
    CFMeshPair aMinusB = a.setMinus(b);
    CFMeshPair bMinusA = b.setMinus(a);

    return tuple(bMinusA, inter, aMinusB);
  }

  
  Array<CFMeshPair> resolveSets(const Array<CFMeshPair>& s)
  {

    if (s.size() == 1) return s;

    Array<CFMeshPair> T = tuple(s[0]);
    
    for (int i=0; i<s.size(); i++)
      {
        CFMeshPair A = s[i];
        Array<CFMeshPair> TN;
        for (int j=0; j<T.size(); j++)
          {
            Array<CFMeshPair> p = resolvePair(A, T[j]);
            if (!p[0].isEmpty())
              {
                TN.append(p[0]);
              }
            if (!p[1].isEmpty())
              {
                TN.append(p[1]);
              }
            A = p[2];
          }
        if (!A.isEmpty())
          {
            TN.append(A);
          }
        T = TN;
      }
    return T;
  }


  Array<CFMeshPair>
  findDisjointFilters(const Array<CellFilter>& filters,
                      const Array<Set<int> >& funcs,
                      const Mesh& mesh)
  {
    TimeMonitor timer(csPartitionTimer() );
    Array<CFMeshPair> cf(filters.size());

    TEST_FOR_EXCEPT(filters.size() != funcs.size());
    for (int i=0; i<filters.size(); i++)
      {
        cf[i] = CFMeshPair(filters[i], mesh, funcs[i]);
      }
    return resolveSets(cf);
  }

}
