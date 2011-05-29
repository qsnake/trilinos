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

#ifndef SUNDANCE_SUBSETMANAGER_H
#define SUNDANCE_SUBSETMANAGER_H

#include "SundanceDefs.hpp"
#include "SundanceMap.hpp"
#include "SundanceSet.hpp"
#include "SundanceCellFilter.hpp"
#include "SundanceCellFilterBase.hpp"

namespace Sundance
{
/** 
 * SubsetManager keeps maps from cell filters to their subsets. This is
 * a device to avoid circular references between subset and superset.
 */

class SubsetManager
{
public:
  /** */
  static void registerSubset(const CellFilter& filter, 
    const CellFilter& subset)
    {
      if (!subsetMap().containsKey(filter))
        subsetMap().put(filter, Sundance::Set<CellFilter>());
      subsetMap().get(filter).put(subset);
    }

  /** */
  static void registerDisjoint(const CellFilter& filter, 
    const CellFilter& subset)
    {
      if (!disjointMap().containsKey(filter))
        disjointMap().put(filter, Sundance::Set<CellFilter>());
      disjointMap().get(filter).put(subset);
    }

  /** */
  static void registerLabeledSubset(const CellFilter& filter, 
    int label, const CellFilter& subset)
    {
      if (!labeledMap().containsKey(filter))
        labeledMap().put(filter, Sundance::Map<int, CellFilter>());
      labeledMap().get(filter).put(label, subset);
    }


  /** */
  static const Sundance::Set<CellFilter>& getSubsets(const CellFilter& filter)
    {
      if (!subsetMap().containsKey(filter))
        subsetMap().put(filter, Sundance::Set<CellFilter>());
      return subsetMap().get(filter);
    }

  /** */
  static const Sundance::Set<CellFilter>& getDisjoints(const CellFilter& filter)
    {
      if (!disjointMap().containsKey(filter))
        disjointMap().put(filter, Sundance::Set<CellFilter>());
      return disjointMap().get(filter);
    }

  /** */
  static const Sundance::Map<int, CellFilter>& 
  getLabeledSubsets(const CellFilter& filter)
    {
      if (!labeledMap().containsKey(filter))
        labeledMap().put(filter, Sundance::Map<int, CellFilter>());
      return labeledMap().get(filter);
    }

private:
  
  /** */
  static Sundance::Map<CellFilter, Sundance::Set<CellFilter> >& subsetMap()
    {
      static Sundance::Map<CellFilter, Sundance::Set<CellFilter> > rtn;
      return rtn;
    }
  
  /** */
  static Sundance::Map<CellFilter, Sundance::Set<CellFilter> >& disjointMap()
    {
      static Sundance::Map<CellFilter, Sundance::Set<CellFilter> > rtn;
      return rtn;
    }
  
  /** */
  static Sundance::Map<CellFilter, Sundance::Map<int, CellFilter> >& labeledMap()
    {
      static Sundance::Map<CellFilter, Sundance::Map<int, CellFilter> > rtn;
      return rtn;
    }
};

}



#endif
