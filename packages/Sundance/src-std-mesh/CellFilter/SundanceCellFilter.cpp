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

#include "SundanceCellFilter.hpp"
#include "SundanceCellFilterBase.hpp"
#include "SundanceExplicitCellSet.hpp"
#include "SundanceBinaryCellFilter.hpp"
#include "SundanceSubsetCellFilter.hpp"
#include "SundanceLabelCellPredicate.hpp"
#include "SundanceNullCellFilterStub.hpp"
#include "SundanceNullCellFilter.hpp"
#include "SundanceTabs.hpp"
#include "SundanceSubsetManager.hpp"

using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;

bool CellFilter::isNullCellFilter() const 
{
  return dynamic_cast<const NullCellFilterStub*>(ptr().get()) != 0;
}

bool CellFilter::isNull() const
{
  return ptr().get() == 0 || isNullCellFilter();
}

void CellFilter::setName(const std::string& name)
{
  nonConstCfbPtr()->setName(name);
}

CellSet CellFilter::getCells(const Mesh& mesh) const
{
  if (isNull() || isNullCellFilter())
  {
    return new ExplicitCellSet(mesh, -1, 
      NullCell);
  }
  return cfbPtr()->getCells(mesh);
}



int CellFilter::dimension(const Mesh& mesh) const
{
  if (isNullCellFilter())
  {
    return -1;
  }
  return cfbPtr()->dimension(mesh);
}



CellFilter CellFilter::operator+(const CellFilter& other) const 
{
  if (isNull())
  {
    return other;
  }
  else if (other.isNull())
  {
    return *this;
  }
  else
  {
    CellFilter rtn 
      = new BinaryCellFilter(*this, other, BinaryCellFilter::Union);
    rtn.registerSubset(*this);
    rtn.registerSubset(other);
    return rtn;
  }
}



CellFilter CellFilter::operator-(const CellFilter& other) const 
{
  if (other.isNull())
  {
    return *this;
  }
  else if (isKnownDisjointWith(other) || other.isKnownDisjointWith(*this))
  {
    return *this;
  }
  else if (isKnownSubsetOf(other))
  {
    CellFilter rtn;
    return rtn;
  }
  else if (*this == other)
  {
    CellFilter rtn;
    return rtn;
  }
  else
  {
    CellFilter rtn 
      = new BinaryCellFilter(*this, other, BinaryCellFilter::Difference);
    rtn.registerDisjoint(other);
    this->registerSubset(rtn);
    return rtn;
  }
}



CellFilter CellFilter::intersection(const CellFilter& other) const 
{
  if (isNull() || other.isNull())
  {
    CellFilter rtn;
    return rtn;
  }
  else if (isKnownDisjointWith(other) || other.isKnownDisjointWith(*this))
  {
    CellFilter rtn;
    return rtn;
  }
  else if (isKnownSubsetOf(other))
  {
    return *this;
  }
  else if (other.isKnownSubsetOf(*this))
  {
    return other;
  }
  else if (*this==other)
  {
    return *this;
  }
  else
  {
    CellFilter rtn 
      = new BinaryCellFilter(*this, other, BinaryCellFilter::Intersection);
    other.registerSubset(rtn);
    this->registerSubset(rtn);
    
    return rtn;
  }
}



CellFilter CellFilter::labeledSubset(int label) const
{
  CellPredicate pred = new LabelCellPredicate(label);
  CellFilter rtn = new SubsetCellFilter(*this, pred);
  this->registerLabeledSubset(label, rtn);
  this->registerSubset(rtn);

  return rtn;
}

CellFilter CellFilter::subset(const CellPredicate& pred) const
{
  CellFilter rtn = new SubsetCellFilter(*this, pred);
  this->registerSubset(rtn);
  return rtn;
}


CellFilter CellFilter::subset(const RCP<CellPredicateFunctorBase>& test) const
{
  CellFilter rtn = new SubsetCellFilter(*this, CellPredicate(test));
  this->registerSubset(rtn);
  return rtn;
}

bool CellFilter::isKnownSubsetOf(const CellFilter& other) const
{
  if (other.knownSubsets().contains(*this)) return true;
  return false;
}

bool CellFilter::isKnownDisjointWith(const CellFilter& other) const
{
  if (other.knownDisjoints().contains(*this)) return true;
  if (this->knownDisjoints().contains(other)) return true;

  return false;
}

bool CellFilter::isSubsetOf(const CellFilter& other,
  const Mesh& mesh) const
{
  if (isKnownSubsetOf(other)) 
  {
    return true;
  }
  else
  {
    CellSet myCells = getCells(mesh);
    CellSet yourCells = other.getCells(mesh);
    CellSet inter = myCells.setIntersection(yourCells);
    if (inter.begin() == inter.end()) return false;
    CellSet diff = myCells.setDifference(inter);
    return (diff.begin() == diff.end());
  }
}



bool CellFilter::operator==(const CellFilter& other) const
{
  if (*this < other) return false;
  if (other < *this) return false;
  return true;
}

bool CellFilter::operator!=(const CellFilter& other) const
{
  return !( *this == other );
}


const Set<CellFilter>& CellFilter::knownSubsets() const
{
  return SubsetManager::getSubsets(*this);
}
const Set<CellFilter>& CellFilter::knownDisjoints() const
{
  return SubsetManager::getDisjoints(*this);
}

void CellFilter::registerSubset(const CellFilter& sub) const
{
  SubsetManager::registerSubset(*this, sub);
  
  for (Set<CellFilter>::const_iterator 
         i=sub.knownSubsets().begin(); i!=sub.knownSubsets().end(); i++)
  {
    SubsetManager::registerSubset(*this, *i);
  }
}


void CellFilter::registerDisjoint(const CellFilter& sub) const
{
  SubsetManager::registerDisjoint(*this, sub);
  
  for (Set<CellFilter>::const_iterator 
         i=sub.knownDisjoints().begin(); i!=sub.knownDisjoints().end(); i++)
  {
    SubsetManager::registerDisjoint(*this, *i);
  }
}

void CellFilter::registerLabeledSubset(int label, const CellFilter& sub) const
{
  SubsetManager::registerLabeledSubset(*this, label, sub);
  
  const Map<int, CellFilter>& subsub = SubsetManager::getLabeledSubsets(sub);

  for (Map<int, CellFilter>::const_iterator 
         iter=subsub.begin(); iter!=subsub.end(); iter++)
  {
    if (iter->first == label) continue;
    SubsetManager::registerDisjoint(sub, iter->second);
    SubsetManager::registerDisjoint(iter->second, sub);
  }
}


XMLObject CellFilter::toXML() const 
{
  return ptr()->toXML();
}

string CellFilter::toString() const 
{
  return cfbPtr()->toString();
}

const CellFilterBase* CellFilter::cfbPtr() const
{
  const CellFilterBase* rtn = dynamic_cast<const CellFilterBase*>(ptr().get());
  TEST_FOR_EXCEPTION(rtn==0, InternalError, "CellFilter::cfbPtr() cast failed");
  return rtn;
}

CellFilterBase* CellFilter::nonConstCfbPtr()
{
  CellFilterBase* rtn = dynamic_cast<CellFilterBase*>(ptr().get());
  TEST_FOR_EXCEPTION(rtn==0, InternalError, "CellFilter::nonConstCfbPtr() cast failed");
  return rtn;
}
