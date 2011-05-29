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

#ifndef SUNDANCE_CELLITERATOR_H
#define SUNDANCE_CELLITERATOR_H


#include "SundanceDefs.hpp"
#include "SundanceSet.hpp"
#include "SundanceMap.hpp"
#include "SundanceCellType.hpp"
#include "SundanceCellReordererImplemBase.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "SundanceMesh.hpp"

namespace Sundance
{
using namespace Teuchos;

/**
 * CellIterator is an iterator for walking through cell sets.
 * It satisfies the requirements for an input iterator in STL.
 *
 * This class design violates the usual rules of good OO style:
 * it has polymorphic behavior packed into a single class. The
 * justification for this decision is to avoid the expense
 * of the clone() operations that would be required by the
 * copy ctor for iterators were 
 * a polymorphic class heirarchy used. 
 *
 * Two cell set types exist: explicit, where the member cells LIDs
 * are enumerated in a physical Set<int> object, and implicit,
 * where no physical Set is made, rather, the sequence of cell LIDs
 * is obtained through some scheme of walking the mesh. 
 *
 * The only cell sets that can represented implicitly are
 * the set of all cells of a given dimension.
 * 
 * \see CellSet, CellFilter, CellIteratorPos
 */
class CellIterator : public std::iterator<std::input_iterator_tag, int>
{
public:

      
  /** 
   * CellIteratorPos is used to specify whether a new CellIterator
   * is positioned at the beginning or end of a set.
   */
  enum CellIteratorPos {Begin, End};

  /** Empty ctor */
  CellIterator();

  /** Copy ctor */
  CellIterator(const CellIterator& other);

  /** Construct an implicit iterator for walking all cells of a given
   * dimension on the given mesh. */
  CellIterator(const Mesh& mesh, int cellDim, CellIteratorPos pos);

  /** Construct an explicit iterator for walking an explicitly
   * enumerated set of cells. */
  CellIterator(const Set<int>* cells, CellIteratorPos pos);

  /** */
  CellIterator& operator=(const CellIterator& other);
      
  /** Dereferencing operator */
  const int& operator*() const 
    {
      if (isImplicit_) return currentLID_;
      else return *iter_;
    }
      
  /** Postfix increment: advances iterator and returns previous value  */
  CellIterator operator++(int) 
    {
      CellIterator old = *this;
      advance();
      return old;
    }
      

  /** Prefix increment: advances iterator, returning new value */
  CellIterator& operator++()
    {
      advance();
      return *this;
    }

  /** */
  bool operator==(const CellIterator& other) const 
    {
      if (isImplicit_)
      {
        return currentLID_ == other.currentLID_;
      }
      else
      {
        return iter_ == other.iter_;
      }
    }

  /** */
  bool operator!=(const CellIterator& other) const 
    {
      return !(*this == other);
    }

      
private:

  /** Advance the iterator */
  void advance()
    {
      CellIterator old = *this;
      if (isImplicit_) 
      {
        if (reorderer_ != 0) 
        {
          currentLID_ = reorderer_->advance(currentLID_);
        }
        else currentLID_++;
      }
      else iter_++;
    }
      
  /** Flag indicating whether this iterator is implicit */
  bool isImplicit_;
      
  /** The LID to which this iterator is currently pointing.
   * Used only for implicit iterators. */
  int currentLID_;

  /** Unmanaged pointer to the reorderer used for walking 
   * implicit cell sets. Used only for implicit iterators. */
  const CellReordererImplemBase* reorderer_; 

  /** iterator for enumerated cells.
   * Used only for explicit iterators. */
  Set<int>::const_iterator iter_;

  /** */
  static Set<int> dummy() {static Set<int> rtn; return rtn;}
};
}


#endif
