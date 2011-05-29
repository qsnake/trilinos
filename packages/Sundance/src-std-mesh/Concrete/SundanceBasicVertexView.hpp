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

#ifndef SUNDANCE_BASICVERTEXVIEW_H
#define SUNDANCE_BASICVERTEXVIEW_H

#include "SundanceDefs.hpp"
#include "Teuchos_Array.hpp"

namespace Sundance
{
using namespace Teuchos;

/**
 * VertexView is a read-only "view" of a cell's vertices, where the 
 * vertices are stored contiguously in a large master array. By working 
 * with views, we can greatly reduce the number of temporary arrays 
 * created during hashtable searches for existing vertex arrays.
 */
class VertexView
{
public:
  /** empty ctor, needed for storing VertexViews in Teuchos hashtables */
  VertexView() : base_(0), offset_(0), length_(0) {;}
  /** Construct a view into an array
   * \param bese pointer to the start of the master data array. By 
   * using double indirection, the master array can be resized or 
   * relocated and VertexViews can remain valid. 
   * \param  offset the index of the vertex subarray being viewed. 
   * \param length the number of vertices included in this view.
   */
  VertexView(int** base, int offset, int length)
    : base_(base), offset_(offset), length_(length) {;}

  /**
   * Return a hash code for the vertex set. 
   */
  int hashCode() const ;

  /** 
   * Test equality between two vertex sets. 
   * Two vertex sets are equal when their vertices are identical.
   */
  bool operator==(const VertexView& other) const ;

  /**
   * Write to a std::string
   */
  std::string toString() const ;


private:
  int** base_;
  int offset_;
  int length_;
};

/*
 * Two vertex sets are equal when their vertices are identical.
 * 
 */
inline bool VertexView::operator==(const VertexView& other) const
{
  /* For efficiency's sake, skip the test for equal lengths because we
   * can assume the caller is only comparing equal length vertex views */
  int* p = *base_ + offset_*length_;
  int* op = *(other.base_) + other.offset_*length_;

  for (int i=0; i<length_; i++)
  {
    if (p[i] != op[i]) return false;
  }
  return true;
}


}

namespace Teuchos
{
/** \relates VertexView */
inline int hashCode(const Sundance::VertexView& v) 
{return v.hashCode();}

/** \relates VertexView */
inline std::string toString(const Sundance::VertexView& v) 
{return v.toString();}
}

#endif
