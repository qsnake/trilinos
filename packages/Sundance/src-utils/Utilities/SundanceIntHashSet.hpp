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

#ifndef SUNDANCE_INTHASHSET_H
#define SUNDANCE_INTHASHSET_H

#include "SundanceDefs.hpp"
#include "SundanceSet.hpp"
#include "Teuchos_Array.hpp"
#include <list>


namespace Sundance
{
  using namespace Teuchos;

  /** 
   *
   */
  class IntHashSet
  {
  public:
    /** */
    IntHashSet();

    /** */
    void setCapacity(int capacity);

    /** */
    inline void put(int x) 
    {
      
      std::list<int>& d = data_[hashFunc(x)];
      for (std::list<int>::const_iterator i=d.begin(); i != d.end(); i++)
        {
          if (x == *i) return;
        }
      d.push_back(x);
      size_++;
    }

    /** */
    bool contains(int x) const ;

    /** */
    int size() const {return size_;}

    /** */
    void fillArray(int* a) const ;

  private:

    inline int hashFunc( const int key ) const { return (seed() ^ key)%capacity_; }

    static unsigned int seed() {static int rtn = (2654435761U); return rtn;}
    unsigned int capacity_;
    Array<std::list<int> > data_;
    int size_;
  };

}


#endif
