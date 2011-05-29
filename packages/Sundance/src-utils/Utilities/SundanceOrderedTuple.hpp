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

#ifndef SUNDANCE_ORDEREDTUPLE_H
#define SUNDANCE_ORDEREDTUPLE_H

#include "SundanceDefs.hpp"


#ifndef DOXYGEN_DEVELOPER_ONLY

namespace Sundance
{
  using namespace Teuchos;

  /** OrderedPair provides a means of lexigraphic comparison of a pair of
   * objects. The pair {a1, b1} is compared to {a2, b2} by first
   * comparing the most significant entries a1 and a2, and if they are
   * equal, comparing the least significant entries b1 and b2. */
  template<class A, class B>
    class OrderedPair
    {
    public:
      /** */
      OrderedPair(const A& _a, const B& _b)
        : a_(_a), b_(_b) {;}

      /** */
      inline bool operator<(const OrderedPair<A, B>& other) const
        {
          if ( a_ < other.a_ ) 
            {
              return true;
            }
          if ( other.a_ < a_) 
            {
              return false;
            }

          bool rtn = b_ < other.b_;
          return rtn;
        }

      /** */
      const A& first() const {return a_;}

      /** */
      const B& second() const {return b_;}

    private:
      A a_;
      B b_;
    };

  /** Lexigraphically-comparable triple of objects. */
  template<class A, class B, class C>
    class OrderedTriple : public OrderedPair<A, OrderedPair<B, C> >
    {
    public:
      /** */
      OrderedTriple(const A& _a, const B& _b, const C& _c)
        : OrderedPair<A, OrderedPair<B, C> >(_a, OrderedPair<B,C>(_b,_c))
        {;}

      const A& a() const {return this->first();}

      const B& b() const {return this->second().first();}

      const C& c() const {return this->second().second();}
    };

  /** Lexigraphically-comparable quartet of objects. */
  template<class A, class B, class C, class D>
    class OrderedQuartet : public OrderedPair<A, OrderedTriple<B, C, D> >
    {
    public:
      /** */
      OrderedQuartet(const A& _a, const B& _b, const C& _c, const D& _d)
        : OrderedPair<A, OrderedTriple<B, C, D> >(_a, OrderedTriple<B,C,D>(_b,_c,_d))
        {;}

      const A& a() const {return this->first();}
      const B& b() const {return this->second().first();}
      const C& c() const {return this->second().second().first();}
      const D& d() const {return this->second().second().second();}
    };

  /** */
  template <class A, class B>
  inline std::ostream& operator<<(std::ostream& os, const OrderedPair<A,B>& p)
  {
    os << "{" << p.first() << ", " << p.second() << "}";
    return os;
  }

}


#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
