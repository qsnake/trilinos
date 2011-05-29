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

#ifndef SUNDANCE_COMBINATORIALUTILS_H
#define SUNDANCE_COMBINATORIALUTILS_H

#include "SundanceDefs.hpp"
#include "Teuchos_Array.hpp"
#include "SundanceMultiSet.hpp"
#include "SundanceMap.hpp"
#include "SundanceOrderedTuple.hpp"

namespace Sundance
{
using Teuchos::Array;
/**
 * Return partitions of an integer
 * @author Kevin Long
 */
Array<Array<int> > partitionInteger(int n);

/** 
 * Return compositions of an integer
 */
Array<Array<Array<int> > > compositions(int n);


/** 
 * Return the non-negative compositions of an integer into J terms
 */
Array<Array<int> > nonNegCompositions(int n, int J);



/** 
 * Return the s-term compositions of a multiset. These are the permutations
 * of the s-term partitions. 
 */
Array<Array<MultiSet<int> > > multisetCompositions(int s,
  const MultiSet<int>& x);



/** Return all subsets of a multiset. */
Set<MultiSet<int> > multisetSubsets(const MultiSet<int>& x);

/** Return all N-tuples of subsets of a multisets such that the
 * union of tuple entries is a subset of the original multiset. 
 * For example, if the input multiset is {a,b,b}, the subsets are
 * \pre
 * {a}, {b}, {a,b}, {a, b, b}, {b, b}.
 * \pre
 * The 2-tuples are all pairs drawn from that collection. Of these, only
 * \pre
 * ({a}, {b})
 * ({a}, {b, b})
 * \pre 
 * satisfy the requirement that the union of tuple entries is a subset
 * of {a,b,b}.
 */
Array<Array<MultiSet<int> > > multisetSubsetNTuples(int n, 
  const MultiSet<int>& x);


/** 
 * Generate the (n-choose-m) distinct index combinations for
 * choosing m entries from an array of n elements.
 */
Array<Array<int> > distinctIndexTuples(int m, int n);

/** 
 * For a given integer vector of length N, find all combinations 
 * of integers (i_1, 1_2, ... i_N) such that
 * \f$0 \le i_j < s_j\f$. 
 * For example, if \f$ s = (2,3) \f$, return
 * \code
 * { {0,0}, {0,1}, {0,2}, {1,0}, {1,1}, {1,2} }
 * \endcode
 */
Array<Array<int> > indexCombinations(const Array<int>& s);

/** Compute the n-th power of 2 */
int pow2(int n);

/** Compute the n bits of an integer x < 2^n. The bits are ordered
 * least-to-most significant.
 */
Array<int> bitsOfAnInteger(int x, int n);


/** 
 *
 */
template <class T> 
class Pair: public OrderedPair<T, T>
{
public:
  Pair(const T& a, const T& b)
    : OrderedPair<T, T>(a, b)
    {;}
};

template <class T> inline
std::ostream& operator<<(std::ostream& os, const Pair<T>& p)
{
  os << "Pair(" << p.first() << ", " << p.second() << ")";
  return os;
}

/** 
 *
 */
template <class T> 
class SortedPair: public OrderedPair<T, T>
{
public:
  SortedPair(const T& a, const T& b)
    : OrderedPair<T, T>(std::min(a, b), std::max(a,b))
    {;}
};

template <class T> inline
std::ostream& operator<<(std::ostream& os, const SortedPair<T>& p)
{
  os << "Pair(" << p.first() << ", " << p.second() << ")";
  return os;
}


/** 
 * Form all pairs of multisets that can be generated from an 
 * initial pair by adding n copies of entry x. This is used in the
 * partitioning of multisets.
 */
Set<Pair<MultiSet<int> > >
loadPartitions(int x, int n, 
  const MultiSet<int>& left, 
  const MultiSet<int>& right);


/** 
 * Return the size-2 partitions of a multiset. This is used in a recursive
 * algorithm to determine all partitions of a multset. 
 */
Set<Pair<MultiSet<int> > >
binaryPartition(const MultiSet<int>& m);

/** 
 * Return the partitions of a multiset. The partitions are sets
 * of subsets such that the union of the members equals the original
 * multiset. 
 */
Set<MultiSet<MultiSet<int> > >
multisetPartitions(const MultiSet<int>& m);

/** 
 * Given a multiset, create a mapping from entry to number of repetitions
 * in the multisets.
 */
Map<int, int> countMap(const MultiSet<int>& m);


/** */
template <class T> inline
Array<Array<Array<T> > > indexArrangements(const MultiSet<T>& mu,
  const Array<int>& k)
{
  int nBins = k.size();
    
  int M = 0;
  for (int i=0; i<nBins; i++)
  {
    M += k[i];
  }
    
  Array<T> I;
  typename MultiSet<T>::const_iterator iter;
  for (iter=mu.begin(); iter!=mu.end(); iter++)
  {
    I.append(*iter);
  }

  Array<Array<Array<T> > > rtn;
    
  do
  {
    Array<Array<T> > bins(nBins);
    int count = 0;
    for (int i=0; i<nBins; i++)
    {
      for (int j=0; j<k[i]; j++)
      {
        bins[i].append(I[count++]);
      }
    }
    rtn.append(bins);
  }
  while (std::next_permutation(I.begin(), I.end()));
  return rtn;
}

/** 
 * Return the distinct arrangements of the given multiset into n bins
 */
template <class T> inline
Array<Array<Array<T> > > binnings(const MultiSet<T>& mu, int n)
{
  int N = mu.size();
  Array<Array<int> > c = compositions(N)[n-1];
  Array<Array<Array<T> > > rtn;

  for (int i=0; i<c.size(); i++)
  {
    Array<Array<Array<T> > > a = indexArrangements(mu, c[i]);
    for (int j=0; j<a.size(); j++)
    {
      rtn.append(a[j]);
    }
  }
  return rtn;
}


inline int factorial(const MultiSet<int>& ms)
{
  Sundance::Map<int, int> counts;
    
  for (MultiSet<int>::const_iterator i=ms.begin(); i!=ms.end(); i++)
  {
    if (counts.containsKey(*i)) counts[*i]++;
    else counts.put(*i, 1);
  }

  int rtn = 1;
  for (Sundance::Map<int, int>::const_iterator
         i=counts.begin(); i!=counts.end(); i++)
  {
    int f = 1;
    for (int j=1; j<=i->second; j++) f *= j;
    rtn *= f;
  }
  return rtn;
}


}



#endif



