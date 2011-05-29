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

#include "SundanceMultiIndex.hpp"
#include "SundanceExceptions.hpp"

using namespace Sundance;
using namespace Sundance;

using namespace Sundance;
using namespace Teuchos;

MultiIndex::MultiIndex()
	: m_(maxDim(), 0)
{;}

MultiIndex::MultiIndex(int x, int y, int z)
	: m_(maxDim(), 0)
{
	m_[0] = x;
	m_[1] = y;
	m_[2] = z;
}

MultiIndex MultiIndex::operator+(const MultiIndex& other) const 
{
	MultiIndex rtn;

	for (int i=0; i<maxDim(); i++)
		{
			rtn.m_[i] = m_[i] + other[i];
		}
	return rtn;
}

MultiIndex MultiIndex::operator-(const MultiIndex& other) const 
{
	MultiIndex rtn;

	for (int i=0; i<maxDim(); i++)
		{
			rtn.m_[i] = m_[i] - other[i];
		}
	return rtn;
}

MultiIndex MultiIndex::operator-() const 
{
	MultiIndex rtn;

	for (int i=0; i<maxDim(); i++)
		{
			rtn.m_[i] = -m_[i];
		}
	return rtn;
}

string MultiIndex::toString() const
{
	return "(" + Teuchos::toString(m_[0]) + ","
		+ Teuchos::toString(m_[1]) + ","
		+ Teuchos::toString(m_[2]) + ")";
} 

XMLObject MultiIndex::toXML() const
{
	XMLObject rtn("MultiIndex");
	rtn.addAttribute("indices", toString());
	return rtn;
}

bool MultiIndex::operator==(const MultiIndex& m) const 
{
	for (int i=0; i<maxDim(); i++)
		{
			if (m_[i] != m[i]) return false;
		}
	return true;
}

bool MultiIndex::operator<(const MultiIndex& m) const 
{
	for (int i=0; i<maxDim(); i++)
		{
			if (m_[i] > m.m_[i]) return false;
      if (m_[i] < m.m_[i]) return true;
		}
	return false;
}



int MultiIndex::order() const 
{
  int h = 0;
	for (int i=0; i<maxDim(); i++)
		{
			h += m_[i];
		}
	return h;
}

bool MultiIndex::isValid() const 
{
	for (int i=0; i<maxDim(); i++)
		{
      if (m_[i] < 0) return false;
		}
	return true;
}

int MultiIndex::firstOrderDirection() const 
{
  TEST_FOR_EXCEPTION(order() != 1, InternalError,
                     "bad order in MultiIndex::firstOrderDirection() const");
	for (int i=0; i<maxDim(); i++)
		{
			if (m_[i] == 1) return i;
		}
  return -1;
}



string MultiIndex::coordForm() const
{
  std::string rtn;
  
  for (int i=0; i<m_[0]; i++)
		{
      rtn += "x";
		}
  for (int i=0; i<m_[1]; i++)
		{
      rtn += "y";
		}
  for (int i=0; i<m_[2]; i++)
		{
      rtn += "z";
		}
	return rtn;
}


