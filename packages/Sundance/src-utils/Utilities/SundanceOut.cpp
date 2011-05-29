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

#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"
#include "Teuchos_Array.hpp"

using namespace Sundance;
using namespace Teuchos;
using namespace std;


namespace Sundance
{

void writeTable(std::ostream& os, const Tabs& tab,
  const Array<double>& a, int cols)
{
  int rows = a.size() / cols;

  for (int i=0; i<rows; i++)
  {
    os << tab << setw(10) << i << ":";
    for (int j=0; j<cols; j++) 
      os << setw(12) << setprecision(6) << a[i*cols+j];
    os << std::endl;
  }
  int n = a.size() - rows * cols;
  if (n==0) return ;
  os << tab << setw(10) << rows << ":" ;
  for (int j=0; j<n; j++) 
    os << setw(12) << setprecision(6) << a[rows*cols+j];
  os << std::endl;
}


void writeTable(std::ostream& os, const Tabs& tab,
  const Array<int>& a, int cols)
{
  int rows = a.size() / cols;

  for (int i=0; i<rows; i++)
  {
    os << tab << setw(10) << i << ":";
    for (int j=0; j<cols; j++) 
      os << setw(10) << a[i*cols+j];
    os << std::endl;
  }
  int n = a.size() - rows * cols;
  if (n==0) return ;
  os << tab << setw(10) << rows << ":" ;
  for (int j=0; j<n; j++) 
    os << setw(10) << a[rows*cols+j];
  os << std::endl;
}


}







