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

#include "Sundance.hpp"
#include "SundanceTriangleQuadrature.hpp"
#include "SundanceFeketeTriangleQuadrature.hpp"
#include "SundanceTetQuadrature.hpp"

using namespace Sundance;
using namespace Teuchos;
using namespace TSFExtended;


int main(int argc, char** argv)
{
  int stat = 0 ;
	try
		{
			GlobalMPISession session(&argc, &argv);

      Array<int> validTetOrders = tuple(1, 2, 4, 6);
      int maxorder = 15;
      Array<int> validTriOrders(maxorder);
      for (int i=0; i<maxorder; i++)
    	  validTriOrders[i] = i+1;

      Array<int> validFeketeTriOrders = tuple(1, 2, 3, 4, 5, 6, 9);


      Array<int> triFailures;
      Array<int> tetFailures;
      Array<int> FeketeTriFailures;

      std::cerr << "------------- testing triangle rules -------------------"  << std::endl;
      for (int i=0; i<validTriOrders.size(); i++)
				{
          int p = validTriOrders[i];
					bool pass = TriangleQuadrature::test(p);
					if (pass) std::cerr << "order " << p << " PASSED" << std::endl;
					else 
            {
              std::cerr << "order " << p << " FAILED" << std::endl;
              triFailures.append(p);
            }
				}
      std::cerr << "------------- testing tet rules -------------------"  << std::endl;
          
      for (int i=0; i<validTetOrders.size(); i++)
				{
          int p = validTetOrders[i];
					bool pass = TetQuadrature::test(p);
					if (pass) std::cerr << "order " << p << " PASSED" << std::endl;
					else 
            {
              std::cerr << "order " << p << " FAILED" << std::endl;
              tetFailures.append(p);
            }
				}

      std::cerr << "--------- testing Fekete triangle rules ----------------" << std::endl;
      for (int i = 0; i < validFeketeTriOrders.size(); i++)
		{
			int p = validFeketeTriOrders[i];
			bool pass = FeketeTriangleQuadrature::test(p);
			if (pass)
				cerr << "order " << p << " PASSED" << std::endl;
			else
			{
				cerr << "order " << p << " FAILED" << std::endl;
				FeketeTriFailures.append(p);
			}
		}


      if (tetFailures.size()>0) 
        {
          cout << "failures detected for tets: orders " << tetFailures << std::endl;
          cout << "tet tests FAILED" << std::endl;
          stat = -1;
        }
      else
        {
          cout << "tet tests PASSED" << std::endl;
        }

      if (triFailures.size()>0) 
        {
          cout << "failures detected for tris: orders " << triFailures << std::endl;
          cout << "tri tests FAILED" << std::endl;
          stat = -1;
        }
      else
        {
          cout << "tri tests PASSED" << std::endl;
        }

		if (FeketeTriFailures.size() > 0)
		{
			cout << "failures detected for Fekete tris: orders "
					<< FeketeTriFailures << std::endl;
			cout << "Fekete tri tests FAILED" << std::endl;
			stat = -1;
		}
		else
		{
			cout << "Fekete tri tests PASSED" << std::endl;
		}

		}
	catch(std::exception& e)
		{
      std::cerr << "Detected exception: " << e.what() << std::endl;
      std::cerr << "Quadrature test FAILED" << std::endl;
      stat = -1;
		}

  return stat;
  
}
