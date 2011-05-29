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

#include "SundanceVerboseFieldWriter.hpp"
#include "SundanceExceptions.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"


using namespace Sundance;
using namespace Teuchos;


void VerboseFieldWriter::write() const 
{
  int nProc = mesh().comm().getNProc();
  int myRank = mesh().comm().getRank();

  RCP<std::ostream> osp;
  if (filename().length()==0)
    {
      osp = rcp(&std::cout, false);
    }
  else 
    {
      std::string f = filename() + ".txt";
      if (nProc > 1) f = f + "." + Teuchos::toString(myRank);
      osp = rcp(new std::ofstream(f.c_str()));
    }
  std::ostream& os = *osp;

  if (myRank==0) os << "VerboseFieldWriter output" << std::endl;
  for (int p=0; p<nProc; p++)
    {
      mesh().comm().synchronize();
      mesh().comm().synchronize();
      mesh().comm().synchronize();
      if (p != myRank) continue;
      os << "======== processor " << p << " ============================ "
         << std::endl;
      Tabs tab0;
      int dim = mesh().spatialDim();
      int nPts = mesh().numCells(0);
      int nElems = mesh().numCells(dim);
      os << tab0 << "spatial dimension = " << dim << std::endl;
      os << tab0 << "num points = " << nPts << std::endl;
      os << tab0 << "num elements = " << nElems << std::endl;
      os << tab0 << "Point list: " << std::endl;

      int dummy;
      for (int i=0; i<nPts; i++)
        {
          Tabs tab1;
          os << tab1 << "L=" << i 
             << " G=" << mesh().mapLIDToGID(0, i) 
             << " x=" << mesh().nodePosition(i) 
             << " owner=" << mesh().ownerProcID(0,i) 
             << " label=" << mesh().label(0,i) << std::endl;
          int nc = mesh().numMaxCofacets(0,i);
          Tabs tab2;
          os << tab2 << "num cofacets=" << nc << " cofs = {";
          for (int c=0; c<nc; c++)
            {
              if (c==0) os << " " ;
              else os << ", ";
              os << mesh().mapLIDToGID(dim, mesh().maxCofacetLID(0,i,c,dummy));
            }
          os << "}" << std::endl;
        }

      
      os << tab0 << "Element list: " << std::endl;

      for (int i=0; i<nElems; i++)
        {
          int facetSign;
          Tabs tab1;
          os << tab1 << "L=" << i 
             << " G=" << mesh().mapLIDToGID(dim, i) 
             << ", nodes L={";
          int numNodes = mesh().numFacets(dim, i, 0);
          for (int j=0; j<numNodes; j++)
            {
              if (j != 0) os << ", ";
              os << mesh().facetLID(dim, i, 0, j, facetSign);
            }
          os << "}, G={";
          for (int j=0; j<numNodes; j++)
            {
              if (j != 0) os << ", ";
              os << mesh().mapLIDToGID(0, mesh().facetLID(dim, i, 0, j, facetSign));
            }
          os << "}, owner=" << mesh().ownerProcID(dim,i)
             << ", label=" << mesh().label(dim,i) << std::endl;
          for (int fd=1; fd<dim; fd++)
            {
              Tabs tab2;
              os << tab2 << "facets of dimension " << fd << std::endl;
              int nf = mesh().numFacets(dim, i, fd);
              for (int f=0; f<nf; f++)
                {

                  Tabs tab3;
                  int flid = mesh().facetLID(dim, i, fd, f, facetSign);
                  int fgid = -1;
                  int fowner = -1;
                  if (mesh().hasIntermediateGIDs(fd))
                    { 
                      fgid = mesh().mapLIDToGID(fd, flid);
                      fowner = mesh().ownerProcID(fd, flid);
                    }
                  os << tab3 << "f#=" << f << " L=" << flid
                     << " G=" << fgid << " owner=" << fowner
                     << " nodes={";
                  int nfn = mesh().numFacets(fd, flid, 0);
                  for (int fn=0; fn<nfn; fn++)
                    {
                      if (fn != 0) os << ", ";
                      os << mesh().facetLID(fd, flid, 0, fn, facetSign);
                    }
                  os << "} sign=" << facetSign << std::endl;
                }
              
            }
        }
      
    }
}



