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

#include "SundanceSpatiallyHomogeneousDOFMapBase.hpp"
#include "SundanceTabs.hpp"
#include "SundanceOut.hpp"
#include "SundanceMaximalCellFilter.hpp"

using namespace Sundance;
using namespace Teuchos;


SpatiallyHomogeneousDOFMapBase
::SpatiallyHomogeneousDOFMapBase(const Mesh& mesh,
  int nTotalFuncs,
  int setupVerb)
  : DOFMapBase(mesh, setupVerb), allowedFuncs_(), funcDomains_()
{
  RCP<Set<int> > f = rcp(new Set<int>());
  for (int i=0; i<nTotalFuncs; i++) f->put(i);
  allowedFuncs_ = f;
  CellFilter cf = new MaximalCellFilter();
  funcDomains_ = Array<CellFilter>(nTotalFuncs, cf);
}


void SpatiallyHomogeneousDOFMapBase::print(std::ostream& os) const
{
  int myRank = mesh().comm().getRank();

  Tabs tabs;
  int dim = mesh().spatialDim();
  RCP<const MapStructure> s = mapStruct();

  for (int p=0; p<mesh().comm().getNProc(); p++)
  {
    mesh().comm().synchronize();
    mesh().comm().synchronize();
    if (p == myRank)
    {
      os << tabs << 
        "========= DOFMap on proc p=" << p << " =============" << std::endl;
      for (int d=dim; d>=0; d--)
      {
        Tabs tabs1;
        os << tabs1 << "dimension = " << d << std::endl;
        for (int c=0; c<mesh().numCells(d); c++)
        {
          Tabs tabs2;
          os << tabs2 << "Cell d=" << d << " LID=" << c << " GID=" 
             << mesh().mapLIDToGID(d, c);
          if (d==0) 
          {
            os << " x=" << mesh().nodePosition(c) << std::endl;
          }
          else 
          {
            Array<int> facetLIDs;
            Array<int> facetDirs;
            mesh().getFacetArray(d, c, 0, facetLIDs, facetDirs);
            Array<int> facetGIDs(facetLIDs.size());
            for (int v=0; v<facetLIDs.size(); v++)
            {
              facetGIDs[v] = mesh().mapLIDToGID(0, facetLIDs[v]);
            }
            os << " nodes LIDs=" << facetLIDs << " GIDs=" << facetGIDs
               << std::endl;
          }
          for (int b=0; b<s->numBasisChunks(); b++)
          {
            for (int f=0; f<s->funcs(b).size(); f++)
            {
              Tabs tabs3;
              Array<int> dofs;
              getDOFsForCell(d, c, s->funcs(b)[f], dofs);
              os << tabs3 << "f=" << s->funcs(b)[f] << " " 
                 << dofs << std::endl;
              if (false)
              {
                os << tabs3 << "{";
                for (int i=0; i<dofs.size(); i++)
                {
                  if (i != 0) os << ", ";
                  if (isLocalDOF(dofs[i])) os << "L";
                  else os << "R";
                }
                os << "}" << std::endl;
              }
            }
          }
        }
      }
    }
    mesh().comm().synchronize();
    mesh().comm().synchronize();
  }
}
