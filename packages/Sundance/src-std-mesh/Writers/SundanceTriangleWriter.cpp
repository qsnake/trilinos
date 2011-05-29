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

#include "SundanceTriangleWriter.hpp"
#include "SundanceExceptions.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"


using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;


void TriangleWriter::write() const 
{
  std::string f = filename();
  if (nProc() > 1) f = f + Teuchos::toString(myRank());

  /* write header information on root proc only */
  if (nProc() > 1 && myRank()==0) writeHeader(filename());

  /* write local mesh on all procs */
  writePoints(f);
  writeCells(f);
  writeEdges(f);
  if (mesh().spatialDim() > 2)
    {
      writeFaces(f);
    }
  if (nProc() > 1) writeParallelInfo(f);
}

void TriangleWriter::writeHeader(const std::string& filename) const 
{
  std::string hdrfile = filename + ".hdr";
  std::ofstream os(hdrfile.c_str());

  os << nProc() << std::endl;
  for (int p=0; p<nProc(); p++) 
    {
      os << filename + Teuchos::toString(p) << std::endl;
    }

  os << pointScalarNames().length() << std::endl;
  for (int i=0; i<pointScalarNames().length(); i++)
    {
      os << i << " " << pointScalarNames()[i] << std::endl;
    }
  os << cellScalarNames().length() << std::endl;
  for (int i=0; i<cellScalarNames().length(); i++)
    {
      os << i << " " << cellScalarNames()[i] << std::endl;
    }
  
  for (int i=0; i<comments().length(); i++)
    {
      os << "# " << comments()[i] << std::endl;
    }
}


void TriangleWriter::writePoints(const std::string& filename) const 
{
  int nPts = mesh().numCells(0);
  int dim = mesh().spatialDim();
  int nAttr = pointScalarFields().length();
  int nBdryMarkers = 0;

  std::string nodefile = filename + ".node";
  std::ofstream os(nodefile.c_str());

  os << nPts << " " << dim << " " << nAttr << " " << nBdryMarkers << std::endl;

  for (int i=0; i<nPts; i++)
    {
      os << i+indexOffset_;
      const Point& x = mesh().nodePosition(i);
      for (int d=0; d<dim; d++)
        {
          os << " " << x[d];
        }
      /*
      for (int f=0; f<nAttr; f++)
        {
          double val = pointScalarFields()[f].probeAtMeshPoint(i);
          os << " " << val;
        }
      */
      for (int b=0; b<nBdryMarkers; b++)
        {
	  //bvbw          SUNDANCE_ERROR("Boundary markers not supported yet");
	  assert(0);
        }
      os << std::endl;
    }
  
  for (int i=0; i<comments().length(); i++)
    {
      os << "# " << comments()[i] << std::endl;
    }
}

void TriangleWriter::writeFaces(const std::string& filename) const 
{
  std::string facefile = filename + ".face";
  std::ofstream os(facefile.c_str());

  int dim = 2;
  int nFaces = mesh().numCells(dim);
  int dummySign;

  os << nFaces << " 0" << std::endl;

  for (int c=0; c<nFaces; c++)
    {
      os << c + indexOffset_;
      int nNodes = 3;
      
      for (int i=0; i<nNodes; i++)
        {
          os << " " << mesh().facetLID(2,c,0,i,dummySign) + indexOffset_;
        }
      os << std::endl;
    }
  
  for (int i=0; i<comments().length(); i++)
    {
      os << "# " << comments()[i] << std::endl;
    }
}

void TriangleWriter::writeEdges(const std::string& filename) const 
{
  std::string edgefile = filename + ".edge";
  std::ofstream os(edgefile.c_str());

  int dim = 1;
  int nEdges = mesh().numCells(dim);
  int nNodes = 2;
  int dummySign;

  os << nEdges << " 0" << std::endl;

  for (int c=0; c<nEdges; c++)
    {
      os << c + indexOffset_;
      for (int i=0; i<nNodes; i++)
        {
          os << " " << mesh().facetLID(1,c,0,i,dummySign) + indexOffset_;
        }
      os << std::endl;
    }
  
  for (int i=0; i<comments().length(); i++)
    {
      os << "# " << comments()[i] << std::endl;
    }
}

void TriangleWriter::writeCells(const std::string& filename) const 
{
  std::string elefile = filename + ".ele";
  std::ofstream os(elefile.c_str());

  int dim = mesh().spatialDim();
  int nCells = mesh().numCells(dim);
  int nAttr = cellScalarFields().length();
  int dummySign;

  os << nCells << " " << dim+1 << " " << nAttr << std::endl;

  for (int c=0; c<nCells; c++)
    {
      os << c + indexOffset_;
      int nNodes = dim+1;
      
      for (int i=0; i<nNodes; i++)
        {
          os << " " << mesh().facetLID(dim,c,0,i,dummySign) + indexOffset_;
        }
      /*
      for (int f=0; f<nAttr; f++)
        {
          os << " " << cellScalarFields()[f].average(cell).value();
        }
      */
      os << std::endl;
    }
  
  for (int i=0; i<comments().length(); i++)
    {
      os << "# " << comments()[i] << std::endl;
    }
}


void TriangleWriter::writeParallelInfo(const std::string& filename) const 
{
  std::string parfile = filename + ".par";
  std::ofstream os(parfile.c_str());

  int dim = mesh().spatialDim();
  int nCells = mesh().numCells(dim);
  int nEdges = mesh().numCells(1);
  int nPts = mesh().numCells(0);

  os << myRank() << " " << nProc() << std::endl;

  os << nPts << std::endl;
  for (int i=0; i<nPts; i++)
    {
      os << i << " " << mesh().mapLIDToGID(0,i) 
         << " " << mesh().ownerProcID(0,i) << std::endl;
    }

  os << nCells << std::endl;
  for (int c=0; c<nCells; c++)
    {
      os << c << " " << mesh().mapLIDToGID(dim,c) 
         << " " << mesh().ownerProcID(dim,c) << std::endl;
    }

  os << nEdges << std::endl;
  for (int c=0; c<nEdges; c++)
    {
      os << c << " " << mesh().mapLIDToGID(1,c) 
         << " " << mesh().ownerProcID(1,c) << std::endl;
    }

  if (dim > 2)
    {
      int nFaces = mesh().numCells(2);
      os << nFaces << std::endl;
      for (int c=0; c<nFaces; c++)
        {
          os << c << " " << mesh().mapLIDToGID(2,c) 
             << " " << mesh().ownerProcID(2,c) << std::endl;
        }
    }
  
  for (int i=0; i<comments().length(); i++)
    {
      os << "# " << comments()[i] << std::endl;
    }
}


