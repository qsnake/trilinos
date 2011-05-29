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

#include "SundanceVTKWriter.hpp"
#include "SundanceExceptions.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"
#include "Teuchos_XMLObject.hpp"



using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;




void VTKWriter::write() const 
{
  lowLevelWrite(filename(), false);
  if (nProc() > 1 && myRank()==0) lowLevelWrite(filename(), true);
}

void VTKWriter::lowLevelWrite(const std::string& filename, bool isPHeader) const 
{
  std::string PHeader = "";
  if (isPHeader) PHeader="P";

  std::string f = filename;
  
  if (isPHeader) f = f + ".pvtu";
  else if (nProc() > 1) 
    {
      f = f + Teuchos::toString(myRank()) + ".vtu";
    }
  else
    {
      f = f + ".vtu";
    }
  
  SUNDANCE_VERB_MEDIUM("writing VTK file " << f);

  std::ofstream os(f.c_str());

  XMLObject head("VTKFile");
  head.addAttribute("type", PHeader + "UnstructuredGrid");
  head.addAttribute("version", "0.1");
  
  os << head.header() << std::endl;

	for (int i=0; i<comments().length(); i++)
		{
			os << "<!-- " << comments()[i] << " -->" << std::endl;
		}

  XMLObject ug(PHeader + "UnstructuredGrid");
  os << ug.header() << std::endl;

  if (isPHeader)
    {
      writePoints(os, isPHeader);
      writePointData(os, isPHeader);
      writeCellData(os, isPHeader);
      for (int p=0; p<nProc(); p++)
        {
          XMLObject pc("Piece");
          std::string pfile = filename + Teuchos::toString(p) + ".vtu";
          size_t pos = pfile.find_last_of("/");
          if (pos != std::string::npos)
          {
            pfile = pfile.substr(pos+1);
          }
          pc.addAttribute("Source", pfile);
          os << pc << std::endl;
        }
    }
  else
    {
      XMLObject pc("Piece");
      pc.addAttribute("NumberOfPoints", Teuchos::toString(mesh().numCells(0)));
      pc.addAttribute("NumberOfCells", Teuchos::toString(mesh().numCells(mesh().spatialDim())));

      os << pc.header() << std::endl;

      writePoints(os, false);
      writeCells(os);
      writePointData(os, false);
      writeCellData(os, false);

      os << pc.footer() << std::endl;
    }

	os << ug.footer() << std::endl;
	os << head.footer() << std::endl;
}

void VTKWriter::writePoints(std::ostream& os, bool isPHeader) const 
{
  std::string PHeader = "";
  if (isPHeader) PHeader="P";
  XMLObject pts(PHeader + "Points");

  os << pts.header() << std::endl;

  XMLObject xml(PHeader + "DataArray");
  xml.addAttribute("NumberOfComponents", "3");
  xml.addAttribute("type", "Float32");
  xml.addAttribute("format", "ascii");

  os << xml.header() << std::endl;

  /* write the points, unless this call is for the dummy header on the root proc */
  if (!isPHeader)
    {
      int np = mesh().numCells(0);
      int dim = mesh().spatialDim();
      
      for (int i=0; i<np; i++)
        {
          const Point& x = mesh().nodePosition(i);
          
          for (int d=0; d<dim; d++)
            {
              os << x[d] << " ";
            }
          for (int d=dim; d<3; d++)
            {
              os << "0.0 ";
            }
          os << std::endl;
        }
    }

  os << xml.footer() << std::endl;

  os << pts.footer() << std::endl;
}


void VTKWriter::writeCells(std::ostream& os) const 
{
  XMLObject cells("Cells");
  os << cells.header() << std::endl;

  XMLObject conn("DataArray");
  conn.addAttribute("type", "Int32");
  conn.addAttribute("Name", "connectivity");
  conn.addAttribute("format", "ascii");

  int dim = mesh().spatialDim();
  int nc = mesh().numCells(dim);
  int dummySign;

  os << conn.header() << std::endl;
  CellType cellType = mesh().cellType(dim);
  
  for (int c=0; c<nc; c++)
    {
      int nNodes = mesh().numFacets( dim , c , 0 );

		switch(cellType)
			{
			case TriangleCell: case LineCell: case TetCell:
			      for (int i=0; i<nNodes; i++)
			        {
			          os << " " << mesh().facetLID(dim,c,0,i,dummySign);
			        }
			      os << std::endl;
				break;
			case QuadCell:  // for quads we have different order
		          os << " " << mesh().facetLID(dim,c,0,0,dummySign);
		          os << " " << mesh().facetLID(dim,c,0,1,dummySign);
		          os << " " << mesh().facetLID(dim,c,0,3,dummySign);
		          os << " " << mesh().facetLID(dim,c,0,2,dummySign);
		          os << std::endl;
				break;
			case BrickCell:  // for quads we have different order
		          os << " " << mesh().facetLID(dim,c,0,0,dummySign);
		          os << " " << mesh().facetLID(dim,c,0,1,dummySign);
		          os << " " << mesh().facetLID(dim,c,0,2,dummySign);
		          os << " " << mesh().facetLID(dim,c,0,3,dummySign);
		          os << " " << mesh().facetLID(dim,c,0,4,dummySign);
		          os << " " << mesh().facetLID(dim,c,0,5,dummySign);
		          os << " " << mesh().facetLID(dim,c,0,6,dummySign);
		          os << " " << mesh().facetLID(dim,c,0,7,dummySign);
		          os << std::endl;
				break;
            default:
            	 TEST_FOR_EXCEPTION(true, RuntimeError, "call type " << cellType <<
            			 " not handled in VTKWriter::writeCells()");
			}
    }
  
  os << conn.footer() << std::endl;


  XMLObject offsets("DataArray");
  offsets.addAttribute("type", "Int32");
  offsets.addAttribute("Name", "offsets");
  offsets.addAttribute("format", "ascii");
  
  os << offsets.header() << std::endl;

  int count = 0;
  for (int c=0; c<nc; c++)
    {
			count += mesh().numFacets(dim, c, 0);
      os << count << std::endl;
    }

  os << offsets.footer() << std::endl;

  XMLObject types("DataArray");
  types.addAttribute("type", "UInt8");
  types.addAttribute("Name", "types");
  types.addAttribute("format", "ascii");

  os << types.header() << std::endl;

  for (int c=0; c<nc; c++)
    {

			int vtkCode = 0;
			switch(cellType)
				{
				case TriangleCell:
					vtkCode = 5;
					break;
				case QuadCell:
					vtkCode = 9;
					break;
				case TetCell:
					vtkCode = 10;
					break;
				case BrickCell:
					vtkCode = 11;
					break;
        default:
          TEST_FOR_EXCEPTION(true, RuntimeError,
                             "call type " << cellType << " not handled "
                             "in VTKWriter::writeCells()");
				}
			os << vtkCode << std::endl;
    }

  os << types.footer() << std::endl;

  os << cells.footer() << std::endl;
}

void VTKWriter::writePointData(std::ostream& os, bool isPHeader) const 
{
  std::string PHeader = "";
  if (isPHeader) PHeader="P";

  XMLObject xml(PHeader + "PointData");

  if (pointVectorNames().length() > 0) xml.addAttribute("Vectors", pointVectorNames()[0]);
  if (pointScalarNames().length() > 0) xml.addAttribute("Scalars", pointScalarNames()[0]);

  os << xml.header() << std::endl;

  for (int i=0; i<pointScalarNames().length(); i++)
    {
      writeDataArray(os, pointScalarNames()[i], pointScalarFields()[i], isPHeader, true);
    }

  for (int i=0; i<pointVectorNames().length(); i++)
    {
      writeDataArray(os, pointVectorNames()[i], pointVectorFields()[i], isPHeader, true);
    }

  os << xml.footer() << std::endl;
}

void VTKWriter::writeCellData(std::ostream& os, bool isPHeader) const 
{
  std::string PHeader = "";
  if (isPHeader) PHeader="P";

  XMLObject xml(PHeader + "CellData");

  if (cellVectorNames().length() > 0) xml.addAttribute("Vectors", cellVectorNames()[0]);
  if (cellScalarNames().length() > 0) xml.addAttribute("Scalars", cellScalarNames()[0]);

  os << xml.header() << std::endl;

  for (int i=0; i<cellScalarNames().length(); i++)
    {
      writeDataArray(os, cellScalarNames()[i], cellScalarFields()[i], isPHeader, false);
    }

  /* ---- vector plot works intermeadaty for VTK's */
  for (int i=0; i<cellVectorNames().length(); i++)
    {
      writeDataArray(os, cellVectorNames()[i], cellVectorFields()[i], isPHeader, false);
    }

  os << xml.footer() << std::endl;
}


void VTKWriter::writeDataArray(std::ostream& os, const std::string& name, 
                               const RCP<FieldBase>& expr, bool isPHeader, bool isPointData) const 
{
  std::string PHeader = "";
  if (isPHeader) PHeader="P";

  XMLObject xml(PHeader + "DataArray");
  xml.addAttribute("type", "Float32");
  xml.addAttribute("Name", name);
  xml.addAttribute("format", "ascii");
  
  if (expr->numElems() > 1)
    {
	  // Since we are plotting always in 3D the vector components have 3 component and not "expr->numElems()"
      xml.addAttribute("NumberOfComponents", "3" );
    }
  
  os << xml.header() << std::endl;

  /* write the point|cell data, unless this is a parallel header */
  if (!isPHeader)
    {

      if (isPointData)
        {
          int numNodes = mesh().numCells(0);

          for (int i=0; i<numNodes; i++)
            {
              for (int j=0; j < expr->numElems(); j++)
                {
                  if (expr->isDefined(0,i,j)){
                	double val = expr->getData(0, i, j);
                    val = (fabs(val) > 1e-16) ? val : 0.0;
                    os << (float) val << std::endl;
                  }else
                    os << undefinedValue() << std::endl;
                }
              // write the rest with 0.0 if it is a zero component
              if (expr->numElems() > 1)
            	  for (int j= expr->numElems(); j < 3 ; j++)
            		  os << "0.0 " << std::endl;
            }
        }
      else
        {
          int dim = mesh().spatialDim();
          int nc = mesh().numCells(dim);
          
          for (int c=0; c<nc; c++)
            {
              for (int j=0; j < expr->numElems(); j++)
                {
                  if (expr->isDefined(dim,c,j)){
                   	double val = expr->getData(dim, c, j);
                    val = (fabs(val) > 1e-16) ? val : 0.0;
                    os << (float) val << std::endl;
                  }else
                    os << undefinedValue() << std::endl;
                }
              // write the rest with 0.0 if it is a zero component
              if (expr->numElems() > 1)
            	  for (int j= expr->numElems(); j < 3 ; j++)
            		  os << "0.0 " << std::endl;
            }
        }
    }

  os << xml.footer() << std::endl;
}

