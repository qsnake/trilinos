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

#include "SundanceFieldWriterBase.hpp"
#include "SundanceExceptions.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"


using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;


FieldWriterBase::FieldWriterBase(const std::string& filename) 
  : filename_(filename),
    mesh_(),
    nProc_(0), 
    myRank_(-1),
    meshID_(-1),
    comments_(),
    pointScalarFields_(),
    cellScalarFields_(),
    pointVectorFields_(),
    cellVectorFields_(),
    pointScalarNames_(),
    cellScalarNames_(),
    pointVectorNames_(),
    cellVectorNames_(),
    undefinedValue_(0.0)
{;}


void FieldWriterBase::impersonateParallelProc(int nProc, int rank)
{
  nProc_ = nProc;
  myRank_ = rank;
}

int FieldWriterBase::nProc() const
{
  if (nProc_ < 1) return mesh().comm().getNProc(); 
  return nProc_;
}

int FieldWriterBase::myRank() const
{
  if (myRank_ < 0) return mesh().comm().getRank(); 
  return myRank_;
}




void FieldWriterBase::addMesh(const Mesh& mesh) 
{
  if (meshID_ < 0)
    {
      mesh_ = mesh;
      meshID_ = mesh.id();
    }
                     
  TEST_FOR_EXCEPTION(meshID_ != mesh.id(), RuntimeError,
                     "FieldWriterBase::setMesh(): inconsistent meshes: "
                     "existing mesh has meshID=" << meshID_ << ", newly "
                     "added mesh has meshID=" << mesh.id());
}

void FieldWriterBase::addField(const std::string& name, 
                               const RCP<FieldBase>& expr) 
{

  std::string fieldName = name;

  if (expr->numElems() > 1)
    {
      //TEST_FOR_EXCEPTION(expr->numElems() > 1, RuntimeError,
      //                   "FieldWriterBase::addField not ready for vector fields");

	  std::cout << "WARNING! : expr->numElems() > 1 , FieldWriterBase::addField only VTK can plot vector field " << std::endl;
	  std::cout << "WARNING! : All expressions(in the list of the expressions) must be of the same kind!!! " << std::endl;
	  /* Vector field plotting should be implemented for VTK files */
	  /* We assume that all the expressions are the same !!!! */
	  if (expr->isPointData()){
		  pointVectorFields_.append(expr);
		  pointVectorNames_.append(fieldName);
	  }else{
		  cellVectorFields_.append(expr);
		  cellVectorNames_.append(fieldName);
	  }

    } 
  else if (expr->isPointData()) 
    {
      /* expr is a single scalar field defined at points */
      pointScalarFields_.append(expr);
      pointScalarNames_.append(fieldName);
    }
  else if (expr->isCellData())
    {
      /* expr is a single scalar field defined on cells */
      cellScalarFields_.append(expr);
      cellScalarNames_.append(fieldName);
    }
}

void FieldWriterBase::addCommentLine(const std::string& line) 
{
  comments_.append(line);
}






