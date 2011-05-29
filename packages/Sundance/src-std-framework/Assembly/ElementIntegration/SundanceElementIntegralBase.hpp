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

#ifndef SUNDANCE_ELEMENTINTEGRALBASE_H
#define SUNDANCE_ELEMENTINTEGRALBASE_H

#include "SundanceDefs.hpp"
#include "SundanceCellJacobianBatch.hpp"
#include "SundanceQuadratureFamily.hpp"
#include "SundanceBasisFamily.hpp"
#include "Teuchos_Array.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceStdFwk
{
  using namespace SundanceUtils;
  using namespace SundanceStdMesh;
  using namespace SundanceStdMesh::Internal;
  using namespace SundanceCore;
  using namespace SundanceCore::Internal;
  
  namespace Internal
  {
    using namespace Teuchos;
    
    /** 
     * ElementIntegralFunctional provides an abstract base class for storing
     * information about the geometry/topology over which the form is defined
     * and the quadrature rule.
     * Computational routines are provided by subclasses; this is for storage
     * and information only, providing no particular interface to integration.
    */
    class ElementIntegralBase
      : public TSFExtended::ParameterControlledObjectWithVerbosity<ElementIntegralBase>,
	public TSFExtended::Printable
    {
    public:
      /** Constructor */
      ElementIntegralBase( int spatialDim,
			   const CellType& maxCellType,
			   int dim, 
			   const CellType& cellType,
			   const QuadratureFamily& quad,
			   const ParameterList& verbParams 
			   = *ElementIntegralBase::defaultVerbParams()):
	ParameterControlledObjectWithVerbosity<ElementIntegralBase>("Integration",verbParams),
	spatialDim_( spatialDim ), 
	maxCellType_( maxCellType ) ,
	dim_( dim ),
	cellType_( cellType ),
	quad_( quad ),
	verbParams_( verbParams ) {;}
	
      
      /** Destructor */
      virtual ~ElementIntegralBase() {;}
      
      virtual int spatialDim() const { return spatialDim_; }
      virtual int dim() const { return dim_; }
      virtual const CellType & maxCellType() const { return maxCellType_; }
      virtual const CellType & cellType() const { return cellType_; }
      virtual const QuadratureFamily & quad() const { return quad_; }
      virtual const ParameterList &verbParams() const { return verbParams_; }
 
      static RefCountPtr<ParameterList> defaultVerbParams()
      {
	static RefCountPtr<ParameterList> rtn = rcp(new ParameterList("Integration"));
	static int first = true;
	if (first)
	  {
	    rtn->set<int>("setup", 0);
	    rtn->set<int>("transformation", 0);
	    rtn->set<int>("integration", 0);
	    rtn->set<int>("extract weak form", 0);
	    rtn->set<int>("find groups", 0);
	    first = false;
	  }
	return rtn;
      }
      
    private:
      const int spatialDim_;
      const CellType &maxCellType_;
      const int dim_;
      const CellType &cellType_;
      const QuadratureFamily & quad_;
      const ParameterList & verbParams_;
      
    };
  }
}
#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif
