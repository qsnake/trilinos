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

#ifndef SUNDANCE_ELEMENTINTEGRALFUNCTIONALCELL_H
#define SUNDANCE_ELEMENTINTEGRALFUNCTIONALCELL_H

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
     * ElementIntegralFunctional provides an abstract base class for
     streaming functional evaluation over a batch of maximal cells (not facets).
     It provides a single pure virtual function for subclasses
     to implement.
    */
    class ElementIntegralFunctionalCell: public ElementIntegralFunctional
    {
    public:
      /** Constructor */
      ElementIntegralFunctionalCell( int spatialDim,
				     const CellType& maxCellType,
				     const QuadratureFamily& quad,
				     const ParameterList& verbParams 
				     = *ElementIntegralBase::defaultVerbParams()) : 
	ElementIntegralFunctional( spatialDim , maxCellType , 
				   spatialDim , maxCellType ,
				   quad , verbParams ){;}
      
      /** Destructor */
      virtual ~ElementIntegralFunctionalCell() {;}
      
      virtual void evaluate( CellJacobianBatch& JTrans,
			     const Array<int>& isLocalFlag,
			     const double* const coeff,
			     RefCountPtr<Array<double> >& A) const = 0;
    };
  }
}
#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif
