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

#ifndef SUNDANCE_IQI_HDIVHDIVBLF_DIVUDIVV_CELL_H
#define SUNDANCE_IQI_HDIVHDIVBLF_DIVUDIVV_CELL_H


#include "SundanceDefs.hpp"
#include "SundanceCellJacobianBatch.hpp"
#include "SundanceQuadratureFamily.hpp"
#include "SundanceBasisFamily.hpp"
#include "Teuchos_Array.hpp"
#include "SundanceElementIntegralBilinearFormCell.hpp"


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
	VecVec bilinear form
	\int div(u) div(v)
	over a cell
	u is an unknown H(div) function
	v is an H(div) test function
    */
    class IQI_HdivHdivBLF_DivUDivV_Cell : public ElementIntegralBilinearFormCell
    {
    public:
      /** Constructor */
      IQI_HdivHdivBLF_DivUDivV_Cell( int spatialDim,
				 const CellType& maxCellType,
				 const BasisFamily &testBasis,
				 const BasisFamily &unkBasis ,
				 const QuadratureFamily& quad,
				 const ParameterList& verbParams 
				 = *ElementIntegralBase::defaultVerbParams() );
    
      /** Destructor */
      virtual ~IQI_HdivHdivBLF_DivUDivV_Cell() {;}
    
      virtual void evaluate( CellJacobianBatch& JTrans,
                             const double* const coeff,
                             RefCountPtr<Array<double> >& A) const = 0;
    



    };
  }
}
#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif
