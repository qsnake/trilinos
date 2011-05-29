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

#ifndef SUNDANCE_IQI_LF_VN_FACET_H
#define SUNDANCE_IQI_LF_VN_FACET_H

#include "SundanceDefs.hpp"
#include "SundanceCellJacobianBatch.hpp"
#include "SundanceQuadratureFamily.hpp"
#include "SundanceBasisFamily.hpp"
#include "Teuchos_Array.hpp"
#include "SundanceElementIntegralLinearFormFacet.hpp"
#include "Intrepid_FieldContainer.hpp"

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
     * 
     IQI: IntrepidQuadratureIntegral
     HdivLF: LinearForm where test function lies in H(div)
     VN: the test function has its normal component taken
     Facet: this only makes sense on codimension 1 facets.

     This class evaluates integrals of the form
     \int_e f v.n
     where e is a facet of co-dimension 1, f is tabulated
     at quadrature points, and v is a test function in an H(div) space.
    */
    class IQI_HdivLF_VN_Facet : public ElementIntegralLinearFormFacet
    {
    public:
      /** Constructor */
      // this is only going to take maxCellType
      // as an argument, then use SundanceStdMesh::facetType(maxCellType,spatialDim-1,0)
      // to pass up the chain since we can only work on facets of codimension 1
      IQI_HdivLF_VN_Facet( int spatialDim,
			   const CellType& maxCellType,
			   const BasisFamily& testBasis,
			   const QuadratureFamily& quad,
			   const ParameterList& verbParams 
			   = *ElementIntegralBase::defaultVerbParams());

      /** Destructor */
      virtual ~IQI_HdivLF_VN_Facet() {;}
      
      /** evaluates integral of coeff against divergence of each basis function */
      virtual void evaluate( CellJacobianBatch& JTrans,
			     CellJacobianBatch& JVol ,
			     const Array<int>& facetIndex ,
			     const double* const coeff,
			     RefCountPtr<Array<double> >& A) const;


    };
  }
}
#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif
