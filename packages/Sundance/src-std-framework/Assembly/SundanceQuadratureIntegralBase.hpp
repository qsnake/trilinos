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

#ifndef SUNDANCE_QUADRATUREINTEGRALBASE_H
#define SUNDANCE_QUADRATUREINTEGRALBASE_H

#include "SundanceDefs.hpp"
#include "SundanceElementIntegral.hpp"

namespace Sundance
{
using namespace Teuchos;
    
/** 
 *  
 *
 */
class QuadratureIntegralBase
  : public ElementIntegral
{
public:

  /** Construct a zero form to be computed by quadrature */
  QuadratureIntegralBase(int spatialDim,
    const CellType& maxCellType,
    int dim, 
    const CellType& cellType,
    const QuadratureFamily& quad,
    bool isInternalBdry,
    const ParametrizedCurve& globalCurve,
    const Mesh& mesh,
    int verb);

  /** Construct a one form to be computed by quadrature */
  QuadratureIntegralBase(int spatialDim,
    const CellType& maxCellType,
    int dim, 
    const CellType& cellType,
    const BasisFamily& testBasis,
    int alpha,
    int testDerivOrder,
    const QuadratureFamily& quad,
    bool isInternalBdry,
    const ParametrizedCurve& globalCurve,
    const Mesh& mesh,
    int verb);
      
  /** Construct a two-form to be computed by quadrature */
  QuadratureIntegralBase(int spatialDim,
    const CellType& maxCellType,
    int dim,
    const CellType& cellType,
    const BasisFamily& testBasis,
    int alpha,
    int testDerivOrder,
    const BasisFamily& unkBasis,
    int beta,
    int unkDerivOrder,
    const QuadratureFamily& quad,
    bool isInternalBdry,
    const ParametrizedCurve& globalCurve,
    const Mesh& mesh,
    int verb);
      
  /** virtual dtor */
  virtual ~QuadratureIntegralBase(){;}
      
     /** */
  virtual void transform(const CellJacobianBatch& JTrans,
    const CellJacobianBatch& JVol,
    const Array<int>& isLocalFlag,
    const Array<int>& facetNum,
    const RCP<Array<int> >& cellLIDs,
    const double* const coeff,
    RCP<Array<double> >& A) const 
    {
      if (order()==2) transformTwoForm(JTrans, JVol, facetNum, cellLIDs,coeff, A);
      else if (order()==1) transformOneForm(JTrans, JVol, facetNum, cellLIDs,coeff, A);
      else transformZeroForm(JTrans, JVol, isLocalFlag, facetNum, cellLIDs,coeff, A);
    }
      
  /** */
  virtual void transformZeroForm(const CellJacobianBatch& JTrans,
    const CellJacobianBatch& JVol,
    const Array<int>& isLocalFlag,
    const Array<int>& facetIndex,
    const RCP<Array<int> >& cellLIDs,
    const double* const coeff,
    RCP<Array<double> >& A) const = 0;
      
  /** */
  virtual void transformTwoForm(const CellJacobianBatch& JTrans,
    const CellJacobianBatch& JVol,
    const Array<int>& facetIndex,
    const RCP<Array<int> >& cellLIDs,
    const double* const coeff,
    RCP<Array<double> >& A) const = 0;
      
  /** */
  virtual void transformOneForm(const CellJacobianBatch& JTrans,
    const CellJacobianBatch& JVol,
    const Array<int>& facetIndex,
    const RCP<Array<int> >& cellLIDs,
    const double* const coeff,
    RCP<Array<double> >& A) const = 0;
      
      
  /** */
  virtual int nQuad() const {return nQuad_;}
      
  static double& totalFlops() {static double rtn = 0; return rtn;}

protected:
  static void addFlops(const double& flops) {totalFlops() += flops;}
      
  const QuadratureFamily& quad() const {return quad_;}
  /** */
  int nQuad_ ;

private:

  /* */
  QuadratureFamily quad_;
      
};

}


#endif
