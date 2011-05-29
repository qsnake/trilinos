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

#ifndef SUNDANCE_REFINTEGRAL_H
#define SUNDANCE_REFINTEGRAL_H

#include "SundanceDefs.hpp"
#include "SundanceElementIntegral.hpp"


namespace Sundance
{

using namespace Teuchos;

/** 
 * RefIntegral represents the integrals of a product 
 * of basis functions (or their derivatives) over a reference cell. 
 * This object can be created once and for all up front, and then 
 * reused to produce the integrals of constant-coefficient weak forms
 * over simplicial cells by means of a linear transformation.
 * 
 * An instance of this object can represent either a one form,
 * \f[
 * W_{i\gamma} = \int_{T_R} D_\gamma \psi_i 
 * \f]  
 * or a two form,
 * \f[
 * W_{(ij)(\gamma\delta)} = \int_{T_R} D_\gamma \psi_i D_\delta \phi_j.
 * \f]  
 * In the notation for the two form, we have grouped together the index
 * pairs \f$(ij)\f$ and \f$(\gamma\delta)\f$ to indicate that we 
 * will think of this 4-tensor as a 2D array, with \f$(ij)\f$ combinations
 * running down rows and \f$(\gamma\delta)\f$ running across columns. 
 * We will consider the node number combination 
 * \f$(ij)\f$ to be a single index \f$k\f$, and 
 * \f$(\gamma\delta)\f$ to be a single index \f$\epsilon\f$. Thus, 
 * the storage and use of one-forms and two-forms is essentially identical.
 * 
 * This object's job in life is to be multiplied with a transformation
 * matrix \f$T_{\epsilon c}\f$to produce an array of local element matrices
 * \f$A_{kc}\f$,
 * \f[
 * A_{kc} = W_{k\epsilon} T_{\epsilon c}
 * \f]
 * The index \f$c\f$ is over cells (cells are processed in batches).
 * 
 * Physical storage is as a 1D vector stored in colum-major order. 
 */
class RefIntegral : public ElementIntegral
{
public:
  /** Construct a reference zero-form */
  RefIntegral(int spatialDim,
    const CellType& maxCellType,
    int dim, 
    const CellType& cellType,
    const QuadratureFamily& quad_in,
    bool isInternalBdry,
    const ParametrizedCurve& globalCurve,
    const Mesh& mesh,
    int verb);

  /** Construct a reference one-form */
  RefIntegral(int spatialDim,
    const CellType& maxCellType,
    int dim, 
    const CellType& cellType,
    const BasisFamily& testBasis,
    int alpha,
    int testDerivOrder,
    const QuadratureFamily& quad_in,
    bool isInternalBdry,
    const ParametrizedCurve& globalCurve,
    const Mesh& mesh,
    int verb);

  /** Construct a reference two-form */
  RefIntegral(int spatialDim,
    const CellType& maxCellType,
    int dim,
    const CellType& cellType,
    const BasisFamily& testBasis,
    int alpha,
    int testDerivOrder,
    const BasisFamily& unkBasis,
    int beta,
    int unkDerivOrder,
    const QuadratureFamily& quad_in,
    bool isInternalBdry,
    const ParametrizedCurve& globalCurve,
    const Mesh& mesh,
    int verb);

  /** virtual dtor */
  virtual ~RefIntegral(){;}
      
  /** */
  void print(std::ostream& os) const ;

  /** */
  void transform(const CellJacobianBatch& JTrans,
    const CellJacobianBatch& JVol,
    const Array<int>& isLocalFlag,
    const Array<int>& facetNum,
    const RCP<Array<int> >& cellLIDs,
    const double& coeff,
    RCP<Array<double> >& A) const
    {
      if (order()==2) transformTwoForm(JTrans, JVol, facetNum, cellLIDs, coeff, A);
      else if (order()==1) transformOneForm(JTrans, JVol, facetNum, cellLIDs, coeff, A);
      else transformZeroForm(JVol, isLocalFlag, cellLIDs, coeff, A);
    }

  /** */
  void transformTwoForm(const CellJacobianBatch& JTrans,
    const CellJacobianBatch& JVol,
    const Array<int>& facetNum, 
    const RCP<Array<int> >& cellLIDs,
    const double& coeff,
    RCP<Array<double> >& A) const ;

  /** */
  void transformOneForm(const CellJacobianBatch& JTrans,
    const CellJacobianBatch& JVol,
    const Array<int>& facetNum, 
    const RCP<Array<int> >& cellLIDs,
    const double& coeff,
    RCP<Array<double> >& A) const ;

  /** */
  void transformZeroForm(const CellJacobianBatch& JVol,
    const Array<int>& isLocalFlag,
    const RCP<Array<int> >& cellLIDs,
    const double& coeff,
    RCP<Array<double> >& A) const ;

  /** */
  inline double& value(int facetCase, int testDerivDir, int testNode,
    int unkDerivDir, int unkNode)
    {return W_[facetCase][unkNode + nNodesUnk()*testNode 
        + nNodes()*(unkDerivDir 
          + nRefDerivUnk()*testDerivDir)];}

  /** */
  inline const double& value(int facetCase, 
    int testDerivDir, int testNode,
    int unkDerivDir, int unkNode) const 
    {
      return W_[facetCase][unkNode + nNodesUnk()*testNode 
        + nNodes()*(unkDerivDir 
          + nRefDerivUnk()*testDerivDir)];
    }
      
  /** */
  inline double& value(int facetCase, int testDerivDir, int testNode)
    {return W_[facetCase][nNodesTest()*testDerivDir + testNode];}

  /** */
  inline const double& value(int facetCase, 
    int testDerivDir, int testNode) const 
    {return W_[facetCase][nNodesTest()*testDerivDir + testNode];}

  static double& totalFlops() {static double rtn = 0; return rtn;}

protected:

  static void addFlops(const double& flops) {totalFlops() += flops;}
      
private:

  Array<Array<double> > W_;

  /** For ACI (ACI = Adaptive Cell Integration), store the reference integral values for one form
   * The indexes facet, quadPoints, nRefDerivTest , nNodesTest */
  Array<Array<Array<Array<double> > > > W_ACI_F1_;

  /** For ACI (ACI = Adaptive Cell Integration), store the reference integral values for two form
   * The indexes facet, quadPoints, nRefDerivTest , nNodesTest , nRefDerivUnk , nNodesUnk */
  Array<Array<Array<Array<Array<Array<double> > > > > > W_ACI_F2_;

  /** The quadrature family needed for special integration */
  QuadratureFamily quad_;

  /** The quadrature points*/
  Array<Point> quadPts_;

  /** The standard weights (in case of ACI these might change)*/
  Array<double> quadWeights_;
};
}

#endif
