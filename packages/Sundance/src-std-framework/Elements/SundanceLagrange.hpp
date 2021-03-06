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

#ifndef SUNDANCE_LAGRANGE_H
#define SUNDANCE_LAGRANGE_H

#include "SundanceDefs.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "SundanceBasisFamilyBase.hpp"

namespace Sundance 
{
/** 
 * Lagrange basis 
 */
class Lagrange : public ScalarBasis
{
public:
  /** */
  Lagrange(int order);

  /**   
   * \brief Inform caller as to whether a given cell type is supported 
   */
  bool supportsCellTypePair(
    const CellType& maximalCellType,
    const CellType& cellType
    ) const ;

  /** */
  void print(std::ostream& os) const ;

  /** */
  int order() const {return order_;}

  /** return the number of nodes for this basis on the given cell type */
  int nReferenceDOFsWithoutFacets(
    const CellType& maximalCellType,
    const CellType& cellType
    ) const ;

  /** */
  void getReferenceDOFs(
    const CellType& maximalCellType,
    const CellType& cellType,
    Array<Array<Array<int> > >& dofs) const ;

  /** */
  void refEval(
    const CellType& cellType,
    const Array<Point>& pts,
    const SpatialDerivSpecifier& deriv,
    Array<Array<Array<double> > >& result,
    int verbosity=0) const ;

  /** constraints for hanging DoFs*/
   void getConstrainsForHNDoF(
    	const int indexInParent,
    	const int maxCellDim,
    	const int maxNrChild,
    	const int facetDim,
    	const int facetIndex,
    	const int nodeIndex,
    	Array<int>& localDoFs,
	    Array<int>& parentFacetDim,
	    Array<int>& parentFacetIndex,
	    Array<int>& parentFacetNode,
    	Array<double>& coefs
    	);

  /* Handleable boilerplate */
  GET_RCP(BasisFamilyBase);

  /** */
  std::string description() const 
    {return "Lagrange(" + Teuchos::toString(order()) + ")";}

private:
  int order_;

  static Array<int> makeRange(int low, int high);

  /** evaluate on a line cell  */
  void evalOnLine(const Point& pt,
    const MultiIndex& deriv,
    Array<double>& result) const ;
    
  /** evaluate on a triangle cell  */
  void evalOnTriangle(const Point& pt,
    const MultiIndex& deriv,
    Array<double>& result) const ;

  /** evaluate on a quad cell  */
  void evalOnquad(const Point& pt,
    const MultiIndex& deriv,
    Array<double>& result) const ;

  /** evaluate on a tet cell  */
  void evalOnTet(const Point& pt,
    const MultiIndex& deriv,
    Array<double>& result) const ;

  /** evaluate on a tet cell  */
  void evalOnBrick(const Point& pt,
    const MultiIndex& deriv,
    Array<double>& result) const ;

  /** get the exact position one DoF on the Ref Element
   * this is needed for the treatment of hanging nodes
   * @param maxCellDim [in] the MaxCell dim of this element
   * @param facetDim   [in] the facet dim which the DoF is on
   * @param facetIndex [in] the facet index which the DoF is on
   * @param nodeIndex  [in] the node index of the DoF
   * @param pos       [out] the position of the DoF
   * */
  void returnDoFPositionOnRef(
	const int maxCellDim,
	const int facetDim,
	const int facetIndex,
	const int nodeIndex,
	Point& pos) const;

  /** This method calls the "getReferenceDOFs" method and then for each DoF
   * extracts the facet dimension, facet index and node in 3 different array <br>
   * This method could be used for general purpose for other basis, where HN is possible
   * @param maxCellDim [in]
   * @param nrDoF  [in] nr of DoF for this element
   * @param facetD [out]
   * @param facetI [out]
   * @param facetN [out] */
  void  getDoFsOrdered(
  		const CellType maxCellDim,
  		int nrDoF,
  		Array<int>& facetD,
  		Array<int>& facetI,
  		Array<int>& facetN);

  /** */
  bool doFInfromationCalculated_;

  /** For each DoF the dimension of the element which is the DoF on*/
  Array<int> facetD_;

  /** For each DoF the facet index of the element*/
  Array<int> facetI_;

  /** For each DoF the facet index of the element*/
  Array<int> facetN_;
};
}

#endif
