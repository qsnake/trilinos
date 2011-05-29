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

#ifndef SUNDANCE_RAVIARTTHOMAS_H
#define SUNDANCE_RAVIARTTHOMAS_H

#include "SundanceDefs.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "SundanceBasisFamilyBase.hpp"

namespace Sundance 
{
/** 
 * Lowest-order Raviart-Thomas basis 
 */
class RaviartThomas : public HDivVectorBasis
{
public:
  /** */
  RaviartThomas(int spatialDim);

  /** */
  int order() const {return 1;}

  /** */
  void print(std::ostream& os) const ;

  /**   
   * \brief Inform caller as to whether a given cell type is supported 
   */
  bool supportsCellTypePair(
    const CellType& maximalCellType,
    const CellType& cellType
    ) const ;

  /** */
  virtual void getReferenceDOFs(
    const CellType& maximalCellType,
    const CellType& cellType,
    Array<Array<Array<int> > >& dofs
    ) const ;

  /** */
  virtual int nReferenceDOFsWithoutFacets(
    const CellType& maximalCellType,
    const CellType& cellType
    ) const ;

  /** */
  virtual bool lessThan(const BasisDOFTopologyBase* other) const ;

  /** */
  std::string description() const ;

  /** */
  virtual void refEval(
    const CellType& cellType,
    const Array<Point>& pts,
    const SpatialDerivSpecifier& deriv,
    Array<Array<Array<double> > >& result,
    int verbosity=0) const ;  


  /* Handleable boilerplate */
  GET_RCP(BasisFamilyBase);
private:

    
};
}

#endif
