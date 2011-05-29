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

#ifndef SUNDANCE_VECTORBASISCOMPONENT_H
#define SUNDANCE_VECTORBASISCOMPONENT_H

#include "SundanceDefs.hpp"
#include "SundanceBasisFamily.hpp"

namespace Sundance {

using namespace Teuchos;

/** 
 * This class is for the representation of a single component
 * of a vector-valued basis family. 
 */
class VectorBasisComponent : public BasisFamilyBase
{
public:
  /** */
  VectorBasisComponent(const BasisFamily& master, int direction);

  /** */
  bool lessThan(const BasisFamilyBase* other) const ;

  /** */
  int order() const {return master_.order();}

  /** */
  int dim() const 
    {return master_.dim();}

  /** */
  bool isCovariantBasis() const 
    {return master_.isCovariantBasis();}

  /** */
  bool isContravariantBasis() const 
    {return master_.isContravariantBasis();}

  /** */
  int direction() const {return direction_;}

  /** */
  bool supportsCellTypePair(
    const CellType& maximalCellType,
    const CellType& cellType
    ) const
    {
      return master_.ptr()->supportsCellTypePair(maximalCellType, 
        cellType);
    }


  /** */
  void getReferenceDOFs(
    const CellType& maximalCellType,
    const CellType& cellType,
    Array<Array<Array<int> > >& dofs
    ) const 
    {
      master_.ptr()->getReferenceDOFs(maximalCellType, 
        cellType, dofs);
    }


  /** */
  int nReferenceDOFsWithFacets(
    const CellType& maximalCellType,
    const CellType& cellType
    ) const
    {
      return master_.ptr()->nReferenceDOFsWithFacets(maximalCellType, 
        cellType);
    }

  /** */
  int nReferenceDOFsWithoutFacets(
    const CellType& maximalCellType,
    const CellType& cellType
    ) const
    {
      return master_.ptr()->nReferenceDOFsWithoutFacets(maximalCellType, 
        cellType);
    }

  /** */
  void refEval(
    const CellType& maximalCellType,
    const CellType& cellType,
    const Array<Point>& pts,
    const MultiIndex& deriv,
    Array<Array<Array<double> > >& result
    ) const
    {
      master_.ptr()->refEval(maximalCellType, cellType,
        pts, deriv, result);
    }


private:
  BasisFamily master_;
  int direction_;
};


} // namespace Sundance


#endif
