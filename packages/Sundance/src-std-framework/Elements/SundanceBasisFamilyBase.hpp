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

#ifndef SUNDANCE_BASISFAMILYBASE_H
#define SUNDANCE_BASISFAMILYBASE_H

#include "SundanceCellJacobianBatch.hpp"
#include "SundanceDefs.hpp"
#include "SundanceBasisDOFTopologyBase.hpp"
#include "SundanceTensorBasisBase.hpp"
#include "SundanceBasisReferenceEvaluationBase.hpp"
#include "SundanceHandleable.hpp"
#include "SundanceMesh.hpp"
#include "SundancePrintable.hpp"
#include "SundanceObjectWithVerbosity.hpp"
#include "SundanceTypeUtils.hpp"
#include "Teuchos_XMLObject.hpp"

namespace Sundance {

using namespace Teuchos;


/** 
 *
 */
class BasisFamilyBase
  : public Sundance::Handleable<BasisFamilyBase>,
    public Sundance::Printable,
    public ObjectWithClassVerbosity<BasisFamilyBase>,
    public BasisDOFTopologyBase,
    public TensorBasisBase,
    public BasisReferenceEvaluationBase
{
public:

  /** */
  virtual int order() const = 0 ;

  /** */
  virtual bool lessThan(const BasisDOFTopologyBase* other) const ;

  /** \brief Indicates whether mapping the basis requires an additional
      correction */
  virtual bool requiresBasisTransformation() const { return false; }
  /** \brief Default transformation does nothing */
  virtual void preApplyTransformation( const CellType &maxCellType ,
				       const Mesh &mesh,
				       const Array<int> &cellLIDs,
				       const CellJacobianBatch& JVol,
				       RCP<Array<double> >& A
				       ) const {;}
  /** \brief Default transformation does nothing */
  virtual void postApplyTransformation( const CellType &maxCellType ,
				       const Mesh &mesh,
				       const Array<int> &cellLIDs,
					const CellJacobianBatch& JVol,
					RCP<Array<double> >& A
					) const {;}
  /** \brief Default transformation does nothing */
  virtual void preApplyTransformationTranspose( const CellType &maxCellType ,
						const Mesh &mesh,
						const Array<int> &cellLIDs,
						const CellJacobianBatch& JVol,
						Array<double>& A ) const {;}

};


/* ----- Subtypes of BasisFamilyBase that specify transformation type --- */

/**     
 * Base class for scalar-valued basis families. Bases for scalar
 * fields living in, e.g., H1 or L2, should derive from this class. 
 */
class ScalarBasis 
  : public BasisFamilyBase
{
public:
  /** Inform caller that my tensor order is zero */
  int tensorOrder() const {return 0;}

  /** Inform caller that I am a scalar basis */
  bool isScalarBasis() const {return true;}
  
   /** Return the dimension of the members of a scalar-valued basis  */
  int dim() const {return 1;}
};

/** */
class VectorBasis
  : public BasisFamilyBase
{
public:
  /** */
  VectorBasis(int dim) : dim_(dim) {}
  /** Inform caller that my tensor order is one */
  int tensorOrder() const {return 1;}

   /** Return the dimension of the members of a scalar-valued basis  */
  int dim() const {return dim_;}

private:
  int dim_;
};


/** 
 * Base class for bases living in H(div) 
 */
class HDivVectorBasis
  : public VectorBasis
{
public:
  /** */
  HDivVectorBasis(int dim) : VectorBasis(dim) {}
  
  /** Inform caller that I am an H(div) basis */
  bool isHDivBasis() const {return true;}

};

/** 
 * Base class for bases living in H(curl) 
 */
class HCurlVectorBasis
  : public VectorBasis
{
public:
  /** */
  HCurlVectorBasis(int dim) : VectorBasis(dim) {}

  /** Inform caller that I am an H(curl) basis */
  bool isHCurlBasis() const {return true;}

};




} // namespace Sundance


#endif
