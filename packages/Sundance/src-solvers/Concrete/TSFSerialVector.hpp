/* ***********************************************************************
// 
//           TSFExtended: Trilinos Solver Framework Extended
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// **********************************************************************/

#ifndef TSF_SERIAL_VECTOR_HPP
#define TSF_SERIAL_VECTOR_HPP

#include "SundanceDefs.hpp"
#include "SundancePrintable.hpp"
#include "TSFIndexableVector.hpp"
#include "TSFVectorDecl.hpp"
#include "TSFRawDataAccessibleVector.hpp"
#include "Thyra_VectorDefaultBase.hpp"
#include "TSFSerialVectorSpace.hpp"


namespace TSFExtended
{
using Teuchos::Range1D;
using namespace Thyra;
using namespace Teuchos;
/**
 * TSF implementation of a serial vector, implementing the LoadableVector
 * interface allowing an application to access elements. This class derives
 * from Thyra::VectorDefaultBase, so it can be used seamlessly in any 
 * Thyra-based code. If created in SPMD, this will be replicated on
 * all processors.
 */
class SerialVector : public Thyra::VectorDefaultBase<double>,
                     public IndexableVector<double>,
                     public RawDataAccessibleVector<double>,
                     public Sundance::Printable
{
public:

  /** Construct with a smart pointer to a vector space. */
  SerialVector(const RCP<const VectorSpaceBase<double> >& vs);

  /** \name VectorBase interface */
  //@{
  /** */
   RCP< const VectorSpaceBase<double> > 
   space() const {return vecSpace_;}

  /** */
  void applyOpImpl(const RTOpPack::RTOpT< double >& op,
		const ArrayView< const Ptr< const VectorBase< double > > > &  	vecs,
		const ArrayView< const Ptr< VectorBase< double > > > &  	targ_vecs,
		const Ptr< RTOpPack::ReductTarget > &  	reduct_obj,
		const OrdType  	global_offset	 
    ) const ;

  /** */
  void acquireDetachedVectorViewImpl(const Range1D& rng,
		RTOpPack::ConstSubVectorView<double>* sub_vec) const ;

  /** */
  void releaseDetachedVectorViewImpl(
    RTOpPack::ConstSubVectorView<double>* sub_vec) const ;

  /** */
  void acquireNonconstDetachedVectorViewImpl(const Range1D& rng,
		RTOpPack::SubVectorView<double> * sub_vec);	 


  /** */
  void commitNonconstDetachedVectorViewImpl(
    RTOpPack::SubVectorView<double>* sub_vec);
  
  //@}

  /** \name IndexableVector interface */
  //@{
  /** read the element at the given global index */
  virtual const double& operator[](OrdType globalIndex) const 
    {return getElement(globalIndex);}

  /** writable access to the element at the given global index */
  virtual double& operator[](OrdType globalIndex) ;
  //@}

  /** \name Raw data access interface */
  //@{
  /** */
  virtual const double* dataPtr() const {return &(data_[0]);}
  /** */
  virtual double* dataPtr() {return &(data_[0]);}
  //@}

  /** \name LoadableVector interface */
  //@{
  /** set a single element */
  void setElement(OrdType globalIndex, const double& value);

  /** add to a single element */
  void addToElement(OrdType globalIndex, const double& value);

  /** set a group of elements */
  void setElements(int numElems, const int* globalIndices, 
    const double* values);


  /** add to a group of elements */
  void addToElements(int numElems, const int* globalIndices, 
    const double* values);

  /** */
  void finalizeAssembly();
  //@}

  /** \name AccessibleVector interface */
  //@{
  /** */
  const double& getElement(OrdType globalIndex) const ;

  /** */
  void getElements(const OrdType* globalIndices, int numElems,
    Teuchos::Array<double>& elems) const ;
  //@}

  /** \name Printable interface */
  //@{
  /** Print to a stream */
  void print(std::ostream& os) const ;
  //@}

  /** */
  static const SerialVector* getConcrete(const Vector<double>& x);
  /** */
  static SerialVector* getConcrete(Vector<double>& x);
      

protected:    
  /** */
  Range1D validateRange(const Range1D& rng) const ;

private:

  RCP<const Thyra::VectorSpaceBase<double> > vecSpace_;

  Array<double> data_;

  int globalDim_;

  mutable bool in_applyOpImpl_;
};
  
}


#endif
