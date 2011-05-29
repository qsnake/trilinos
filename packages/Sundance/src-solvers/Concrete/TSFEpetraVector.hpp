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

#ifndef TSFEPETRAVECTOR_HPP
#define TSFEPETRAVECTOR_HPP

#include "SundanceDefs.hpp"
#include "SundancePrintable.hpp"
#include "TSFIndexableVector.hpp"
#include "TSFVectorDecl.hpp"
#include "TSFRawDataAccessibleVector.hpp"
#include "Thyra_VectorDefaultBase.hpp"
#include "Epetra_FEVector.h"
#include "Epetra_Vector.h"
#include "TSFEpetraVectorSpace.hpp"


namespace TSFExtended
{
using Teuchos::Range1D;
using namespace Thyra;
using namespace Teuchos;
/**
 * TSF extension of Thyra::EpetraVector, implementing the LoadableVector
 * interface allowing an application to access elements. This class derives
 * from Thyra::VectorDefaultBase, so it can be used seamlessly in any 
 * Thyra-based code.
 */
class EpetraVector : public Thyra::VectorDefaultBase<double>,
                     public IndexableVector<double>,
                     public RawDataAccessibleVector<double>,
                     public Sundance::Printable
{
public:

  /** Construct with a smart pointer to an Epetra vector space. */
  EpetraVector(const RCP<const VectorSpaceBase<double> >& vs);

  /** Construct with smart pointers to an Epetra vector space
      and an existing Epetra vector. */
  EpetraVector(const RCP<const VectorSpaceBase<double> >& vs,
    const RCP<Epetra_Vector>& vec);


  /** \name VectorBase interface */
  //@{
  /** */
   RCP< const VectorSpaceBase<double> > 
   space() const {return vecSpace_;}

#ifndef TRILINOS_8
  /** */
  void applyOpImpl(const RTOpPack::RTOpT< double >& op,
		const ArrayView< const Ptr< const VectorBase< double > > > &  	vecs,
		const ArrayView< const Ptr< VectorBase< double > > > &  	targ_vecs,
		const Ptr< RTOpPack::ReductTarget > &  	reduct_obj,
		const OrdType  	global_offset	 
    ) const ;
#else
  virtual void applyOp(
    const RTOpPack::RTOpT<double> &op,
    const int num_vecs,
    const VectorBase<double>*const vecs[],
    const int num_targ_vecs,
    VectorBase<double>*const targ_vecs[],
    RTOpPack::ReductTarget *reduct_obj,
    const OrdType first_ele_offset,
    const OrdType sub_dim,
    const OrdType global_offset
    ) const ;
#endif

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
  virtual const double* dataPtr() const {return &(epetraVec_->operator[](0));}
  /** */
  virtual double* dataPtr() {return &(epetraVec_->operator[](0));}
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
  const RCP<Epetra_Vector>& epetraVec() const 
    {return epetraVec_;}

  /** */
  RCP<Epetra_Vector>& epetraVec() {return epetraVec_;}

  /** Get a read-only Epetra_Vector */
  static const Epetra_Vector& getConcrete(const TSFExtended::Vector<double>& tsfVec);
  /** Get a read-write Epetra_Vector */
  static Epetra_Vector& getConcrete(TSFExtended::Vector<double>& tsfVec);
  /** Get a read-write Epetra_Vector pointer */
  static Epetra_Vector* getConcretePtr(TSFExtended::Vector<double>& tsfVec);

  

    

protected:    
  /** */
  const RCP<const Epetra_Map>& epetraMap() const {return epetraMap_;}

  /** */
  Range1D validateRange(const Range1D& rng) const ;

private:

  RCP<Epetra_Vector> epetraVec_;

  RCP<const Thyra::VectorSpaceBase<double> > vecSpace_;

  RCP<const EpetraVectorSpace> epetraVecSpace_;

  RCP<const Epetra_Map> epetraMap_;

  int localOffset_;

  int localSubDim_;

  int globalDim_;

  mutable bool in_applyOpImpl_;
};
  
}


#endif
