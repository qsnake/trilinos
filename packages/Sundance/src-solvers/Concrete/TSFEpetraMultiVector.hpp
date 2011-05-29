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

#ifndef TSFEPETRAMULTIVECTOR_HPP
#define TSFEPETRAMULTIVECTOR_HPP

#include "SundanceDefs.hpp"
#include "SundancePrintable.hpp"
#include "TSFIndexableVector.hpp"
#include "TSFRawDataAccessibleVector.hpp"
#include "TSFVectorDecl.hpp"
#include "Epetra_FEVector.h"
#include "Epetra_Vector.h"

#ifdef TRILINOS_6
#include "Thyra_MPIVectorStdDecl.hpp"
#else
#define MPIMultiVectorStd DefaultMPIMultiVector
#include "Thyra_DefaultMPIMultiVectorDecl.hpp"
#endif

namespace TSFExtended
{
  using namespace Teuchos;
  using namespace Thyra;
  /**
   * TSF extension of Thyra::EpetraMultiVector, implementing the 
   * LoadableMultiVector
   * interface allowing an application to access elements. This class derives
   * from Thyra::EpetraVector, so it can be used seamlessly in any 
   * Thyra-based code.
   */
  class EpetraMultiVector : public MPIMultiVectorStd<double>,
                            public Sundance::Handleable<MultiVectorBase<double> >,
                            public IndexableVector<double>,
                            public RawDataAccessibleVector<double>,
                            public Printable
  {
  public:
    GET_RCP(MultiVectorBase<double>);

    /** Construct with a smart pointer to an Epetra vector space. */
    EpetraMultiVector(const RCP<const VectorSpaceBase<double> >& vs);

    /** Construct with smart pointers to an Epetra vector space
        and an existing Epetra vector. */
    EpetraMultiVector(const RCP<const VectorSpaceBase<double> >& vs,
                      const RCP<Epetra_MultiVector>& vec);

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
    void setElements(int numElems, const OrdType* globalIndices, 
                     const double* values);


    /** add to a group of elements */
    void addToElements(int numElems, const OrdType* globalIndices, 
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
                     vector<double>& elems) const ;
    //@}
      

    /** \name Printable interface */
    //@{
    /** Write to a stream  */
    void print(std::ostream& os) const 
    {
      epetraVec()->Print(os);
    }
    //@}

    /** Get a read-only Epetra_Vector */
    static const Epetra_Vector& getConcrete(const TSFExtended::Vector<double>& tsfVec);
    /** Get a read-write Epetra_Vector */
    static Epetra_Vector& getConcrete(TSFExtended::Vector<double>& tsfVec);
    /** Get a read-write Epetra_Vector pointer */
    static Epetra_Vector* getConcretePtr(TSFExtended::Vector<double>& tsfVec);



    
    /** */
    const RCP<Epetra_Vector>& epetraVec() const {return epetraVec_;}
    
    /** */
    RCP<Epetra_Vector>& epetraVec() {return epetraVec_;}

  protected:    
    /** */
    const RCP<const Epetra_Map>& epetraMap() const {return epetraMap_;}

  private:

    RCP<Epetra_Vector> epetraVec_;

    RCP<const MPIVectorSpaceBase<double> > mpiVecSpace_;

    RCP<const Epetra_Map> epetraMap_;
  };
  
}

#endif
