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

#ifndef TSFEPETRAVECTORSPACE_HPP
#define TSFEPETRAVECTORSPACE_HPP

#include "SundanceDefs.hpp"
#include "Epetra_Map.h"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_Comm.hpp"
#include "Thyra_ScalarProdVectorSpaceBase.hpp"
#include "Thyra_SpmdVectorSpaceBase.hpp"

namespace TSFExtended
{
using namespace Teuchos;
using namespace Thyra;


/**
 * Adaptor wrapping Epetra map in the Thyra vector space system.
 */
class EpetraVectorSpace 
  : virtual public Thyra::ScalarProdVectorSpaceBase<double>,
    virtual public Thyra::SpmdVectorSpaceBase<double>
{
public:

  /** */
  EpetraVectorSpace(const RCP<const Epetra_Map>& map);
    

  /** @name Overridden form Teuchos::Describable */
  //@{
  /** \brief . */
  std::string description() const;
  //@}

  /** @name Public overridden from VectorSpace */
  //@{

  /** */
  OrdType dim() const ;

  /** */
  bool isCompatible(const VectorSpaceBase<double>& other) const ;

  /** */
  RCP<const VectorSpaceFactoryBase<double> > 
  smallVecSpcFcty() const ;

  /** \brief clone the space */
  RCP< const VectorSpaceBase<double> > clone() const ;

  //@}




  /** \name Overidden from SpmdVectorSpaceBase */
  //@{
  /** */
  Teuchos::RCP<const Teuchos::Comm<OrdType> > getComm() const
    {return comm_;}

  /** */
  OrdType localSubDim() const {return localSubDim_;}

  /** */
  OrdType localOffset() const {return localOffset_;}

  /** */
  OrdType mapCode() const {return -1;}
  //@}


 /** */
  const RCP<const Epetra_Map>& epetraMap() const 
    {return epetraMap_;}


protected:
  /** */
  Teuchos::RCP<const Teuchos::Comm<OrdType> > 
  epetraCommToTeuchosComm(const Epetra_Comm& epComm) const ;
  
  /** @name Protected overridden from VectorSpace */
  //@{
  /** \brief create a vector */
  RCP<VectorBase<double> > createMember() const;
  /** \brief create a multivector */
  RCP<MultiVectorBase<double> > createMembers(int numVecs) const;
  //@}
private:
  /** */
  RCP<const VectorSpaceFactoryBase<double> > smallVecSpcFactory_;
  /** */
  RCP<const Epetra_Map> epetraMap_;

  Teuchos::RCP<const Teuchos::Comm<OrdType> > comm_;

  OrdType localSubDim_;

  OrdType localOffset_;
};
  
}

#endif
