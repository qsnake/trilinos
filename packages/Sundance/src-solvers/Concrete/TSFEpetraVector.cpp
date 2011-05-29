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

#include "TSFEpetraVector.hpp"
#include "TSFEpetraVectorSpace.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_dyn_cast.hpp"
#include "Teuchos_Workspace.hpp"
#include "Thyra_DefaultSpmdVector.hpp"
#ifdef HAVE_MPI
#include "Teuchos_DefaultMpiComm.hpp"
#endif
#include "Teuchos_DefaultSerialComm.hpp"

#include "Thyra_apply_op_helper.hpp"
#include "RTOpPack_SPMD_apply_op.hpp"
#include "RTOp_parallel_helpers.h"


#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "TSFVectorImpl.hpp"
#include "TSFLinearOperatorImpl.hpp"
#endif

using namespace Teuchos;
using namespace TSFExtended;
using namespace Thyra;

using Teuchos::Range1D;


EpetraVector::EpetraVector(
  const RCP<const VectorSpaceBase<double> >& vs)
  : VectorDefaultBase<double>(), 
    epetraVec_(), 
    vecSpace_(vs), 
    epetraVecSpace_(), 
    epetraMap_(),
    localOffset_(),
    localSubDim_(),
    globalDim_(vs->dim()),
    in_applyOpImpl_(false)
{
  const EpetraVectorSpace* epvs 
    = dynamic_cast<const EpetraVectorSpace*>(vs.get());
  TEST_FOR_EXCEPTION(epvs==0, std::runtime_error,
    "could not cast vector space to EpetraVectorSpace in "
    "EpetraVector ctor");

  epetraVecSpace_ = rcp_dynamic_cast<const EpetraVectorSpace>(vs);
  epetraMap_ = epvs->epetraMap();
  epetraVec_ = rcp(new Epetra_Vector(*epetraMap_, true));

  localOffset_ = epetraMap_->MinMyGID();
  localSubDim_ = epetraMap_->NumMyElements();
}



EpetraVector
::EpetraVector(const RCP<const VectorSpaceBase<double> >& vs,
  const RCP<Epetra_Vector>& vec)
  : VectorDefaultBase<double>(), 
    epetraVec_(vec), 
    vecSpace_(vs), 
    epetraVecSpace_(), 
    epetraMap_(),
    localOffset_(),
    localSubDim_(),
    globalDim_(vs->dim()),
    in_applyOpImpl_(false)
{
  const EpetraVectorSpace* epvs 
    = dynamic_cast<const EpetraVectorSpace*>(vs.get());
  TEST_FOR_EXCEPTION(epvs==0, std::runtime_error,
    "could not cast vector space to EpetraVectorSpace in "
    "EpetraVector ctor");


  epetraVecSpace_ = rcp_dynamic_cast<const EpetraVectorSpace>(vs);
  epetraMap_ = epvs->epetraMap();
  if (vec.get() == 0) epetraVec_ = rcp(new Epetra_Vector(*epetraMap_, true));

  localOffset_ = epetraMap_->MinMyGID();
  localSubDim_ = epetraMap_->NumMyElements();
}





#ifndef TRILINOS_8

void EpetraVector::applyOpImpl(const RTOpPack::RTOpT< double >& op,
		const ArrayView< const Ptr< const VectorBase< double > > > &  	vecs,
		const ArrayView< const Ptr< VectorBase< double > > > &  	targ_vecs,
		const Ptr< RTOpPack::ReductTarget > &  	reduct_obj,
		const OrdType  	global_offset_in	 
  ) const 
{
  
  // ToDo: Remove this!
  const OrdType	first_ele_offset_in = 0;
  const OrdType	sub_dim_in = -1;

  using Teuchos::null;
  using Teuchos::dyn_cast;
  using Teuchos::Workspace;

  const int num_vecs = vecs.size();
  const int num_targ_vecs = targ_vecs.size();

#ifdef THYRA_SPMD_VECTOR_BASE_DUMP
  Teuchos::RCP<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();
  Teuchos::OSTab tab(out);
  if(show_dump) {
    *out << "\nEntering SpmdVectorBase<double>::applyOp(...) ...\n";
    *out
      << "\nop = " << typeName(op)
      << "\nnum_vecs = " << num_vecs
      << "\nnum_targ_vecs = " << num_targ_vecs
      << "\nreduct_obj = " << reduct_obj
      << "\nfirst_ele_offset_in = " << first_ele_offset_in
      << "\nsub_dim_in = " << sub_dim_in
      << "\nglobal_offset_in = " << global_offset_in
      << "\n"
      ;
  }
#endif // THYRA_SPMD_VECTOR_BASE_DUMP

  Ptr<Teuchos::WorkspaceStore> wss = Teuchos::get_default_workspace_store().ptr();

#ifdef TEUCHOS_DEBUG
  // ToDo: Validate input!
  TEST_FOR_EXCEPTION(
    in_applyOpImpl_, std::invalid_argument
    ,"SpmdVectorBase<>::applyOp(...): Error, this method is being entered recursively which is a "
    "clear sign that one of the methods acquireDetachedView(...), releaseDetachedView(...) or commitDetachedView(...) "
    "was not implemented properly!"
    );
  Thyra::apply_op_validate_input(
    "SpmdVectorBase<>::applyOp(...)",*space(),
    op, vecs, targ_vecs, reduct_obj,
    global_offset_in
    );
#endif

  const Teuchos::RCP<const Teuchos::Comm<OrdType> >& comm
    = epetraVecSpace_->getComm();

  const SerialComm<OrdType>* serialComm = 
    dynamic_cast<const SerialComm<OrdType>* >(comm.get());


  // Flag that we are in applyOp()
  in_applyOpImpl_ = true;

  // First see if this is a locally replicated vector in which case
  // we treat this as a local operation only.
  const bool locallyReplicated = ( serialComm != 0 && localSubDim_ == globalDim_ );

  // Get the overlap in the current process with the input logical sub-vector
  // from (first_ele_offset_in,sub_dim_in,global_offset_in)
  OrdType  overlap_first_local_ele_off = 0;
  OrdType  overlap_local_sub_dim = 0;
  OrdType  overlap_global_off = 0;
  if(localSubDim_) {
    RTOp_parallel_calc_overlap(
      globalDim_, localSubDim_, localOffset_,
      first_ele_offset_in, sub_dim_in, global_offset_in,
      &overlap_first_local_ele_off, &overlap_local_sub_dim, &overlap_global_off
      );
  }
  const Range1D local_rng = (
    overlap_first_local_ele_off>=0
    ? Range1D( localOffset_ + overlap_first_local_ele_off, localOffset_
      + overlap_first_local_ele_off + overlap_local_sub_dim - 1 )
    : Range1D::Invalid
    );

#ifdef THYRA_SPMD_VECTOR_BASE_DUMP
  if(show_dump) {
    *out
      << "\noverlap_first_local_ele_off = " << overlap_first_local_ele_off
      << "\noverlap_local_sub_dim = " << overlap_local_sub_dim
      << "\noverlap_global_off = " << overlap_global_off
      << "\nlocal_rng = ["<<local_rng.lbound()<<","<<local_rng.ubound()<<"]"
      << "\n"
      ;
  }
#endif // THYRA_SPMD_VECTOR_BASE_DUMP

  // Create sub-vector views of all of the *participating* local data
  Workspace<RTOpPack::ConstSubVectorView<double> > sub_vecs(wss.get(),num_vecs);
  Workspace<RTOpPack::SubVectorView<double> > sub_targ_vecs(wss.get(),num_targ_vecs);
  if( overlap_first_local_ele_off >= 0 ) {
    {for(int k = 0; k < num_vecs; ++k ) {
      vecs[k]->acquireDetachedView( local_rng, &sub_vecs[k] );
      sub_vecs[k].setGlobalOffset( overlap_global_off );
    }}
    {for(int k = 0; k < num_targ_vecs; ++k ) {
      targ_vecs[k]->acquireDetachedView( local_rng, &sub_targ_vecs[k] );
      sub_targ_vecs[k].setGlobalOffset( overlap_global_off );
    }}
  }

  

  // Apply the RTOp operator object (all processors must participate)
  RTOpPack::SPMD_apply_op(
    locallyReplicated ? NULL : &*comm,     // comm
    op,                                    // op
    num_vecs,                              // num_vecs
    sub_vecs.getRawPtr(),                  // sub_vecs
    num_targ_vecs,                         // num_targ_vecs
    sub_targ_vecs.getRawPtr(),             // targ_sub_vecs
    reduct_obj.get()                       // reduct_obj
    );

  // Free and commit the local data
  if (overlap_first_local_ele_off >= 0) {
    for (int k = 0; k < num_vecs; ++k ) {
      sub_vecs[k].setGlobalOffset(local_rng.lbound());
      vecs[k]->releaseDetachedView( &sub_vecs[k] );
    }
    for (int k = 0; k < num_targ_vecs; ++k ) {
      sub_targ_vecs[k].setGlobalOffset(local_rng.lbound());
      targ_vecs[k]->commitDetachedView( &sub_targ_vecs[k] );
    }
  }

  // Flag that we are leaving applyOp()
  in_applyOpImpl_ = false;

#ifdef THYRA_SPMD_VECTOR_BASE_DUMP
  if(show_dump) {
    *out << "\nLeaving SpmdVectorBase<double>::applyOp(...) ...\n";
  }
#endif // THYRA_SPMD_VECTOR_BASE_DUMP

}

#else

void EpetraVector::applyOp(
  const RTOpPack::RTOpT<double> &op,
  const int num_vecs,
  const VectorBase<double>*const vecs[],
  const int num_targ_vecs,
  VectorBase<double>*const targ_vecs[],
  RTOpPack::ReductTarget *reduct_obj,
  const OrdType first_ele_offset_in,
  const OrdType sub_dim_in,
  const OrdType global_offset_in
  ) const 
{
#ifdef THYRA_SPMD_VECTOR_BASE_DUMP
  Teuchos::RCP<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();
  Teuchos::OSTab tab(out);
  if(show_dump) {
    *out << "\nEntering SpmdVectorBase<double>::applyOp(...) ...\n";
    *out
      << "\nop = " << typeName(op)
      << "\nnum_vecs = " << num_vecs
      << "\nnum_targ_vecs = " << num_targ_vecs
      << "\nreduct_obj = " << reduct_obj
      << "\nfirst_ele_offset_in = " << first_ele_offset_in
      << "\nsub_dim_in = " << sub_dim_in
      << "\nglobal_offset_in = " << global_offset_in
      << "\n"
      ;
  }
#endif // THYRA_SPMD_VECTOR_BASE_DUMP
  using Teuchos::dyn_cast;
  using Teuchos::Workspace;
  Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();

#ifdef TEUCHOS_DEBUG
  // ToDo: Validate input!
  TEST_FOR_EXCEPTION(
    in_applyOpImpl_, std::invalid_argument
    ,"SpmdVectorBase<>::applyOp(...): Error, this method is being entered recursively which is a "
    "clear sign that one of the methods acquireDetachedView(...), releaseDetachedView(...) or commitDetachedView(...) "
    "was not implemented properly!"
    );
  Thyra::apply_op_validate_input(
    "SpmdVectorBase<>::applyOp(...)",*space()
    ,op,num_vecs,vecs,num_targ_vecs,targ_vecs,reduct_obj,first_ele_offset_in,sub_dim_in,global_offset_in
    );
#endif  

  const Teuchos::RCP<const Teuchos::Comm<OrdType> >& comm
    = epetraVecSpace_->getComm();

  const SerialComm<OrdType>* serialComm = 
    dynamic_cast<const SerialComm<OrdType>* >(comm.get());


  // Flag that we are in applyOp()
  in_applyOpImpl_ = true;

  // First see if this is a locally replicated vector in which case
  // we treat this as a local operation only.
  const bool locallyReplicated = ( serialComm != 0 && localSubDim_ == globalDim_ );
  
  // Get the overlap in the current process with the input logical sub-vector
  // from (first_ele_offset_in,sub_dim_in,global_offset_in)
  OrdType  overlap_first_local_ele_off  = 0;
  OrdType  overlap_local_sub_dim        = 0;
  OrdType  overlap_global_off           = 0;
  if(localSubDim_) {
    RTOp_parallel_calc_overlap(
      globalDim_, localSubDim_, localOffset_, first_ele_offset_in, sub_dim_in, global_offset_in
      ,&overlap_first_local_ele_off, &overlap_local_sub_dim, &overlap_global_off
      );
  }
  const Range1D local_rng = (
    overlap_first_local_ele_off>=0
    ? Range1D( localOffset_ + overlap_first_local_ele_off, localOffset_ + overlap_first_local_ele_off + overlap_local_sub_dim - 1 )
    : Range1D::Invalid
    );
#ifdef THYRA_SPMD_VECTOR_BASE_DUMP
  if(show_dump) {
    *out
      << "\noverlap_first_local_ele_off = " << overlap_first_local_ele_off
      << "\noverlap_local_sub_dim = " << overlap_local_sub_dim
      << "\noverlap_global_off = " << overlap_global_off
      << "\nlocal_rng = ["<<local_rng.lbound()<<","<<local_rng.ubound()<<"]"
      << "\n"
      ;
  }
#endif // THYRA_SPMD_VECTOR_BASE_DUMP
  // Create sub-vector views of all of the *participating* local data
  Workspace<RTOpPack::ConstSubVectorView<double> > sub_vecs(wss,num_vecs);
  Workspace<RTOpPack::SubVectorView<double> > sub_targ_vecs(wss,num_targ_vecs);
  if( overlap_first_local_ele_off >= 0 ) {
    {for(int k = 0; k < num_vecs; ++k ) {
      vecs[k]->acquireDetachedView( local_rng, &sub_vecs[k] );
      sub_vecs[k].setGlobalOffset( overlap_global_off );
    }}
    {for(int k = 0; k < num_targ_vecs; ++k ) {
      targ_vecs[k]->acquireDetachedView( local_rng, &sub_targ_vecs[k] );
      sub_targ_vecs[k].setGlobalOffset( overlap_global_off );
    }}
  }
  // Apply the RTOp operator object (all processors must participate)
  const bool validSubVecs = ( localSubDim_ > 0 && overlap_first_local_ele_off >= 0 );
  RTOpPack::SPMD_apply_op(
    locallyReplicated ? NULL : comm.get()                       // comm
    ,op                                                         // op
    ,num_vecs                                                   // num_vecs
    ,num_vecs && validSubVecs ? &sub_vecs[0] : NULL             // sub_vecs
    ,num_targ_vecs                                              // num_targ_vecs
    ,num_targ_vecs && validSubVecs ? &sub_targ_vecs[0] : NULL   // targ_sub_vecs
    ,reduct_obj                                                 // reduct_obj
    );
  // Free and commit the local data
  if( overlap_first_local_ele_off >= 0 ) {
    {for(int k = 0; k < num_vecs; ++k ) {
      sub_vecs[k].setGlobalOffset(local_rng.lbound());
      vecs[k]->releaseDetachedView( &sub_vecs[k] );
    }}
    {for(int k = 0; k < num_targ_vecs; ++k ) {
      sub_targ_vecs[k].setGlobalOffset(local_rng.lbound());
      targ_vecs[k]->commitDetachedView( &sub_targ_vecs[k] );
    }}
  }
  // Flag that we are leaving applyOp()
  in_applyOpImpl_ = false;
#ifdef THYRA_SPMD_VECTOR_BASE_DUMP
  if(show_dump) {
    *out << "\nLeaving SpmdVectorBase<double>::applyOp(...) ...\n";
  }
#endif // THYRA_SPMD_VECTOR_BASE_DUMP
}



#endif




void EpetraVector::acquireDetachedVectorViewImpl(
  const Teuchos::Range1D& rng_in,
  RTOpPack::ConstSubVectorView<double>* sub_vec) const
{
  /* if range is empty, return a null view */
  if (rng_in == Range1D::Invalid)
  {
    sub_vec->initialize(rng_in.lbound(), 0, Teuchos::ArrayRCP<double>(), 1);
    return ;
  }

  /* */
  const Range1D rng = validateRange(rng_in);

  /* If any requested elements are off-processor, fall back to slow 
   * VDB implementation */ 
  if ( 
    rng.lbound() < localOffset_ 
    || localOffset_ + localSubDim_ - 1 < rng.ubound())
  {
    VectorDefaultBase<double>::acquireDetachedVectorViewImpl(rng_in,sub_vec);
    return ;
  }

  /* All requested elements are on-processor. */
  const double* localValues = &(epetraVec_->operator[](0));
  Teuchos::ArrayRCP<const double> locVals(localValues, rng.lbound()-localOffset_,
    rng.size(), false);
//    rng.ubound()-localOffset_, false);

  OrdType stride = 1;
  
  sub_vec->initialize(rng.lbound(), rng.size(),
    locVals, stride);
//    localValues+(rng.lbound()-localOffset_), stride);
}




void EpetraVector::releaseDetachedVectorViewImpl(
  RTOpPack::ConstSubVectorView<double>* sub_vec) const
{
  if(
    sub_vec->globalOffset() < localOffset_ 
    || localOffset_+localSubDim_ < sub_vec->globalOffset()+sub_vec->subDim()
    )
  {
    // Let the default implementation handle it!
    VectorDefaultBase<double>::releaseDetachedVectorViewImpl(sub_vec);
    return;
  }
  // Nothing to deallocate!
  sub_vec->uninitialize();
}




void EpetraVector::acquireNonconstDetachedVectorViewImpl(
  const Teuchos::Range1D& rng_in, 
  RTOpPack::SubVectorView<double>* sub_vec) 
{
  if( rng_in == Range1D::Invalid ) {
    // Just return an null view
    sub_vec->initialize(rng_in.lbound(), 0, Teuchos::ArrayRCP<double>(), 1);
    return;
  }
  
  const Range1D rng = validateRange(rng_in);
  if(
    rng.lbound() < localOffset_ 
    ||
    localOffset_+localSubDim_-1 < rng.ubound()
    )
  {
    /* rng consists of off-processor elements so use the 
     * default implementation! */
    VectorDefaultBase<double>::acquireNonconstDetachedVectorViewImpl(rng_in,
      sub_vec);
    return;
  }

  /* All requested elements are on-processor. */
  double* localValues = &(epetraVec_->operator[](0));
  Teuchos::ArrayRCP<double> locVals(localValues, rng.lbound()-localOffset_,
//    rng.ubound()-localOffset_, false);
    rng.size(), false);
  OrdType stride = 1;

  sub_vec->initialize(rng.lbound(), rng.size(),
    locVals, stride);
//    localValues+(rng.lbound()-localOffset_),stride);
}



void EpetraVector::commitNonconstDetachedVectorViewImpl(
  RTOpPack::SubVectorView<double>* sub_vec
  )
{
  if(
    sub_vec->globalOffset() < localOffset_
    ||
    localOffset_+localSubDim_ < sub_vec->globalOffset()+sub_vec->subDim()
    )
  {
    // Let the default implementation handle it!
    VectorDefaultBase<double>::commitNonconstDetachedVectorViewImpl(sub_vec);
    return;
  }
  sub_vec->uninitialize();  // Nothing to deallocate!
}




Teuchos::Range1D EpetraVector::validateRange(const Teuchos::Range1D& rng_in) const 
{
  const Range1D rng = Teuchos::full_range(rng_in,0,globalDim_-1);
  TEST_FOR_EXCEPTION(
    !( 0 <= rng.lbound() && rng.ubound() < globalDim_ ), std::invalid_argument
    ,"EpetraVector::validateRange(rowRng): Error, the range rng = ["
    <<rng.lbound()<<","<<rng.ubound()<<"] is not "
    "in the range [0,"<<(globalDim_-1)<<"]!"
    );
  return rng;
}





double& EpetraVector::operator[](OrdType globalIndex) 
{
  const Epetra_BlockMap& myMap = epetraVec()->Map();
  return (*epetraVec())[myMap.LID(globalIndex)];
}

void EpetraVector::setElement(OrdType index, const double& value)
{
  int loc_index[1] = { index };
  epetraVec()->ReplaceGlobalValues(1, const_cast<double*>(&value), 
    loc_index);
}

void EpetraVector::addToElement(OrdType index, const double& value)
{
//  cout << "adding (" << index << ", " << value << ")" << std::endl;
  int loc_index[1] = { index };
  epetraVec()->SumIntoGlobalValues(1, const_cast<double*>(&value), 
    loc_index);
}

const double& EpetraVector::getElement(OrdType index) const 
{
  const Epetra_BlockMap& myMap = epetraVec()->Map();
  return (*epetraVec())[myMap.LID(index)];
}

void EpetraVector::getElements(const OrdType* globalIndices, int numElems,
  Teuchos::Array<double>& elems) const
{
  elems.resize(numElems);
  const Epetra_BlockMap& myMap = epetraVec()->Map();
  RCP<const Epetra_Vector> epv = epetraVec();

  for (int i=0; i<numElems; i++)
  {
    elems[i] = (*epv)[myMap.LID(globalIndices[i])];
  }
}

void EpetraVector::setElements(int numElems, const int* globalIndices,
  const double* values)
{
  Epetra_FEVector* vec = dynamic_cast<Epetra_FEVector*>(epetraVec().get());
  int ierr = vec->ReplaceGlobalValues(numElems, globalIndices, values);
  TEST_FOR_EXCEPTION(ierr < 0, std::runtime_error, "ReplaceGlobalValues returned "
    "ierr=" << ierr << " in EpetraVector::setElements()");
}

void EpetraVector::addToElements(int numElems, const int* globalIndices,
  const double* values)
{
  Epetra_FEVector* vec = dynamic_cast<Epetra_FEVector*>(epetraVec().get());
  int ierr = vec->SumIntoGlobalValues(numElems, globalIndices, values);
  TEST_FOR_EXCEPTION(ierr < 0, std::runtime_error, "SumIntoGlobalValues returned "
    "ierr=" << ierr << " in EpetraVector::addToElements()");
}

void EpetraVector::finalizeAssembly()
{
  Epetra_FEVector* vec = dynamic_cast<Epetra_FEVector*>(epetraVec().get());
  vec->GlobalAssemble();
}


void EpetraVector::print(std::ostream& os) const 
{
  epetraVec()->Print(os);
}


const Epetra_Vector& EpetraVector::getConcrete(const TSFExtended::Vector<double>& tsfVec)
{
  const EpetraVector* epv 
    = dynamic_cast<const EpetraVector*>(tsfVec.ptr().get());
  TEST_FOR_EXCEPTION(epv==0, std::runtime_error,
    "EpetraVector::getConcrete called on a vector that "
    "could not be cast to an EpetraVector");
  return *(epv->epetraVec());
}

Epetra_Vector& EpetraVector::getConcrete(TSFExtended::Vector<double>& tsfVec)
{
  EpetraVector* epv 
    = dynamic_cast<EpetraVector*>(tsfVec.ptr().get());
  TEST_FOR_EXCEPTION(epv==0, std::runtime_error,
    "EpetraVector::getConcrete called on a vector that "
    "could not be cast to an EpetraVector");
  return *(epv->epetraVec());
}


Epetra_Vector* EpetraVector::getConcretePtr(TSFExtended::Vector<double>& tsfVec)
{
  EpetraVector* epv 
    = dynamic_cast<EpetraVector*>(tsfVec.ptr().get());
  TEST_FOR_EXCEPTION(epv==0, std::runtime_error,
    "EpetraVector::getConcrete called on a vector that "
    "could not be cast to an EpetraVector");
  return epv->epetraVec().get();
}




