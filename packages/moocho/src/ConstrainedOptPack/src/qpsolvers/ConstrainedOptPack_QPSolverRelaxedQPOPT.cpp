// @HEADER
// ***********************************************************************
// 
// Moocho: Multi-functional Object-Oriented arCHitecture for Optimization
//                  Copyright (2003) Sandia Corporation
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
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifdef CONSTRAINED_OPTIMIZATION_PACK_USE_QPOPT

#include <assert.h>

#include <algorithm>
#include <ostream>
#include <iomanip>

#include "ConstrainedOptPack_QPSolverRelaxedQPOPT.hpp"
#include "ConstrainedOptPack_QPOPT_CppDecl.hpp"
#include "AbstractLinAlgPack_MatrixOp.hpp"
#include "AbstractLinAlgPack_SpVectorClass.hpp"
#include "AbstractLinAlgPack_EtaVector.hpp"
#include "AbstractLinAlgPack_LinAlgOpPack.hpp"
#include "AbstractLinAlgPack_VectorMutableDense.hpp"
#include "DenseLinAlgPack_LinAlgOpPack.hpp"
#include "DenseLinAlgPack_DMatrixOut.hpp"
#include "DenseLinAlgPack_DVectorOut.hpp"

// //////////////////////////////////////////////////////////
// Local implementation functions.

namespace {

typedef FortranTypes::f_int       f_int;
typedef FortranTypes::f_dbl_prec  f_dbl_prec;
typedef FortranTypes::f_logical   f_logical;

// Compute:
//
// HESS * x = [ G  0 ] * [ X(1,N-1) ] = [ G * X(1,N-1) ]
//            [ 0  M ]   [   X(N)   ]   [ M * X(N)     ]
//
// The matrix vector product is implemented through the MatrixOp interface.
//
inline
void qphess_server_relax( const f_int& N, const f_int& LDH
  , const f_int& JTHCOL, const f_dbl_prec* H, const f_dbl_prec* X, f_dbl_prec* HX
  , f_int* IW, const f_int& LENIW, f_dbl_prec* W, const f_int& LENW )
{
  using DenseLinAlgPack::DVectorSlice;
  using AbstractLinAlgPack::SpVector;
  using LinAlgOpPack::V_MtV;
  using ConstrainedOptPack::QPSolverRelaxedQPOPT;

  // Here we have used some casting to pass on information about the qp solver
  // that called QPSOL.
  const QPSolverRelaxedQPOPT* qp_solver = reinterpret_cast<const QPSolverRelaxedQPOPT*>(H);
  const AbstractLinAlgPack::MatrixOp &G = *qp_solver->G();

  const DVectorSlice x(const_cast<f_dbl_prec*>(X),N);   // Just a view!
  DVectorSlice hx(HX,N);                                // Just a view!

  if( JTHCOL == 0 ) {
    // We are performing a dense mat-vec
    // hx(1,N-1) = G * x(1,N-1)
    AbstractLinAlgPack::VectorMutableDense x_n(x(1,N-1),Teuchos::null);
    AbstractLinAlgPack::VectorMutableDense hx_n(hx(1,N-1),Teuchos::null);
    V_MtV( &hx_n, G, BLAS_Cpp::no_trans, x_n );
    // hx(N) = bigM * x(N)
    hx(N) = qp_solver->use_as_bigM() * x(N);
  }
  else {
    // we are extracting the JTHCOL column of G so use sparse x
    if(JTHCOL == N) {
      // 0
      hx(1,N-1) = 0.0;
      // bigM
      hx(N) = qp_solver->use_as_bigM();
    }
    else {
      // G(:,JTHCOL)
      AbstractLinAlgPack::EtaVector e_j(JTHCOL,N-1);
      AbstractLinAlgPack::VectorMutableDense hx_n(hx(1,N-1),Teuchos::null);
      V_MtV( &hx_n, G, BLAS_Cpp::no_trans, e_j() );
      // 0
      hx(N) = 0.0;
    }
  }
}

}	// end namespace

// ///////////////////////////////////////////////////////////////////////////
// Fortran declarations.

extern "C" {

// These are declarations for the subroutines that preform the communication
// between C++ and Fortran.  There is no use in putting them in a
// namespace since the namespace name will not be used by the linker since
// we are using extern "C".

//
FORTRAN_FUNC_DECL_UL_( void, QPHESS_SERVER_RELAX2, qphess_server_relax2 ) ( const f_int& N, const f_int& LDH
  , const f_int& JTHCOL, const f_dbl_prec* H, const f_dbl_prec* X, f_dbl_prec* HX
  , f_int* IW, const f_int& LENIW, f_dbl_prec* W, const f_int& LENW )
{
  qphess_server_relax( N, LDH, JTHCOL, H, X, HX, IW, LENIW, W, LENW );
}

}	// end extern "C"

namespace LinAlgOpPack {
  using AbstractLinAlgPack::Mp_StM;
  using AbstractLinAlgPack::Vp_StMtV;
}

// ///////////////////////////////////////
// QPSolverRelaxedQPOPT

namespace ConstrainedOptPack {

QPSolverRelaxedQPOPT::QPSolverRelaxedQPOPT(
  value_type        max_qp_iter_frac
  )
  :max_qp_iter_frac_(max_qp_iter_frac)
  ,ITMAX_(100)
  ,FEATOL_(1.0e-10)
  ,LDH_(1)
{}

QPSolverRelaxedQPOPT::~QPSolverRelaxedQPOPT()
{
  release_memory();
}

// Overridden from QPSolverRelaxed

void QPSolverRelaxedQPOPT::release_memory()
{
  // ToDo: Resize all arrays to zero!
  QPSolverRelaxedQPOPTSOL::release_memory();
}

// Overridden protected members

QPSolverRelaxedQPOPT::f_int QPSolverRelaxedQPOPT::liwork(f_int N, f_int NCLIN) const
{	return 2* N + 3; }

QPSolverRelaxedQPOPT::f_int QPSolverRelaxedQPOPT::lrwork(f_int N, f_int NCLIN) const
{	return NCLIN > 0 ? 2*N*N + 8*N + 5*NCLIN : N*N + 8 *N; }

QPSolverRelaxedQPOPTSOL::EInform QPSolverRelaxedQPOPT::call_qp_solver(bool warm_start)
{
  
  // Set the rest of the QPOPT inputs that could not be set in the constructor.

  BIGBND_ = this->infinite_bound();
  ITMAX_ = max_qp_iter_frac() * N_;
  LDA_	= ( A_.rows() > 0 ? A_.rows() : 1 );
  H_		= reinterpret_cast<f_dbl_prec*>(this);	// used to implement QPHESS
  AX_.resize( NCLIN_ > 0 ? NCLIN_ : 1 );

  // Set option parameters
  {
    namespace ns = QPOPT_CppDecl;
    namespace ft = FortranTypes;
    ns::reset_defaults();
    ns::set_logical_option( ns::WARM_START            , warm_start ? ft::F_FALSE : ft::F_TRUE		);
    ns::set_real_option(    ns::FEASIBILITY_TOLERANCE , FEATOL_                                 );
    ns::set_real_option(    ns::INFINITE_BOUND_SIZE   , BIGBND_                                 );
    ns::set_int_option(     ns::ITERATION_LIMIT       , ITMAX_                                  );
    ns::set_int_option(     ns::PRINT_FILE            , 0                                       );
    ns::set_int_option(     ns::SUMMARY_FILE          , 0                                       );
    ns::set_int_option(     ns::PRINT_LEVEL           , 0                                       );
    ns::set_int_option(     ns::PROBLEM_TYPE          , ns::QP2                                 );
  }

  QPOPT_CppDecl::qpopt(
    N_, NCLIN_, LDA_, LDH_, NCLIN_ ? &A_(1,1) : NULL, &BL_(1), &BU_(1)
    , &CVEC_(1), H_
    , FORTRAN_NAME_UL_(QPHESS_SERVER_RELAX2,qphess_server_relax2)
    , &ISTATE_[0], &X_(1), INFORM_, ITER_, OBJ_, &AX_(1)
    , &CLAMDA_(1), &IWORK_[0], LIWORK_, &WORK_[0], LWORK_ );

  EInform return_inform;
  typedef QPSolverRelaxedQPOPTSOL bc;
  switch(INFORM_) {
    case STRONG_LOCAL_MIN:
      return_inform = bc::STRONG_LOCAL_MIN;
      break;
    case WEAK_LOCAL_MIN:
      return_inform = bc::WEAK_LOCAL_MIN;
      break;
    case UNBOUNDED:
      TEST_FOR_EXCEPTION(
        true, Unbounded
        ,"QPSolverRelaxedQPOPT::call_qp_solver() : Error,"
        " QPOPT returned that the QP is unbounded!" );
    case INFEASIBLE:
      TEST_FOR_EXCEPTION(
        true, Infeasible
        ,"QPSolverRelaxedQPOPT::call_qp_solver() : Error,"
        " QPOPT returned that the QP is infeasible!" );
    case ITMAX_EXCEEDED:
      return_inform = bc::MAX_ITER_EXCEEDED;
      break;
    case MAX_DOF_TOO_SMALL:
      TEST_FOR_EXCEPTION(
        true, std::runtime_error,
        "QPSolverRelaxedQPOPT::call_qp_solver() : Error,"
        " QPOPT says that the max dof is too small" );
    case INVALID_INPUT:
      TEST_FOR_EXCEPTION(
        true, InvalidInput,
        "QPSolverRelaxedQPOPT::call_qp_solver() : Error,"
        " QPOPT returned that the input is invalid" );
    case PROB_TYPE_NOT_REGOG:
      TEST_FOR_EXCEPTION(
        true, std::logic_error,
        "QPSolverRelaxedQPOPT::call_qp_solver() : Error,"
        " QPOPT says that the problem type is not recognized" );
      break;
    default:
      TEST_FOR_EXCEPT(true); // Should not happen
  }
  
  return return_inform;
}

} // end namespace ConstrainedOptPack

#endif // CONSTRAINED_OPTIMIZATION_PACK_USE_QPOPT
