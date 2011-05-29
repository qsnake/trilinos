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

#include "MoochoPack_EvalNewPointTailoredApproachCoordinate_Step.hpp"
#include "ConstrainedOptPack_MatrixIdentConcatStd.hpp"
#include "NLPInterfacePack_NLPDirect.hpp"
#include "AbstractLinAlgPack_MatrixOp.hpp"
#include "AbstractLinAlgPack_MatrixZero.hpp"
#include "AbstractLinAlgPack_VectorMutable.hpp"
#include "Teuchos_dyn_cast.hpp"

namespace MoochoPack {

EvalNewPointTailoredApproachCoordinate_Step::EvalNewPointTailoredApproachCoordinate_Step(
  const deriv_tester_ptr_t& 	  deriv_tester
  ,const bounds_tester_ptr_t&	  bounds_tester
  ,EFDDerivTesting              fd_deriv_testing
  )
  :EvalNewPointTailoredApproach_Step(deriv_tester,bounds_tester,fd_deriv_testing)
{}

// protected

void EvalNewPointTailoredApproachCoordinate_Step::uninitialize_Y_Uy(
  MatrixOp         *Y
  ,MatrixOp        *Uy
  )
{
  // Nothing to free
}

void EvalNewPointTailoredApproachCoordinate_Step::calc_py_Y_Uy(
  const NLPDirect       &nlp
  ,const D_ptr_t        &D
  ,VectorMutable        *py
  ,MatrixOp             *Y
  ,MatrixOp             *Uy
  ,EJournalOutputLevel  olevel
  ,std::ostream         &out
  )
{
  namespace rcp = MemMngPack;
  using Teuchos::dyn_cast;

  MatrixIdentConcatStd
    &cY = dyn_cast<MatrixIdentConcatStd>(*Y);
  //
  // Y = [      I     ] space_xD  
  //     [    Zero    ] space_xI
  //        space_xD
  //
  VectorSpace::space_ptr_t
    space_x  = nlp.space_x(),
    space_xD = space_x->sub_space(nlp.var_dep())->clone(),
    space_xI = space_x->sub_space(nlp.var_indep())->clone();
  cY.initialize(
    space_x                                                // space_cols
    ,space_xD                                              // space_rows
    ,MatrixIdentConcatStd::BOTTOM                          // top_or_bottom
    ,1.0                                                   // alpha
    ,Teuchos::rcp(
      new MatrixZero(
        space_xI    // space_cols
        ,space_xD   // space_rows
        ) )                                            // D_ptr
    ,BLAS_Cpp::no_trans                                    // D_trans
    );
  // py is not altered here!
}

void EvalNewPointTailoredApproachCoordinate_Step::recalc_py(
  const MatrixOp          &D
  ,VectorMutable          *py
  ,EJournalOutputLevel    olevel
  ,std::ostream           &out
  )
{
  // py is not altered here!
}

void EvalNewPointTailoredApproachCoordinate_Step::print_calc_py_Y_Uy(
  std::ostream& out, const std::string& L
  ) const
{
  out
    << L << "*** Coordinate decomposition\n"
    << L << "py_k = py_k\n"
    << L << "Y = [ I ; 0 ] <: R^(n x m) [0 represented using MatrixZero]\n"
    << L << "Uy = Gc(var_dep,con_undecomp)\'\n"
    ;
}

}	// end namespace MoochoPack 
