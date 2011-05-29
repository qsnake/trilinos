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

#include <ostream>
#include <typeinfo>
#include <iostream>

#include "AbstractLinAlgPack_MatrixSymDiagStd.hpp"
#include "AbstractLinAlgPack_MatrixSymOpNonsing.hpp"
#include "AbstractLinAlgPack_MatrixOpOut.hpp"
#include "AbstractLinAlgPack_MultiVectorMutable.hpp"
#include "AbstractLinAlgPack_VectorStdOps.hpp"
#include "AbstractLinAlgPack_VectorOut.hpp"
#include "AbstractLinAlgPack_VectorAuxiliaryOps.hpp"
#include "AbstractLinAlgPack_LinAlgOpPack.hpp"
#include "AbstractLinAlgPack_assert_print_nan_inf.hpp"
#include "ConstrainedOptPack_MatrixIdentConcat.hpp"
#include "IterationPack_print_algorithm_step.hpp"
#include "MoochoPack_IpState.hpp"
#include "MoochoPack_UpdateReducedSigma_Step.hpp"
#include "MoochoPack_moocho_algo_conversion.hpp"

#include "OptionsFromStreamPack_StringToIntMap.hpp"

#include "Teuchos_dyn_cast.hpp"

namespace MoochoPack {

UpdateReducedSigma_Step::UpdateReducedSigma_Step(
  const e_update_methods update_method
  )
  :
  update_method_(update_method)
  {}

bool UpdateReducedSigma_Step::do_step(
  Algorithm& _algo, poss_type step_poss, IterationPack::EDoStepType type
  ,poss_type assoc_step_poss
  )
  {
  using Teuchos::dyn_cast;
  using IterationPack::print_algorithm_step;

  NLPAlgo            &algo   = dyn_cast<NLPAlgo>(_algo);
  IpState             &s      = dyn_cast<IpState>(_algo.state());
  
  EJournalOutputLevel olevel  = algo.algo_cntr().journal_output_level();
  std::ostream        &out    = algo.track().journal_out();
  
  // print step header.
  if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) 
    {
    using IterationPack::print_algorithm_step;
    print_algorithm_step( _algo, step_poss, type, assoc_step_poss, out );
    }

  switch (update_method_)
    {
    case ALWAYS_EXPLICIT:
      {
      if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) 
        {
        out << "\nupdate_method is always_explicit, forming Reduced Sigma Explicitly ...\n";
        }
      FormReducedSigmaExplicitly(algo,s,olevel,out);
      }
      break;
    case BFGS_PRIMAL:
    case BFGS_DUAL_NO_CORRECTION:
    case BFGS_DUAL_EXPLICIT_CORRECTION:
    case BFGS_DUAL_SCALING_CORRECTION:
      {
      TEST_FOR_EXCEPTION(true, std::logic_error, "Option BFGS_primal not handled yet in UpdateReducedSigma\n");
      }
      break;
    default:
      TEST_FOR_EXCEPT(true); // local error ?
    };

  if( (int)olevel >= (int)PRINT_ITERATION_QUANTITIES ) 
    {
    out << "\nrHB_k =\n" << s.rHB().get_k(0);
    }

  return true;
  }


void UpdateReducedSigma_Step::print_step(
  const Algorithm& _algo, poss_type step_poss, IterationPack::EDoStepType type
  ,poss_type assoc_step_poss, std::ostream& out, const std::string& L
  ) const
  {
  out << L << "*** Update Z'*Sigma*Z\n"
    << L << "if (update_method == always_explicit) then\n"
    << L << "  Sigma_k = invXl*Vl-invXu*Vu\n"
    << L << "  Sigma_I = Sigma_k.sub_view(Z.I_rng)\n"
    << L << "  Sigma_D_sqrt = (Sigma_k.sub_view(Z.D_rng))^1/2\n"
    << L << "  J = Sigma_D_sqrt*Z\n"
    << L << "  rHB_k = J'*J + Sigma_I\n"
    << L << "elsif ( update_method == BFGS_??? ) then\n"
    << L << "  NOT IMPLEMENTED YET!\n"
    << L << "end\n";
  }

void UpdateReducedSigma_Step::FormReducedSigmaExplicitly(
  NLPAlgo& algo, IpState& s, EJournalOutputLevel olevel,  std::ostream& out
  )
  {
  namespace mmp = MemMngPack;
  using Teuchos::dyn_cast;
  using AbstractLinAlgPack::ele_wise_prod;
  using AbstractLinAlgPack::ele_wise_sqrt;
   using LinAlgOpPack::Mp_M;
   using LinAlgOpPack::Mp_MtM;
   using LinAlgOpPack::M_MtM;
   using LinAlgOpPack::M_StM;
  using LinAlgOpPack::V_MtV;
  using LinAlgOpPack::assign;

  // Calculate Reduced Sigma directly from
  // Sigma = invXl*Vl + invXu*Vu
  // Z_kT*Sigma*Z_k

  // Get the iteration quantities
  const MatrixIdentConcat     &Z     = dyn_cast<MatrixIdentConcat>(s.Z().get_k(0));
  const MatrixSymDiagStd  &invXl = s.invXl().get_k(0);
  const MatrixSymDiagStd  &invXu = s.invXu().get_k(0);
  const MatrixSymDiagStd  &Vl    = s.Vl().get_k(0);
  const MatrixSymDiagStd  &Vu    = s.Vu().get_k(0);
  
  MatrixSymDiagStd  &Sigma = s.Sigma().set_k(0);

  MatrixSymOpNonsing& rHB = dyn_cast<MatrixSymOpNonsing>(s.rHB().set_k(0));
  if (rHB.cols() != Z.cols())
    {
    // Initialize space in rHB
    dyn_cast<MatrixSymInitDiag>(rHB).init_identity(Z.space_rows(), 0.0);
    }
  
  // Calculate Sigma = invXl*Vl + invXu*Vu

  ele_wise_prod(1.0, invXl.diag(), Vl.diag(), &(Sigma.diag() = 0.0));
  
  if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS )
    out << "\n||Sigma_l||inf = " << Sigma.diag().norm_inf() << std::endl;
  if( (int)olevel >= (int)PRINT_VECTORS )
    out << "\nSigma_l =\n" << Sigma.diag();
  
  ele_wise_prod(1.0, invXu.diag(), Vu.diag(), &Sigma.diag() );

  if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS )
    out << "\n||Sigma_k||inf = ||Sigma_l + Sigma_u||inf = " << Sigma.diag().norm_inf() << std::endl;
  if( (int)olevel >= (int)PRINT_VECTORS ) 
    out << "\nSigma_k = Sigma_l + Sigma_u =\n" << Sigma.diag();

  // Calculate the cross term (Z'*Sigma*Ypy) first
  VectorSpace::vec_mut_ptr_t temp = Z.space_cols().create_member(0.0);
  ele_wise_prod(1.0, s.Ypy().get_k(0), Sigma.diag(), temp.get());
  VectorMutable& w_sigma = s.w_sigma().set_k(0);
  V_MtV(&w_sigma, Z, BLAS_Cpp::trans, *temp);

  if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS )
    out << "\n||w_sigma_k||inf = " << w_sigma.norm_inf() << std::endl;
  if( (int)olevel >= (int)PRINT_VECTORS ) 
    out << "\nw_sigma_k = \n" << w_sigma;
  
  // Calculate Reduced Sigma
  // Try sigma^1/2 making use of dependent and independent variables

  Teuchos::RCP<const VectorMutable>
    Sigma_D_diag = Sigma.diag().sub_view(Z.D_rng()),
    Sigma_I_diag = Sigma.diag().sub_view(Z.I_rng());
  const size_type
    Sigma_D_nz = Sigma_D_diag->nz(),
    Sigma_I_nz = Sigma_I_diag->nz();

  if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS )
    {
    out << "\nSigma_D.diag().nz() = " << Sigma_D_nz;
    out << "\nSigma_I.diag().nz() = " << Sigma_I_nz << std::endl;
    }

  if( Sigma_D_nz || Sigma_I_nz )
    {
    if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS )
      {
      out << "\nForming explicit, nonzero rHB_k = Z_k'*Sigma_k*Z_k ...\n";
      }
    if( Sigma_D_nz )
      {

      MatrixSymDiagStd Sigma_D_sqrt(Sigma_D_diag->clone());

      ele_wise_sqrt(&Sigma_D_sqrt.diag());

      if( static_cast<int>(olevel) >= static_cast<int>(PRINT_VECTORS) ) 
        {
        out << "\nSigma_D_sqrt =\n" << Sigma_D_sqrt;
        }
  
      Teuchos::RCP<MultiVectorMutable>
        J = Sigma_D_sqrt.space_cols().create_members(Z.cols());
      M_MtM(
        static_cast<MatrixOp*>(J.get())
        ,Sigma_D_sqrt, BLAS_Cpp::no_trans, Z.D(), BLAS_Cpp::no_trans);

      LinAlgOpPack::syrk( *J, BLAS_Cpp::trans, 1.0, 0.0, &rHB );

      if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ITERATION_QUANTITIES) ) 
        {
        out << "\nJ =\n" << *J;
        out << "\nJ'*J =\n" << rHB;
        }

      }

    if( Sigma_I_nz )
      {

      const MatrixSymDiagStd Sigma_I(
        Teuchos::rcp_const_cast<VectorMutable>(Sigma_I_diag)
        );

      if(Sigma_D_nz)
        {
        Mp_M( &rHB, Sigma_I, BLAS_Cpp::no_trans );
        }
      else
        {
        assign( &rHB, Sigma_I, BLAS_Cpp::no_trans );
        }

      }
    }
  else
    {
    if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS )
      {
      out << "\nSigma_k is zero so setting rHB_k = 0.0 ...\n";
      }
    rHB.zero_out();
    }

  /*
   // Try the full unspecialised calculation, this is expensive, but will
  // serve as a debug for the more efficient calculations.

  VectorSpace::multi_vec_mut_ptr_t Sigma_Z = Z.space_cols().create_members(Z.cols());
  M_MtM((MatrixOp*)Sigma_Z.get(), Sigma, BLAS_Cpp::no_trans, Z, BLAS_Cpp::no_trans);

  //std::cout << "Sigma_Z\n";
  //Sigma_Z->output(std::cout);

  VectorSpace::multi_vec_mut_ptr_t ZT_Sigma_Z = Z.space_rows().create_members(Z.cols());
  M_MtM((MatrixOp*)ZT_Sigma_Z.get(), (MatrixOp&)Z, BLAS_Cpp::trans, (MatrixOp&)*Sigma_Z, BLAS_Cpp::no_trans);

  std::cout << "ZT_Sigma_Z=\n";
  ZT_Sigma_Z->output(std::cout);
  */
  }


namespace {

const int local_num_options = 1;

enum local_EOptions 
  {
  UPDATE_METHOD
  };

const char* local_SOptions[local_num_options] = 
  {
  "update_method",
  };

const int num_update_methods = 5;

const char* s_update_methods[num_update_methods] = 
  {
  "always_explicit",
  "BFGS_primal",
  "BFGS_dual_no_correction",
  "BFGS_dual_explicit_correction",
  "BFGS_dual_scaling_correction"
  };

}

 
UpdateReducedSigma_StepSetOptions::UpdateReducedSigma_StepSetOptions(
  UpdateReducedSigma_Step* target
  , const char opt_grp_name[] )
  :
  OptionsFromStreamPack::SetOptionsFromStreamNode(
    opt_grp_name, local_num_options, local_SOptions ),
  OptionsFromStreamPack::SetOptionsToTargetBase< UpdateReducedSigma_Step >( target )
  {
  }

void UpdateReducedSigma_StepSetOptions::setOption( 
  int option_num, const std::string& option_value )
  {
  using OptionsFromStreamPack::StringToIntMap;

  typedef UpdateReducedSigma_Step target_t;
  switch( (local_EOptions)option_num ) 
    {
    case UPDATE_METHOD:
      {
      StringToIntMap	config_map( UpdateReducedSigma_opt_grp_name, num_update_methods, s_update_methods );
      target().update_method( (UpdateReducedSigma_Step::e_update_methods) config_map( option_value.c_str() ) );
      }
      break;
    default:
      TEST_FOR_EXCEPT(true);	// Local error only?
    }
  }

} // end namespace MoochoPack
