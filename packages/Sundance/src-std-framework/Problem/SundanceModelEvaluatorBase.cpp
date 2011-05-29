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




#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"
#include "SundanceModelEvaluatorBase.hpp"

#ifdef HAVE_SUNDANCE_MOOCHO
#include "Thyra_DefaultSpmdVectorSpace.hpp"


SundanceModelEvaluator
::SundanceModelEvaluator(const VectorType<double>& vecType)
  : vecType_(vecType), 
    objectiveSpace_(rcp(new Thyra::DefaultSpmdVectorSpace<double>(1)))
{;}



ModelEvaluatorBase::InArgs<double> SundanceModelEvaluator::createInArgs() const
{
  ModelEvaluatorBase::InArgsSetup<double> inArgs;
  
  inArgs.setModelEvalDescription(description());
  /* Note: Np is the number of parameter *blocks*, not the number of parameters */
  inArgs.set_Np(1);
  inArgs.setSupports(IN_ARG_x, true);

  return inArgs;
}




ModelEvaluatorBase::OutArgs<double> SundanceModelEvaluator::createOutArgsImpl() const
{
  ModelEvaluatorBase::OutArgsSetup<double> outArgs;

  outArgs.setModelEvalDescription(this->description());
  outArgs.set_Np_Ng(1,1);
  outArgs.setSupports(OUT_ARG_f,true);
  outArgs.setSupports(OUT_ARG_W_op,true);
  outArgs.set_W_properties(
    DerivativeProperties(
      DERIV_LINEARITY_NONCONST
      ,DERIV_RANK_FULL
      ,true // supportsAdjoint
      )
    );
  outArgs.setSupports(OUT_ARG_DfDp,0,DERIV_MV_BY_COL);
  outArgs.set_DfDp_properties(
    0,DerivativeProperties(
      DERIV_LINEARITY_NONCONST
      ,DERIV_RANK_DEFICIENT
      ,true // supportsAdjoint
      )
    );
  outArgs.setSupports(OUT_ARG_DgDx,0,DERIV_TRANS_MV_BY_ROW);
  outArgs.set_DgDx_properties(
    0,DerivativeProperties(
      DERIV_LINEARITY_NONCONST
      ,DERIV_RANK_DEFICIENT
      ,true // supportsAdjoint
      )
    );
  outArgs.setSupports(OUT_ARG_DgDp,0,0,DERIV_TRANS_MV_BY_ROW);
  outArgs.set_DgDp_properties(
    0,0,DerivativeProperties(
      DERIV_LINEARITY_NONCONST
      ,DERIV_RANK_DEFICIENT
      ,true // supportsAdjoint
      )
    );
  return outArgs;
}

ModelEvaluatorBase::InArgs<double>
SundanceModelEvaluator::getNominalValues() const
{
  typedef ModelEvaluatorBase MEB;
  MEB::InArgs<double> nominalValues = this->createInArgs();
  nominalValues.set_x(getInitialState().ptr());
  nominalValues.set_p(0,getInitialParameters().ptr());
  return nominalValues;
}

void SundanceModelEvaluator
::evalModelImpl(const ModelEvaluatorBase::InArgs<double>& inArgs,
            const ModelEvaluatorBase::OutArgs<double>& outArgs) const
{
  Tabs tabs;
  if (MPIComm::world().getRank()==0) SUNDANCE_VERB_MEDIUM(tabs << 
                       "------------------calling eval model -----------------------");
  /* read input args. The const casts are needed until ConstVector is ready */
  TSFExtended::Vector<double> x = rcp_const_cast<VectorBase<double> >(inArgs.get_x());

  TSFExtended::Vector<double> p = rcp_const_cast<VectorBase<double> >(inArgs.get_p(0));


  /* Get objects into which the output value will be written */

  /* constraint residual */
  TSFExtended::Vector<double> f = outArgs.get_f();

  /* objective function value */
  TSFExtended::Vector<double> g = outArgs.get_g(0);

  /* df/dx */
  TSFExtended::LinearOperator<double> df_dx = 
    rcp_dynamic_cast<LinearOpBase<double> >(outArgs.get_W_op());
  
  if (outArgs.get_W_op().get()!=0) 
    {
      TEST_FOR_EXCEPTION(df_dx.ptr().get()==0, RuntimeError,  
                         "W_op is non-null but could not be cast to "
                         "a SingleScalarTypeOpBase<double>");
    }
  

  /* df/dp */
  RCP<MultiVectorBase<double> > df_dp_mv 
    = outArgs.get_DfDp(0).getDerivativeMultiVector().getMultiVector();
  Array<Vector<double> > df_dp(get_p_space(0)->dim());
  if (df_dp_mv.get() != 0)
    {
      for (int i=0; i<df_dp.size(); i++)
        {
          df_dp[i] = df_dp_mv->col(0);
        }
    }

  

  /* dg/dp */
  RCP<MultiVectorBase<double> > dg_dp_trans_mv 
    = outArgs.get_DgDp(0,0).getDerivativeMultiVector().getMultiVector();
  Vector<double> dg_dp_T;
  if (dg_dp_trans_mv.get() != 0)
    {
      dg_dp_T = dg_dp_trans_mv->col(0);
    }

  /* dg/dx */
  RCP<MultiVectorBase<double> > dg_dx_trans_mv 
    = outArgs.get_DgDx(0).getDerivativeMultiVector().getMultiVector();
  Vector<double> dg_dx_T;
  if (dg_dx_trans_mv.get() != 0)
    {
      dg_dx_T = dg_dx_trans_mv->col(0);
    }


  /* do the evaluation */
  double gVal;
  internalEvalModel(x, p, f, gVal, df_dx, df_dp, dg_dp_T, dg_dx_T);

  if (outArgs.get_W_op().get()!=0) 
    {
      TEST_FOR_EXCEPTION(df_dx.ptr().get()==0, RuntimeError,  
                         "W_op is non-null but could not be cast to a "
                         "SingleScalarTypeOpBase<double>");
    }

  /* Fill in objective function value */
  if ( g.ptr().get() != 0 )
    {
      Thyra::set_ele(0,gVal,g.ptr().get());
      // bvbw g[0] = gVal;
    }
  /* f, g, and df_dx are already in the right form. 
   * We need to fill in the multivectors df_dp, dg_dp, and dg_dx */
  
  /* if requested, copy the vectors df_dp into the columns of df_dp_mv */
  if (df_dp_mv.get() != 0)
    {
      Tabs tab2;
      if (MPIComm::world().getRank()==0) SUNDANCE_VERB_MEDIUM(tab2 << "df_dp multivector was requested");
      for (int i=0; i<df_dp.size(); i++)
        {
          Vector<double> col = df_dp_mv->col(i);
          col.acceptCopyOf(df_dp[i]);
        }
    }
  else
    {
      Tabs tab2;
      if (MPIComm::world().getRank()==0) SUNDANCE_VERB_MEDIUM(tab2 << "df_dp multivector was NOT requested");
    }

  /* 
   * If requested, set the dg_dx return value 
   */
  if (dg_dx_trans_mv.get() != 0)
    {
      Tabs tab2;
      if (MPIComm::world().getRank()==0) SUNDANCE_VERB_MEDIUM(tab2 << "dg_dx multivector was requested");
      Vector<double> dg_dx_mv_row =  dg_dx_trans_mv->col(0);
      if (MPIComm::world().getRank()==0) SUNDANCE_VERB_HIGH(tab2 << "extracted zeroth column: " 
                         << dg_dx_mv_row.description());
      TEST_FOR_EXCEPT(dg_dx_T.ptr().get()==0);
      if (MPIComm::world().getRank()==0) SUNDANCE_VERB_MEDIUM(tab2 << "output from internal eval: " 
                           << dg_dx_T.description());
      dg_dx_mv_row.acceptCopyOf(dg_dx_T);
      if (MPIComm::world().getRank()==0) SUNDANCE_VERB_MEDIUM(tab2 << "got dg_dx multivector");
    }
  else
    {
      Tabs tab2;
      if (MPIComm::world().getRank()==0) SUNDANCE_VERB_MEDIUM(tab2 << "dg_dx multivector was NOT requested");
    }

  /* 
   * If requested, set the dg_dp return value 
   */
  if (dg_dp_trans_mv.get() != 0)
    {
      Tabs tab2;
      if (MPIComm::world().getRank()==0) SUNDANCE_VERB_MEDIUM(tab2 << "dg_dp multivector was requested");
      Vector<double> dg_dp_mv_row =  dg_dp_trans_mv->col(0);
      dg_dp_mv_row.acceptCopyOf(dg_dp_T);
      if (MPIComm::world().getRank()==0) SUNDANCE_VERB_MEDIUM(tab2 << "got dg_dp multivector");

    }
  else
    {
      Tabs tab2;
      if (MPIComm::world().getRank()==0) SUNDANCE_VERB_MEDIUM(tab2 << "dg_dp multivector was NOT requested");
    }

  if (MPIComm::world().getRank()==0) SUNDANCE_VERB_MEDIUM(tabs << "done evaluation!");
}



#endif

