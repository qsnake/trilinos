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

#ifndef NLPIP_NLP_FIRST_ORDER_THYRA_MODEL_EVALUATOR_HPP
#define NLPIP_NLP_FIRST_ORDER_THYRA_MODEL_EVALUATOR_HPP

#include <vector>

#include "NLPInterfacePack_NLPThyraModelEvaluatorBase.hpp"

namespace NLPInterfacePack {

/** \brief Implement the %NLPFirstOrder interface using a
 * <tt>Thyra::ModelEvaluator</tt> object.
 *
 * ToDo: Finish documentation!
 */
class NLPFirstOrderThyraModelEvaluator
  : virtual public NLPFirstOrder
  , virtual public NLPThyraModelEvaluatorBase
{
public:

  /** \brief Initialize to uninitialized */
  NLPFirstOrderThyraModelEvaluator();

  /** \brief Calls <tt>initialize()</tt>. */
  NLPFirstOrderThyraModelEvaluator(
    const Teuchos::RCP<Thyra::ModelEvaluator<value_type> >  &model
    ,const int                                                      p_idx
    ,const int                                                      g_idx
    );

  /** \brief .Initialize given a <tt>Thyra::ModelEvaluator</tt> and
   * a description of how to interpret it.
   *
   * ToDo: Finish documentation!
   *
   * Todo: Add arguments for auxiliary inequalites and equalities
   */
  void initialize(
    const Teuchos::RCP<Thyra::ModelEvaluator<value_type> >  &model
    ,const int                                                      p_idx
    ,const int                                                      g_idx
    );

  /** @name Overridden public members from NLP */
  //@{

  /** \brief . */
  void initialize(bool test_setup);
  /** \brief . */
  void unset_quantities();

  //@}

  /** @name Overridden public members from NLPFirstOrder */
  //@{

  /** \brief Overridden to check the concrete type of Gc */
  void set_Gc(MatrixOp* Gc);
  /** \brief . */
  const NLPFirstOrder::mat_fcty_ptr_t factory_Gc() const;
  /** \brief Returns an ExampleBasisSystem */
  const basis_sys_ptr_t basis_sys() const;

  //@}

protected:

  /** @name Overridden protected members from NLPFirstOrder */
  //@{

  /** \brief . */
  void imp_calc_Gc(
    const Vector& x, bool newx
    ,const FirstOrderInfo& first_order_info) const;

  //@}

private:

  // /////////////////////////////////////////
  // Private member functions

  /** \brief . */
  void evalModel( 
    const Vector            &x
    ,bool                   newx
    ,const ZeroOrderInfo    *zero_order_info  // != NULL if only zero-order info
    ,const ObjGradInfo      *obj_grad_info    // != NULL if obj-grad and below info
    ,const FirstOrderInfo   *first_order_info // != NULL if first-order and below info
    ) const;

};	// end class NLPFirstOrderThyraModelEvaluator

}	// end namespace NLPInterfacePack

#endif	// NLPIP_NLP_FIRST_ORDER_THYRA_MODEL_EVALUATOR_HPP
