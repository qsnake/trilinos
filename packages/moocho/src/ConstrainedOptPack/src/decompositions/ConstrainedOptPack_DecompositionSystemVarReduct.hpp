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

#ifndef DECOMPOSITION_SYSTEM_VAR_REDUCT_H
#define DECOMPOSITION_SYSTEM_VAR_REDUCT_H

#include "ConstrainedOptPack_DecompositionSystem.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"
#include "AbstractLinAlgPack_BasisSystem.hpp"

namespace ConstrainedOptPack {

/** \brief Specialization of \c DecompositionSystem for variable reduction decompositions.
 *
 * This interface abstracts a variable reduction decomposition where:
 *
 \verbatim
  
  Gc' = [ C  N ] 
        [ E  F ]

  Z   = [ D ]
        [ I ]

  Uz  = F + E * D

      where:
           C = Gc(var_dep,con_decomp)'     [nonsingular]
           N = Gc(var_indep,con_decomp)'
           E = Gc(var_dep,con_undecomp)'
           F = Gc(var_indep,con_undecomp)'
           D = -inv(C) * N
 \endverbatim
 *
 * This interface simply allows clients to determine how \c D and \c Uz
 * are implemented (implicitly or explicity).
 */
class DecompositionSystemVarReduct : public DecompositionSystem {
public:

  /** @name Public types */
  //@{

  /** \brief . */
  enum EExplicitImplicit {
    MAT_IMP_EXPLICIT
    ,MAT_IMP_IMPLICIT
    ,MAT_IMP_AUTO
  };

  //@}

  /** @name Matrix representations */
  //@{

  /// Set whether to use explicit or implicit <tt>D = -inv(C)*N</tt> matrix.
  STANDARD_MEMBER_COMPOSITION_MEMBERS(EExplicitImplicit,D_imp);
  /// Set whether to use explicit or implicit <tt>Uz = F + E * D</tt> matrix.
  STANDARD_MEMBER_COMPOSITION_MEMBERS(EExplicitImplicit,Uz_imp);

    // ToDo: The above could be implemented as pure virtual funtions if needed later!

  //@}
    
  /** @name Constructors / initializers */
  //@{

  /** \brief . */
  DecompositionSystemVarReduct(
    EExplicitImplicit     D_imp    = MAT_IMP_AUTO
    ,EExplicitImplicit    Uz_imp   = MAT_IMP_AUTO
    )
    :D_imp_(D_imp), Uz_imp_(Uz_imp)
  {}

  //@}

  /** @name Variable partitions. */
  //@{

  /** \brief . */
  virtual Range1D var_indep() const = 0;
  /** \brief . */
  virtual Range1D var_dep() const = 0;

  //@}

private:

  // not defined and not to be called!
  DecompositionSystemVarReduct(const DecompositionSystemVarReduct&);
  DecompositionSystemVarReduct& operator=(const DecompositionSystemVarReduct&);

};	// end class DecompositionSystemVarReduct

}	// end namespace ConstrainedOptPack

#endif	// DECOMPOSITION_SYSTEM_VAR_REDUCT_H
