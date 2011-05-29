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

#ifndef CONSTRAINED_OPTIMIZATION_PACK_TYPES_H
#define CONSTRAINED_OPTIMIZATION_PACK_TYPES_H

#include "NLPInterfacePack_Types.hpp"
#include "NLPInterfacePack_NLP.hpp"

namespace ConstrainedOptPack {

#include "NLPInterfacePack_PublicTypes.ud"

/// Bounds type
enum EBounds { FREE, UPPER, LOWER, EQUALITY };

// concrete classes

class VariableBoundsTester;

// abstract classes

class MatrixSymAddDelUpdateableWithOpFactorized;
class MatrixIdentConcat;
class MeritFuncCalc1D;
class MeritFuncCalc;
class MeritFuncNLP;
class MeritFuncNLE;
class MeritFuncNLF;
class MeritFuncNLPDirecDeriv;
class MeritFuncPenaltyParam;
class MeritFuncPenaltyParams;
class DirectLineSearch_Strategy;

// concrete subclasses

class MeritFuncCalc1DQuadratic;
class MeritFuncCalcNLP;
class MeritFuncNLPL1;
class MeritFuncNLPModL1;
//class MeritFuncCalcNLE;
//class MeritFuncCalcNLF;
//class MatrixHessianSuperBasic;
//class MatrixHessianSuperBasicInitDiagonal;
//class MatrixSymPosDefInvCholFactor;
class MatrixSymPosDefLBFGS;
class MatrixSymAddDelBunchKaufman;
class MatrixSymHessianRelaxNonSing;
class MatrixIdentConcatStd;
class DirectLineSearchArmQuad_Strategy;
class DirectLineSearchArmQuad_StrategySetOptions;
class VarReductOrthogDenseStd_Strategy;

// decomposition classes

class DecompositionSystem;
class DecompositionSystemVarReduct;
class DecompositionSystemVarReductPerm;
class DecompositionSystemVarReductPermStd;
class DecompositionSystemVarReductImp;
class DecompositionSystemCoordinate;
class DecompositionSystemOrthogonal;
class DecompositionSystemTester;
class DecompositionSystemTesterSetOptions;

// Abstract QP solvers

class QPSolverRelaxed;
class QPSolverRelaxedTester;
class QPSolverRelaxedTesterSetOptions;

// Concrete QP solvers

//class QPSchur;
//class QPSolverRelaxedQPSchurRangeSpace;

}	// end namespace ConstrainedOptPack 

#endif // CONSTRAINED_OPTIMIZATION_PACK_TYPES_H
