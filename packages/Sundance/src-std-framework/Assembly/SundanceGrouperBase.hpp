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

#ifndef SUNDANCE_GROUPERBASE_H
#define SUNDANCE_GROUPERBASE_H

#include "SundanceDefs.hpp"
#include "SundanceCellType.hpp"
#include "SundanceParametrizedCurve.hpp"
#include "SundanceMesh.hpp"
#include "Teuchos_RefCountPtr.hpp"

namespace Sundance
{
using namespace Teuchos;

class EquationSet;
class SparsitySuperset;
class MultiIndex;
class MultipleDeriv;
class QuadratureFamily;
class BasisFamily;
class IntegralGroup;


/** 
 * Grouper
 */
class GrouperBase
{
public:
  /** */
  GrouperBase() {}

  /** */
  virtual ~GrouperBase(){;}

  /** */
  virtual void findGroups(const EquationSet& eqn,
    const CellType& maxCellType,
    int spatialDim,
    const CellType& cellType,
    int cellDim,
    const QuadratureFamily& quad,
    const RCP<SparsitySuperset>& sparsity,
    bool isInternalBdry,
    Array<RCP<IntegralGroup> >& groups,
    const ParametrizedCurve& globalCurve,
    const Mesh& mesh ) const = 0 ;


  /** */
  void setVerbosity(
    int setupVerb,
    int integrationVerb,
    int transformVerb);
    

  /** */
  int setupVerb() const {return setupVerb_;}
    
  /** */
  int integrationVerb() const {return integrationVerb_;}
    
  /** */
  int transformVerb() const {return transformVerb_;}

protected:
  void extractWeakForm(const EquationSet& eqn,
    const MultipleDeriv& functionalDeriv,
    BasisFamily& testBasis, 
    BasisFamily& unkBasis,
    MultiIndex& miTest, MultiIndex& miUnk,
    int& rawVarID, int& rawUnkID,  
    int& reducedTestID, int& reducedUnkID, 
    int& testBlock, int& unkBlock, 
    int& rawParamID, int& reducedParamID,
    bool& isOneForm, bool& hasParam) const ;
                              
private:
  int setupVerb_;
  int integrationVerb_;
  int transformVerb_;
};

}


#endif
