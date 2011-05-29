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

#ifndef SUNDANCE_CELLPREDICATE_H
#define SUNDANCE_CELLPREDICATE_H



#include "SundanceDefs.hpp"
#include "SundanceMesh.hpp"
#include "SundanceCellPredicateBase.hpp"
#include "SundancePositionalCellPredicate.hpp"
#include "SundanceHandle.hpp"
#include "SundanceParametrizedCurve.hpp"
#include "SundanceCellCurvePredicate.hpp"

namespace Sundance
{
using namespace Teuchos;
  
  
/** 
 * User-level handle for predicates (deriving from CellPredicateBase)
 * used to decide whether
 * a given cell passes a CellFilter.
 */
class CellPredicate : public Sundance::Handle<CellPredicateBase>
{
public:
    
  /* handle boilerplate */
  HANDLE_CTORS(CellPredicate, CellPredicateBase);

  /** construct from a positional cell predicate functor */
  CellPredicate(const RCP<CellPredicateFunctorBase>& func);

  /** construct from a positional cell predicate functor */
  CellPredicate(Sundance::Handleable<CellPredicateFunctorBase>* func);

  /** construct from a positional cell predicate functor */
  CellPredicate(ParametrizedCurve& curve, CurveCellFilterMode filterMode);

  /** write to XML */
  XMLObject toXML() const {return ptr()->toXML();}

  /** */
  std::string description() const {return ptr()->description();}



  /** set the mesh on which cells are to be tested */
  void setMesh(const Mesh& mesh, int cellDim) const 
    {ptr()->setMesh(mesh, cellDim);}

  /** compare to another predicate, used for placement in STL containers */
  bool operator<(const CellPredicate& other) const ;


};

}

namespace std
{
inline ostream& operator<<(std::ostream& os, const Sundance::CellPredicate& pred)
{
  os << pred.toXML() << std::endl;
  return os;
}
}


#endif
