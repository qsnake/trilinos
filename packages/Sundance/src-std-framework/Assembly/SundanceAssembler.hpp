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

#ifndef SUNDANCE_ASSEMBLER_H
#define SUNDANCE_ASSEMBLER_H

#include "SundanceDefs.hpp"
#include "TSFLoadableVector.hpp"
#include "TSFLoadableMatrix.hpp"
#include "TSFLinearOperatorDecl.hpp"
#include "TSFVectorDecl.hpp"
#include "TSFVectorType.hpp"
#include "Teuchos_HashSet.hpp"
#include "Teuchos_ParameterList.hpp"
#include "TSFIncrementallyConfigurableMatrixFactory.hpp"
#include "TSFCollectivelyConfigurableMatrixFactory.hpp"
#include "TSFPartitionedMatrixFactory.hpp"
#include "TSFPartitionedToMonolithicConverter.hpp"
#include "SundanceRegionQuadCombo.hpp"
#include "SundanceMesh.hpp"
#include "SundanceEvalContext.hpp"
#include "SundanceIntegrationCellSpecifier.hpp"
#include "SundanceComputationType.hpp"


namespace Sundance
{
using namespace Teuchos;

class EquationSet;
class EvaluatableExpr;
class EvalManager;
class EvalVector;
class DiscreteSpace;
class DiscreteFunction;
class CellFilter;
class DOFMapBase;
class IntegralGroup;
class StdFwkEvalMediator;
class AssemblyKernelBase;


typedef std::set<int> ColSetType;

/** 
 * 
 */
class Assembler 
{
public:
  /** */
  Assembler(
    const Mesh& mesh, 
    const RCP<EquationSet>& eqn,
    const Array<VectorType<double> >& rowVectorType,
    const Array<VectorType<double> >& colVectorType,
    bool partitionBCs);


  /** */
  Assembler(
    const Mesh& mesh, 
    const RCP<EquationSet>& eqn);
      
  /** */
  const Array<RCP<DOFMapBase> >& rowMap() const 
    {return rowMap_;}

  /** */
  const Array<RCP<DOFMapBase> >& colMap() const 
    {return colMap_;}

  /** */
  const Array<RCP<DiscreteSpace> >& solutionSpace() const 
    {return externalColSpace_;}

  /** */
  const Array<RCP<DiscreteSpace> >& rowSpace() const 
    {return externalRowSpace_;}

  /** */
  VectorSpace<double> solnVecSpace() const ;

  /** */
  VectorSpace<double> rowVecSpace() const ;

  /** */
  const Array<RCP<Set<int> > >& bcRows() {return bcRows_;}

  /** Allocate, but do not fill, the matrix */
  TSFExtended::LinearOperator<double> allocateMatrix() const ;

  /** */
  void assemble(TSFExtended::LinearOperator<double>& A,
    Array<Vector<double> >& b) const ;

  /** */
  void assembleSensitivities(TSFExtended::LinearOperator<double>& A,
    Array<Vector<double> >& b) const ;


  /** */
  void assemble(Array<Vector<double> >& b) const ;

  /** */
  void evaluate(double& value,
    Array<Vector<double> >& gradient) const ;

  /** */
  void evaluate(double& value) const ;

  /** */
  static int& workSetSize() ;

      
  /** */
  void getGraph(int br, int bc,
    Array<int>& graphData,
    Array<int>& rowPtrs,
    Array<int>& nnzPerRow) const ;
      
  /** */
  void incrementalGetGraph(int br, int bc, 
    IncrementallyConfigurableMatrixFactory* mf) const ;

  /** */
  void flushConfiguration() 
    {
      numConfiguredColumns_ = 0;
      matNeedsConfiguration_ = true;
    }

  /** */
  Vector<double> convertToMonolithicVector(
    const Array<Vector<double> >& internalBlock,
    const Array<Vector<double> >& bcBlock) const ;

  /** */
  static int& numAssembleCalls() {static int rtn=0; return rtn;}

  /** */
  static bool& matrixEliminatesRepeatedCols() {static bool x = false; return x;}

  /** */
  const RCP<EquationSet>& eqnSet() const 
    {return eqn_;}

  /** */
  int maxWatchFlagSetting(const std::string& param) const ;

  /** */
  static Time& assemblyTimer() 
    {
      static RCP<Time> rtn 
        = TimeMonitor::getNewTimer("assembly"); 
      return *rtn;
    }

  /** */
  static Time& configTimer() 
    {
      static RCP<Time> rtn 
        = TimeMonitor::getNewTimer("matrix config"); 
      return *rtn;
    }
  
  /** */
  static Time& fillTimer() 
    {
      static RCP<Time> rtn 
        = TimeMonitor::getNewTimer("matrix/vector fill"); 
      return *rtn;
    }
  

private:

  /** */
  void init(const Mesh& mesh, const RCP<EquationSet>& eqn);

  /** */
  bool detectInternalBdry(int cellDim, const CellFilter& filter) const ;

  /** */
  void displayEvaluationResults(
    const EvalContext& context, 
    const EvaluatableExpr* evalExpr, 
    const Array<double>& constantCoeffs, 
    const Array<RCP<EvalVector> >& vectorCoeffs) const ;

  /** */
  void assemblyLoop(const ComputationType& compType,
    RCP<AssemblyKernelBase> kernel) const ;


  /** */
  void configureMatrix(LinearOperator<double>& A,
    Array<Vector<double> >& b) const ;

  /** */
  void configureVector(Array<Vector<double> >& b) const ;

  /** */
  void configureMatrixBlock(int br, int bc, 
    LinearOperator<double>& A) const ;

  /** */
  void configureVectorBlock(int br, Vector<double>& b) const ;

  /** */
  Array<Array<int> > findNonzeroBlocks() const ;

  /** */
  IntegrationCellSpecifier whetherToUseCofacets(
    const Array<RCP<IntegralGroup> >& groups,
    const EvaluatableExpr* ee,
    bool isMaximalCell,
    int verb) const ;

  
  
  /** */
  static int defaultWorkSetSize() {static int rtn=100; return rtn;}

  bool partitionBCs_;
      
  mutable bool matNeedsConfiguration_;
      
  mutable bool matNeedsFinalization_;

  mutable int numConfiguredColumns_;

  Mesh mesh_;

  RCP<EquationSet> eqn_;

  Array<RCP<DOFMapBase> > rowMap_;

  Array<RCP<DOFMapBase> > colMap_;

  Array<RCP<DiscreteSpace> > externalRowSpace_;

  Array<RCP<DiscreteSpace> > externalColSpace_;

  Array<RCP<DiscreteSpace> > privateRowSpace_;

  Array<RCP<DiscreteSpace> > privateColSpace_;

  Array<RCP<Set<int> > > bcRows_;

  Array<RCP<Set<int> > > bcCols_;

  Array<RegionQuadCombo> rqc_;

  Map<ComputationType, Array<EvalContext> > contexts_;

  Map<ComputationType, Array<int> > skipRqc_;

  Array<int> isBCRqc_;

  Array<int> isInternalBdry_;

  Map<ComputationType, Array<Array<RCP<IntegralGroup> > > > groups_;

  Array<RCP<StdFwkEvalMediator> > mediators_;

  Map<ComputationType, Array<const EvaluatableExpr*> > evalExprs_;

  RCP<EvalManager> evalMgr_;

  Array<RCP<Array<int> > > isBCRow_;

  Array<RCP<Array<int> > > isBCCol_;

  Array<RCP<std::set<int> > > remoteBCCols_;

  Array<int> lowestRow_;

  Array<int> lowestCol_;

  Array<VectorType<double> > rowVecType_;

  Array<VectorType<double> > colVecType_;

  Map<int, int> testIDToBlockMap_;

  Map<int, int> unkIDToBlockMap_;

  Map<int, int> fixedParamIDToVectorNumber_;

  Map<ComputationType, Array<IntegrationCellSpecifier> > rqcRequiresMaximalCofacets_;

  Array<RCP<PartitionedToMonolithicConverter> > converter_;

};

}



#endif
