#ifndef SUNDANCE_RIVARADRIVER_HPP
#define SUNDANCE_RIVARADRIVER_HPP


#include "SundanceMeshType.hpp"
#include "SundanceMeshTransformationBase.hpp"
#include "SundanceRivaraMesh.hpp"
#include "SundanceExpr.hpp"

namespace Sundance
{

class Mesh;
using Sundance::Expr;

class RefinementTransformation : public MeshTransformationBase
{
public:
  /** */
  RefinementTransformation(const MeshType& meshType, const Expr& errExpr,
    const double& reqErr, const double& minArea)
    : MeshTransformationBase(meshType), meshType_(meshType),
      errExpr_(errExpr), reqErr_(reqErr), minArea_(minArea),
      numRefined_(-1) {}

  /** */
  Mesh apply(const Mesh& inputMesh) const ;

  /** */
  int numRefined() const {return numRefined_;}

  /* */
  GET_RCP(MeshTransformationBase);
private:
  /** */
  void meshToRivara(
    const Mesh& mesh,
    Array<int>& lidMap,
    RCP<Rivara::RivaraMesh>& rivMesh) const ;

  /** */
  Mesh rivaraToMesh(const RCP<Rivara::RivaraMesh>& rivMesh,
    const MPIComm& comm) const ;

  MeshType meshType_;
  Expr errExpr_;
  double reqErr_;
  double minArea_;
  mutable int numRefined_;
};


}


#endif
