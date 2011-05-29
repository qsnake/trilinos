#ifndef SUNDANCERIVARAFACE_H
#define SUNDANCERIVARAFACE_H

#include "SundanceDefs.hpp"
#include "SundanceRivaraTreeNode.hpp"
#include "SundanceRivaraNode.hpp"
#include "SundanceSet.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_RefCountPtr.hpp"

using namespace Teuchos;
namespace Sundance
{
namespace Rivara
{
class Element;
class RivaraMesh;

using Sundance::Set;
using Sundance::makeSet;
using Teuchos::Array;
using Teuchos::RefCountPtr;

/** */
class FaceNodes
{
public:
  /** */
  FaceNodes(
    const RCP<Node>& a, 
    const RCP<Node>& b,
    const RCP<Node>& c)
    : nodeLIDSet_(rcp(new Set<int>(makeSet(a->localIndex(), b->localIndex(), c->localIndex()))))
    {}


  /** */
  bool operator<(const FaceNodes& other) const
    {
      return (*nodeLIDSet_) < (*(other.nodeLIDSet_));
    }

  /** */
  const Set<int>& nodes() const {return *nodeLIDSet_;}

private:
  RCP<Set<int> > nodeLIDSet_;
};

/**
 * class Face is a two-dimensional face in a simplicial mesh.
 */

class Face
{
public:
  /**
   * Construct with three nodes
   */
  Face(
    const RCP<Node>& a, 
    const RCP<Node>& b,
    const RCP<Node>& c
    )
    : label_(-1), 
      globalIndex_(-1), 
      nodes_(a,b,c), 
      ownerProc_(std::max(std::max(a->ownerProc(), b->ownerProc()), c->ownerProc()))
    {}

  /** */
  const FaceNodes& nodes() const {return nodes_;}

  /** */
  int ownerProc() const {return ownerProc_;}

  /**
   * Return the global index of this edge
   */
  int globalIndex() const {return globalIndex_;}

  /**
   * Set the global index of this edge
   */
  void setGlobalIndex(int globalIndex) {globalIndex_ = globalIndex;}

  /**
   * Set the label of this edge
   */
  void setLabel(int label) {label_=label;}

  /** 
   * Get the label
   */
  int label() const {return label_;}
private:

  int label_;
  int globalIndex_;
  FaceNodes nodes_;
  int ownerProc_;
};
}
}

#endif
