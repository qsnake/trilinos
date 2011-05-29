#ifndef SUNDANCERIVARANODE_H
#define SUNDANCERIVARANODE_H

#include "SundanceDefs.hpp"
#include "SundancePoint.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_Array.hpp"


namespace Sundance
{
namespace Rivara
{
class Element;
class Edge;
using Sundance::Point;
using Teuchos::Array;

/**
 * Class Node is a vertex in a simplicial mesh.
 */
class Node
{
public:
  Node(int gid, const Point& x, int ownerProc, int label=-1);

  /**
   * Return the spatial position of this node
   */
  const Point& pt() const ;

  /**
   * Return the local index of this node
   */
  int localIndex() const {return localIndex_;}


  /**
   * Return the global index of this node
   */
  int globalIndex() const {return globalIndex_;}

  /**
   * Set the local index of this node
   */
  void setLocalIndex(int localIndex) {localIndex_ = localIndex;}

  /**
   * Add an element to the list of elements containing this node.
   */
  void addConnectingElement(Element* elem);

  /**
   * Add an edge to the list of edges containing this node
   */
  void addConnectingEdge(Edge* edge);

  /**
   * Return the rank of the proc that owns this node
   */
  int ownerProc() const {return ownerProc_;}


  /**
   * Return the label of this node
   */
  int label() const {return label_;}


private:
  int label_;
  int localIndex_;
  int globalIndex_;
  Point x_;

  Array<Element*> elements_;
  Array<Edge*> edges_;

  int ownerProc_;
};
}
}



#endif
