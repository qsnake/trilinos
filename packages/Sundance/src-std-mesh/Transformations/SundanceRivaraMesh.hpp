#ifndef SUNDANCERIVARAMESH_HPP
#define SUNDANCERIVARAMESH_HPP

#include "SundanceDefs.hpp"
#include "SundanceRivaraElement.hpp"
#include <stack>
#include "SundanceMap.hpp"
#include "SundanceIncrementallyCreatableMesh.hpp"
#include "SundanceRivaraElementIterator.hpp"


namespace Sundance
{
namespace Rivara
{
using Sundance::Map;


class RivaraMesh 
{
public:
  RivaraMesh(int dim, const MPIComm& comm);

  int addNode(const RCP<Node>& node);
  int addVertex(int globalIndex, const Point& x, int ownerProcID, int label);

  void addElement(const RCP<Element>& tri);
  int addElement(int globalIndex, const Array<int>& vertexGIDs, int ownerProc,
    int label);

  RCP<Edge> tryEdge(const RCP<Node>& a,
    const RCP<Node>& b,
    int& edgeSign);

  RCP<Face> tryFace(const RCP<Node>& a,
    const RCP<Node>& b,
    const RCP<Node>& c);

  const RCP<Face>& getFace(const RCP<Node>& a,
    const RCP<Node>& b,
    const RCP<Node>& c) const ;

  const RCP<Node>& node(int i) const {return nodes_[i];}

  int numNodes() const {return nodes_.length();}

  std::stack<Element*>& refinementSet()
    {return refinementSet_;}

  std::stack<double>& refinementAreas()
    {return refinementAreas_;}

  void refine();

  ElementIterator iterator() const ;

  friend class ElementIterator;

  RCP<Element> element(int i) const {return elements_[i];}

  int numRootElements() const {return elements_.length();}

  int numElements() const ;

  int spatialDim() const ;

  int& nextGID() {return nextGID_;}

  int nextGID() const {return nextGID_;}
private:

  int spatialDim_;
  
  int nextGID_;

  Array<RCP<Node> > nodes_;

  Array<RCP<Edge> > edges_;

  Array<RCP<Face> > faces_;

  Array<RCP<Element> > elements_;

  Array<Map<int, int> > nodeToEdgeMap_;

  Map<FaceNodes, int> faceToLIDMap_;

  std::stack<Element*> refinementSet_;

  std::stack<double> refinementAreas_;
};
}
}

#endif
