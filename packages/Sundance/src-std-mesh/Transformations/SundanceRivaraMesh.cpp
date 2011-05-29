#include "SundanceRivaraMesh.hpp"
#include "SundanceRivaraEdge.hpp"
#include "SundanceRivaraElement.hpp"
#include "SundanceRivaraNode.hpp"
#include "Teuchos_StrUtils.hpp"
#include <fstream>
#include "Teuchos_TimeMonitor.hpp"
#include "SundanceOut.hpp"

using namespace Sundance::Rivara;
using namespace Sundance;
using Sundance::Map;

using namespace Teuchos;
using std::endl;




static Time& refKernelTimer() 
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("mesh refinement kernel"); 
  return *rtn;
}

RivaraMesh::RivaraMesh(int dim, const MPIComm& comm)
  : spatialDim_(dim), nextGID_(0),
    nodes_(), edges_(), elements_(), nodeToEdgeMap_()
{;}

void RivaraMesh::refine()
{
  TimeMonitor timer(refKernelTimer());
  while (refinementSet().size() > 0)
    {
      Element* next = refinementSet().top();
      refinementSet().pop();
      double maxArea = refinementAreas().top();
      refinementAreas().pop();
      next->refine(this, maxArea);
    }
}

int RivaraMesh::addNode(const RCP<Node>& node)
{
  int lid = nodes_.length();
  node->setLocalIndex(lid);
  nodes_.append(node);
  nodeToEdgeMap_.append(Map<int, int>());
  nextGID()++;
  return lid;
}


int RivaraMesh::addVertex(
  int globalIndex, const Point& x, 
  int ownerProcID, int label)
{
  RCP<Node> node = rcp(new Node(globalIndex, x, ownerProcID, label));
  return addNode(node);
}

void RivaraMesh::addElement(const RCP<Element>& tri)
{
  elements_.append(tri);
}


int RivaraMesh::addElement(
  int globalIndex, 
  const Array<int>& vertexGIDs, 
  int ownerProc,
  int label)
{
  int lid = elements_.size();
  RCP<Element> elem;

  switch(vertexGIDs.size())
  {
    case 3:
      elem = rcp(new Element(this, 
          nodes_[vertexGIDs[0]],
          nodes_[vertexGIDs[1]],
          nodes_[vertexGIDs[2]],
          ownerProc,label));
      break;
    case 4:
      elem = rcp(new Element(this, 
          nodes_[vertexGIDs[0]],
          nodes_[vertexGIDs[1]],
          nodes_[vertexGIDs[2]],
          nodes_[vertexGIDs[3]],
          ownerProc,label));
      break;
    default:
      TEST_FOR_EXCEPT(1);
  }
  elements_.append(elem);
  return lid;
  
}

RCP<Edge> RivaraMesh::tryEdge(const RCP<Node>& a,
                                          const RCP<Node>& b,
                                          int& edgeSign)
{
  int i = a->localIndex();
  int j = b->localIndex();

  if (nodeToEdgeMap_[i].containsKey(j))
    {
      int k = nodeToEdgeMap_[i].get(j);
      edgeSign = 1;
      return edges_[k];
    }
  else if (nodeToEdgeMap_[j].containsKey(i))
    {
      int k = nodeToEdgeMap_[j].get(i);
      edgeSign = -1;
      return edges_[k];
    }
  else
    {
      RCP<Edge> rtn = rcp(new Edge(a,b));
      edgeSign = 1;
      int k = edges_.length();
      edges_.append(rtn);
      nodeToEdgeMap_[i].put(j, k);
      return rtn;
    }
}


RCP<Face> RivaraMesh::tryFace(
  const RCP<Node>& a,
  const RCP<Node>& b,
  const RCP<Node>& c)
{
  FaceNodes f(a,b,c);

  int faceLID;
  if (faceToLIDMap_.containsKey(f))
  {
    faceLID = faceToLIDMap_[f];
  }
  else
  {
    faceLID = faces_.size();
    RCP<Face> newFace = rcp(new Face(a,b,c));
    faceToLIDMap_.put(newFace->nodes(), faceLID);
    faces_.append(newFace);
  }
  
  return faces_[faceLID];
}



const RCP<Face>& RivaraMesh::getFace(
  const RCP<Node>& a,
  const RCP<Node>& b,
  const RCP<Node>& c) const 
{
  FaceNodes f(a,b,c);
  return faces_[faceToLIDMap_.get(f)];
}



int RivaraMesh::numElements() const 
{
	int numLeaves = 0;
  for (int i=0; i<elements_.length(); i++)
    {
      numLeaves += elements_[i]->numLeaves();
    }
	return numLeaves;
}

int RivaraMesh::spatialDim() const
{
	return spatialDim_;
}


