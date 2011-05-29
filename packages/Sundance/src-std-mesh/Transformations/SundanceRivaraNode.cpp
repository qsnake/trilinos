#include "SundanceRivaraNode.hpp"
#include "SundanceRivaraElement.hpp"
#include "SundanceRivaraEdge.hpp"
#include "SundanceRivaraMesh.hpp"

using namespace Sundance;
using namespace Sundance::Rivara;
using namespace Sundance;

Node::Node(int globalIndex, const Point& x, int ownerProc, int label)
  : label_(label), localIndex_(-1), globalIndex_(globalIndex), 
    x_(x), elements_(), edges_(), ownerProc_(ownerProc)
{}

void Node::addConnectingElement(Element* tri)
{
  elements_.append(tri);
}

void Node::addConnectingEdge(Edge* edge)
{
  edges_.append(edge);
}

const Point& Node::pt() const {return x_;}

