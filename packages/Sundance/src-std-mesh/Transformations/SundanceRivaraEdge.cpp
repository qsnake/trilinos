#include "SundanceRivaraEdge.hpp"
#include "SundanceRivaraElement.hpp"
#include "SundanceRivaraNode.hpp"
#include "SundanceRivaraMesh.hpp"

using namespace Sundance::Rivara;
using namespace Teuchos;



Edge::Edge(const RCP<Node>& a,
           const RCP<Node>& b)
  : label_(-1),nodes_(tuple(a,b)), elements_(), midpoint_(),
		ownerProc_()
{
	if (a->ownerProc() > b->ownerProc())
		{
			ownerProc_ = a->ownerProc();
		}
	else
		{
			ownerProc_ = b->ownerProc();
		}
}

void Edge::addConnectingElement(Element* tri)
{
  elements_.append(tri);
}

double Edge::length() const 
{
  const Point& x1 = nodes_[0]->pt();
  const Point& x2 = nodes_[1]->pt();

  return sqrt((x1-x2)*(x1-x2));
}

void Edge::getUnrefinedCofacets(Array<Element*>& c) const 
{
  for (int i=0; i<elements_.length(); i++)
    {
      if (!elements_[i]->hasChildren()) c.append(elements_[i]);
    }
}

RCP<Node> Edge::bisect(RivaraMesh* mesh)
{
  /* if we've already been bisected, return the existing midpoint node */
  if (!(midpoint_.get() == 0))
    {
      return midpoint_;
    }

  const Point& x1 = nodes_[0]->pt();
  const Point& x2 = nodes_[1]->pt();

  int nextGID = mesh->nextGID();
  midpoint_ = rcp(new Node(nextGID, 0.5*(x1 + x2), ownerProc_));
  mesh->addNode(midpoint_);

  int s;
  RCP<Edge> sub1 = mesh->tryEdge(nodes_[0], midpoint_, s);
  RCP<Edge> sub2 = mesh->tryEdge(midpoint_, nodes_[1], s);
  sub1->setParent(this);
  sub2->setParent(this);

  sub1->setLabel(label_);
  sub2->setLabel(label_);

  setChildren(sub1.get(), sub2.get());

  return midpoint_;
}

