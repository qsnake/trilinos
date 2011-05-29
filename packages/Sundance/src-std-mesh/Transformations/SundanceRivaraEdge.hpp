#ifndef SUNDANCERIVARAEDGE_H
#define SUNDANCERIVARAEDGE_H

#include "SundanceDefs.hpp"
#include "SundanceRivaraTreeNode.hpp"
#include "SundanceRivaraNode.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_RefCountPtr.hpp"

namespace Sundance
{
  namespace Rivara
  {
    class Element;
    class RivaraMesh;

    using Teuchos::Array;
    using Teuchos::RefCountPtr;

    /**
     * class Edge is a one-dimensional edge in a simplicial mesh.
     */

    class Edge : public TreeNode
    {
    public:
      /**
       * Construct with two nodes
       */
      Edge(const RCP<Node>& a,
        const RCP<Node>& b);

      /**
       * Add an element to the list of elements containing this edge
       */
      void addConnectingElement(Element* elem);

      /**
       * Return the length of the edge.
       */
      double length() const ;

      /**
       * Return a list of the cofacets of this edge that still need refinement
       */
      void getUnrefinedCofacets(Array<Element*>& elements) const ;

      /**
       * Bisect the edge.
       * @return a new node created at the midpoint of the edge
       */
      RCP<Node> bisect(RivaraMesh* mesh);

      const RCP<Node>& node(int i) const {return nodes_[i];}

      int ownerProc() const {return ownerProc_;}

      /**
       * Return the global index of this edge
       */
      int globalIndex() const ;

      /**
       * Set the global index of this edge
       */
      void setGlobalIndex(int globalIndex) ;

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
      Array<RCP<Node> > nodes_;
      Array<Element*> elements_;

      RCP<Node> midpoint_;

      int ownerProc_;
    };
  }
}

#endif
