#ifndef SUNDANCERIVARATREENODE_H
#define SUNDANCERIVARATREENODE_H

#include "SundanceDefs.hpp"

namespace Sundance
{
  namespace Rivara
    {
      /**
       * Class TreeNode represents a node in a Rivara mesh refinement tree.
       * Each node has either zero or two children depending on whether
       * it has been refined. All non-root nodes have a pointer back to their
       * parent.
       *
       * Only maximal elements will be responsible for deleting their children.
       * Therefore, the TreeNode dtor does not delete children; subtypes
       * that need to delete children should call the deleteChildren()
       * method.
       */
      class TreeNode 
        {
        public:
          /** Empty ctor */
          TreeNode();

          virtual ~TreeNode(){;}

          /** Delete the node's children */
          void deleteChildren();


          /** set the parent of this node */
          void setParent(TreeNode* parent) {parent_ = parent;}

          /** Set the two children of this node */
          void setChildren(TreeNode* left, TreeNode* right)
            {left_ = left; right_ = right;}

          /** return the leftmost leaf beneath this node */
          const TreeNode* first() const ;

          /** return the rightmost leaf beneath this node */
          const TreeNode* last() const ;

          /** Indicate whether this is the leftward child of another node */
          bool isLeftChild() const ;

          /** Indicate whether this is the rightward child of another node */
          bool isRightChild() const ;

          /** Return the next leaf in a left-to-right walk of the tree. If this
           * is the last leaf, return 0. */
          const TreeNode* next() const ;

          /** Indicate whether this node has children */
          bool hasChildren() const {return left_ != 0;}

          /** Return a count of the number of leaves */
          int numLeaves() const ;


        protected:
          TreeNode* left() {return left_;}

          TreeNode* right() {return right_;}
        private:

          TreeNode* parent_;

          TreeNode* left_;

          TreeNode* right_;
        };
    }
}

#endif
