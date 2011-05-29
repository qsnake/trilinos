#include "SundanceRivaraTreeNode.hpp"
#include "SundanceExceptions.hpp"

using namespace Sundance::Rivara;
using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;

TreeNode::TreeNode()
  : parent_(0), left_(0), right_(0)
{;}

void TreeNode::deleteChildren()
{
  if (left_ != 0) delete left_;
  if (right_ != 0) delete right_;
}

const TreeNode* TreeNode::first() const
{
  /* keep walking left as long as possible */
  if (left_ != 0) return left_->first();
  return this;
}

const TreeNode* TreeNode::last() const
{
  /* keep walking right as long as possible */
  if (right_ != 0) return right_->last();
  return this;
}

bool TreeNode::isRightChild() const
{
  /* return true if I am the right child of my parent */
  if (parent_!=0 && parent_->right_ == this) return true;
  return false;
}

bool TreeNode::isLeftChild() const
{
  /* return true if I am the left child of my parent */
  if (parent_!=0 && parent_->left_ == this) return true;
  return false;
}

const TreeNode* TreeNode::next() const 
{

  /* walk up the tree until we are at either the root or are the
   * left child */
  TreeNode* pos = const_cast<TreeNode*>(this);

  while (pos->isRightChild())
    {
      pos = pos->parent_;
    }

  /* if we are the left child, begin at the leftmost leaf of
   * the right sibling tree */
  if (pos->isLeftChild())
    {
      return pos->parent_->right_->first();
    }

  /* if we are at the root, there are no more unwalked leaves */
  if (pos->parent_==0) return pos->parent_;

  /* if we get to this point, there's a bug in this code */
  TEST_FOR_EXCEPT(true);

  return pos->parent_;
}

int TreeNode::numLeaves() const
{
  if (hasChildren()) 
		{
			return left_->numLeaves() + right_->numLeaves();
		}
  return 1;
}





