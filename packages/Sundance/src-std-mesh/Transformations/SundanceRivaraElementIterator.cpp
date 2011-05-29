#include "SundanceRivaraElementIterator.hpp"
#include "SundanceRivaraMesh.hpp"
#include "SundanceRivaraElement.hpp"
#include "SundanceRivaraNode.hpp"


using namespace Sundance::Rivara;
using namespace Teuchos;


ElementIterator::ElementIterator(const RivaraMesh* mesh)
	: mesh_(const_cast<RivaraMesh*>(mesh)), rootIndex_(-1), current_(0),
		startingNewTree_(true)
{}

bool ElementIterator::hasMoreElements() const 
{
	if (rootIndex_ < mesh_->numRootElements()-1)
		{
			return true;
		}
	else if (current_->next() != 0)
		{
			return true;
		}
	else return false;
}

const Element* ElementIterator::getNextElement() const 
{
	if (current_==0)
		{
			current_ = const_cast<Element*>(dynamic_cast<const Element*>(mesh_->element(0)->first()));
			rootIndex_ = 0;
		}
	else if (current_->next() != 0)
		{
			current_ = const_cast<Element*>(dynamic_cast<const Element*>(current_->next()));
		}
	else
		{
			rootIndex_++;
			current_ = const_cast<Element*>(dynamic_cast<const Element*>(mesh_->element(rootIndex_)->first()));
		}
	
	return current_;
	
}

	
