/*
 * SundanceRefinementClass.hpp
 *
 *  Created on: Apr 29, 2010
 *      Author: benk
 */

#ifndef SUNDANCEREFINEMENTCLASS_HPP_
#define SUNDANCEREFINEMENTCLASS_HPP_

#include "SundanceRefinementBase.hpp"
#include "SundanceHandle.hpp"

namespace Sundance{

using namespace Teuchos;

class RefinementClass : public Sundance::Handle<RefinementBase> {
public:

	/* Handle constructors */
	HANDLE_CTORS(RefinementClass, RefinementBase);

	/** see RefinementBase for Docu */
	int refine( const int cellLevel ,  const Point& cellPos ,
			     const Point& cellDimameter) const {
		return ptr()->refine(cellLevel , cellPos , cellDimameter );
	}

	/** see RefinementBase for Docu */
	int estimateRefinementLevel( const Point& cellPos, const Point& cellDimameter) const {
		return ptr()->estimateRefinementLevel( cellPos, cellDimameter);
	}

private:

};
}

#endif /* SUNDANCEREFINEMENTCLASS_HPP_ */
