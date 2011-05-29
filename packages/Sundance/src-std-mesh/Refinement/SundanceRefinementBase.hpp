/*
 * SundanceRefinementBase.hpp
 *
 *  Created on: Apr 29, 2010
 *      Author: benk
 */

#ifndef SUNDANCEREFINEMENTBASE_HPP_
#define SUNDANCEREFINEMENTBASE_HPP_

#include "SundanceDefs.hpp"
#include "SundanceHandleable.hpp"
#include "SundancePoint.hpp"

namespace Sundance {

using namespace Teuchos;

/** define the predicate */

#define REFINE_MESH_(name, code) \
  class name : public RefinementBase, \
               public Sundance::Handleable<RefinementBase> \
  { \
  public:\
    name() : RefinementBase(){;}            \
    virtual ~name(){;} \
    virtual int refine( const int cellLevel , \
						 const Point& cellPos , \
						 const Point& cellDimameter) const code \
    GET_RCP(RefinementBase);\
  }

#define REFINE_MESH(name, code) REFINE_MESH_(name, code);


/** define the predicate also with estimation */
#define REFINE_MESH_ESTIMATE_(name, code , codeEst ) \
  class name : public RefinementBase, \
               public Sundance::Handleable<RefinementBase> \
  { \
  public:\
    name() : RefinementBase(){;}            \
    virtual ~name(){;} \
    virtual int refine( const int cellLevel , \
						 const Point& cellPos , \
						 const Point& cellDimameter) const code \
	virtual int estimateRefinementLevel( const Point& cellPos , \
						 const Point& cellDimameter) const codeEst \
    GET_RCP(RefinementBase);\
  }

#define REFINE_MESH_ESTIMATE(name, code, codeEst) REFINE_MESH_ESTIMATE_(name, code, codeEst);


/** Base class for mesh refinement , but also to define computational domain
 * (which must not be the same as the mesh domain) <br>
 * It is important to take a look at the refinement protocol defined:<br>
 *  0 -> no action , 1 -> refine , 2 -> coarse  <br>
 *  <br> <br>
 *  */
class RefinementBase  {

public:

	RefinementBase() {;}

	virtual ~RefinementBase() {;}

	/** Function which has to answer the question per Cell if the
	 * cell needs to be refined or not
	 * @param cellLevel [in] refinement level of the actual cell
	 * @param cellPos [in] the position of the cell (the middle point)
	 * @param cellDimameter [in] resolution of the cell in each direction, this could also be CellDiameter
	 * @returns refinementAction [out] , 0 -> no action , 1 -> refine , 2 -> coarse */
	virtual int refine(const int cellLevel ,
			            const Point& cellPos ,
			            const Point& cellDimameter) const { return 1; }

	/** Function meant to be used for load estimation of a given cell <br>
	 * This function has a "dummy" body so the sub classes must not overwrite this function
	 * @param cellPos [in] the position of the cell (the middle point)
	 * @param cellDimameter [in] resolution of the cell in each direction (this parameter might not be used)
	 * @returns averageLevel [out] returns the estimated (average) refinement level of this cell */
	virtual int estimateRefinementLevel(
			             const Point& cellPos,
			             const Point& cellDimameter ) const { return 1; }

};
}

#endif /* SUNDANCEREFINEMENTBASE_HPP_ */
