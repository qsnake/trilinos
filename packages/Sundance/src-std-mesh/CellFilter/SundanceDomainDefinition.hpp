/*
 * SundanceDomainDefinition.hpp
 *
 *  Created on: Apr 30, 2010
 *      Author: benk
 */

#ifndef SUNDANCEDOMAINDEFINITION_HPP_
#define SUNDANCEDOMAINDEFINITION_HPP_

#include "SundanceDefs.hpp"
#include "SundanceHandleable.hpp"
#include "SundanceHandle.hpp"
#include "SundancePoint.hpp"

#include "SundanceParametrizedCurve.hpp"

namespace Sundance {

using namespace Teuchos;

/** define the predicate also with estimation */
#define MESH_DOMAIN_(name, code) \
  class name : public MeshDomainBase, \
               public Sundance::Handleable<MeshDomainBase> \
  { \
  public:\
    name() : MeshDomainBase(){;}            \
	virtual bool isInsideComputationalDomain(const Point& x) const code \
    GET_RCP(MeshDomainBase);\
  }

#define MESH_DOMAIN(name, code) MESH_DOMAIN_(name, code);


/** Base class for mesh refinement , but also to define computational domain
 * (which must not be the same as the mesh domain)*/
class MeshDomainBase  {

public:

	MeshDomainBase() {;}

	virtual ~MeshDomainBase() {;}

	/** Function to answer if the domain is inside the computational domain <br>
	 * The strategy should be that if one point of a cell is inside the domain, then the whole
	 * cell should be considered as in the computational domain.
	 * @param x [in] coordinate of the point
	 * @return true if the point is inside the domain, false otherwise */
	virtual bool isInsideComputationalDomain(const Point& x) const { return true;}
};

// ---------------

/**  Class defines mesh domain based on parametrized curve */
class CurveDomain : public MeshDomainBase ,
                    public Sundance::Handleable<MeshDomainBase>{
public:

	/** Ctor with the 2 necessary input arguments */
	CurveDomain(const ParametrizedCurve& curve ,
			    CurveCellFilterMode mode) :
		MeshDomainBase() , curve_(curve) , mode_(mode){;}
    /** empty Dtor */
	virtual ~CurveDomain() {;}

	/**  in or outside the domain */
	virtual bool isInsideComputationalDomain(const Point& x) const {
		if (mode_ == Outside_Curve){
            return (curve_.curveEquation(x) >= -1e-8 );
		}
		else{
			return (curve_.curveEquation(x) <= 1e-8 );
		}
	}

	GET_RCP(MeshDomainBase);

private:

	const ParametrizedCurve& curve_;

	const CurveCellFilterMode mode_;
};

// ---------------

class MeshDomainDef : public Sundance::Handle<MeshDomainBase> {
public:

	/* Handle constructors */
	HANDLE_CTORS(MeshDomainDef, MeshDomainBase);

	/** see MeshDomainBase for Docu */
	bool isInsideComputationalDomain(const Point& x) const {
		return ptr()->isInsideComputationalDomain(x);
	}

private:

};


}

#endif /* SUNDANCEDOMAINDEFINITION_HPP_ */
