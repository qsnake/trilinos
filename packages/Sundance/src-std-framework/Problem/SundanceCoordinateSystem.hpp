/* ========================================================================
 *
 *                           Catfish 
 * 
 *           A Finite Element Schrodinger-Poisson Solver
 *
 * Copyright (C) 2008, Kevin Long, Texas Tech University. 
 *
 * This library is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2.1 of the
 * License, or (at your option) any later version.
 *  
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *  
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 * USA
 * 
 * ======================================================================== */


#ifndef SUNDANCE_COORDINATESYSTEM_HPP
#define SUNDANCE_COORDINATESYSTEM_HPP

#include "SundanceExpr.hpp"
#include "SundanceCoordExpr.hpp"
#include "SundanceExceptions.hpp"
#include <cmath>


namespace Sundance
{
using Sundance::Expr;
using namespace Teuchos;

class CoordinateSystemBase : 
  public Sundance::Handleable<CoordinateSystemBase>
{
public:
  /** */
  virtual ~CoordinateSystemBase(){}

  /** */
  virtual Expr jacobian() const = 0 ;

  /** */
  virtual bool supportsMeshDimension(int dim) const = 0 ;

  /** */
  virtual std::ostream& print(std::ostream& os) const = 0 ;

  /** */
  double pi() const {return 4.0*::atan(1.0);}
};


class CoordinateSystem : public Sundance::Handle<CoordinateSystemBase>
{
public:
  /* boilerplate handle ctors */
  HANDLE_CTORS(CoordinateSystem, CoordinateSystemBase);

  /** */
  Expr jacobian() const {return ptr()->jacobian();}

  /** */
  bool supportsMeshDimension(int dim) const 
    {return ptr()->supportsMeshDimension(dim);}

};


class CartesianCoordinateSystem : public CoordinateSystemBase
{
public:
  /** */
  CartesianCoordinateSystem(){}

  /** */
  Expr jacobian() const {return Expr(1.0);}

  /** */
  bool supportsMeshDimension(int dim) const {return dim > 0;}

  /** */
  std::ostream& print(std::ostream& os) const 
    {
      os << "Cartesian";
      return os;
    };

  
  /** */
  virtual RCP<CoordinateSystemBase> getRcp() {return rcp(this);}
};


class MeridionalCylindricalCoordinateSystem : public CoordinateSystemBase
{
public:
  /** */
  MeridionalCylindricalCoordinateSystem()
    : r_(new CoordExpr(0)) {}

  /** */
  Expr jacobian() const {return pi()*r_;}

  /** */
  bool supportsMeshDimension(int dim) const {return dim > 0 && dim <= 2;}

  /** */
  std::ostream& print(std::ostream& os) const 
    {
      os << "Meridional Cylindrical";
      return os;
    };
  
  /** */
  virtual RCP<CoordinateSystemBase> getRcp() {return rcp(this);}

private:
  Expr r_;
};

class RadialSphericalCoordinateSystem : public CoordinateSystemBase
{
public:
  /** */
  RadialSphericalCoordinateSystem()
    : r_(new CoordExpr(0)) {}

  /** */
  Expr jacobian() const {return pi()*r_*r_;}

  /** */
  bool supportsMeshDimension(int dim) const {return dim == 1;}
  
  /** */
  std::ostream& print(std::ostream& os) const 
    {
      os << "Radial Spherical";
      return os;
    };



  
  /** */
  virtual RCP<CoordinateSystemBase> getRcp() {return rcp(this);}
private:
  Expr r_;
};


class CoordinateSystemBuilder
{
public:
  /** */
  static CoordinateSystem makeCoordinateSystem(const std::string& name);
};


inline std::ostream& operator<<(std::ostream& os, const CoordinateSystem& cs)
{
  return cs.ptr()->print(os);
}

}




#endif
