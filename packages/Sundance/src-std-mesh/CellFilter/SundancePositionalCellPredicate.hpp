/* @HEADER@ */
// ************************************************************************
// 
//                              Sundance
//                 Copyright (2005) Sandia Corporation
// 
// Copyright (year first published) Sandia Corporation.  Under the terms 
// of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government 
// retains certain rights in this software.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//                                                                                 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA                                                                                
// Questions? Contact Kevin Long (krlong@sandia.gov), 
// Sandia National Laboratories, Livermore, California, USA
// 
// ************************************************************************
/* @HEADER@ */


#ifndef SUNDANCE_POSITIONALCELLPREDICATE_H
#define SUNDANCE_POSITIONALCELLPREDICATE_H


#include "SundanceDefs.hpp"
#include "SundanceCellPredicateBase.hpp"

namespace Sundance
{
using namespace Teuchos;

#define NEW_CELL_PREDICATE(name)  \
  class name : public CellPredicateFunctorBase,  \
               public Sundance::Handleable<CellPredicateFunctorBase>  \
  {  \
  public:  \
    name() : CellPredicateFunctorBase(#name) {} \
    virtual ~name() {}  \
    virtual bool operator()(const Point& x) const; \
    GET_RCP(CellPredicateFunctorBase);  \
  }; \
  \
  bool name::operator()(const Point& x) const


#define CELL_PREDICATE_(name, code) \
  class name : public CellPredicateFunctorBase, \
               public Sundance::Handleable<CellPredicateFunctorBase> \
  { \
  public:\
    name() : CellPredicateFunctorBase(#name){;}            \
    virtual ~name(){;}\
    virtual bool operator()(const Point& x) const code \
    GET_RCP(CellPredicateFunctorBase);\
  }


#define CELL_PREDICATE(name, code) CELL_PREDICATE_(name, code);




/** */
class CellPredicateFunctorBase
{
public:
  /** */
  CellPredicateFunctorBase(const std::string& name="Functor(" + Teuchos::toString(topID()) + ")")
    : name_(name) {;}

  /** */
  virtual ~CellPredicateFunctorBase(){;}

  /** */
  virtual bool operator()(const Point& x) const = 0 ;

  /** */
  virtual std::string description() const {return name_;}
private:
  static int& topID() {static int rtn=0; rtn++; return rtn;}
  std::string name_;
};

  
/** 
 * PositionalCellPredicate tests whether the cell's nodes satisfy
 * a condition on their positions.
 */
class PositionalCellPredicate : public CellPredicateBase 
{
public:
      
  /** Construct with a function of positions */
  PositionalCellPredicate(const RCP<CellPredicateFunctorBase>& func) 
    : CellPredicateBase(), func_(func) 
    {;}

  /** virtual dtor */
  virtual ~PositionalCellPredicate(){;}
      
  /** Test the predicate on a batch of cells */
  virtual void testBatch(const Array<int>& cellLID,
    Array<int>& results) const ;

  /** Write to XML */
  virtual XMLObject toXML() const ;

  /** comparison */
  virtual bool lessThan(const CellPredicateBase* other) const ;

  /** */
  virtual std::string description() const {return func_->description();}

  /* */
  GET_RCP(CellPredicateBase);

private:
  RCP<CellPredicateFunctorBase> func_;
};





/** */
class PointCellPredicateFunctor 
  : public CellPredicateFunctorBase
{
public:
  /** */
  PointCellPredicateFunctor(const Point& x, const double& tol=1.0e-12)
    : x_(x), tol_(tol){}

  /** */
  bool operator()(const Point& x) const ;

private:
  Point x_;
  double tol_;
};



/** */
class CoordinateValueCellPredicateFunctor 
  : public CellPredicateFunctorBase
{
public:
  /** */
  CoordinateValueCellPredicateFunctor(
    int direction, const double& value, const double& tol=1.0e-12)
    : direction_(direction), value_(value), tol_(tol) {}

  /** */
  bool operator()(const Point& x) const ;

private:
  int direction_;
  double value_;
  double tol_;
};

/** */
class PointCellPredicate : public PositionalCellPredicate
{
public:
  /** */
  PointCellPredicate(const Point& x, const double& tol=1.0e-12)
    : PositionalCellPredicate(rcp(new PointCellPredicateFunctor(x,tol)))
    {}

  /* */
  GET_RCP(CellPredicateBase);
};

/** */
class CoordinateValueCellPredicate : public PositionalCellPredicate
{
public:
  /** */
  CoordinateValueCellPredicate(int direction,
    const double& value, const double& tol=1.0e-12)
    : PositionalCellPredicate(
      rcp(new CoordinateValueCellPredicateFunctor(direction,value,tol)))
    {}

  /* */
  GET_RCP(CellPredicateBase);
};


}


#endif
