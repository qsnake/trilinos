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

#ifndef SUNDANCE_CELLREORDERERIMPLEMBASE_H
#define SUNDANCE_CELLREORDERERIMPLEMBASE_H


#include "SundanceDefs.hpp"
#include "SundanceNoncopyable.hpp"
#include "SundanceObjectWithVerbosity.hpp"
#include <typeinfo>

namespace Sundance
{
class MeshBase;

/**
 * Abstract interface for the low-level objects that 
 * implement cell reordering. 
 * 
 * <h4> Adding a new reordering algorithm </h4>
 *
 * To add a new reordering algorithm, you should create a new
 * subclass of CellReordererImplemBase. The only method you will
 * need to implement is
 * \code
 * virtual int advance(int currentLID) const 
 * \endcode
 * which should provide the maximal cell LID found after
 * the <tt>currentLID.</tt>
 * Depending on the algorithm , you may also want to override
 * the methods
 * \code
 * virtual int begin() const 
 * virtual int end() const 
 * \endcode
 * which return the index of the first cell to be processed,
 * and a past-the-end index. 
 */
class CellReordererImplemBase 
  : public ObjectWithClassVerbosity<CellReordererImplemBase>
{
public:
  /** Construct with a pointer to a mesh */
  CellReordererImplemBase(const MeshBase* mesh);
      
  /** virtual dtor */
  virtual ~CellReordererImplemBase(){;}

  /** return a descriptive std::string */
  virtual std::string typeName() const {return typeid(*this).name();}
    
  /** */
  virtual int advance(int currentLID) const = 0 ;
      
  /** */
  virtual int begin() const {return 0;}
      
  /** */
  virtual int end() const ;
protected:
  /** */
  const MeshBase* mesh() const {return mesh_;}

private:
  /** Unmanaged pointer to a mesh. The mesh will contain a smart
   * pointer to this reorderer, so to avoid closed reference
   * graphs we store a raw pointer here.*/
  const MeshBase* mesh_;
      
};

}



#endif
