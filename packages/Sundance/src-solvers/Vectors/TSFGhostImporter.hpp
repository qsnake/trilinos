/* ***********************************************************************
// 
//           TSFExtended: Trilinos Solver Framework Extended
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// **********************************************************************/

#ifndef TSFGHOSTIMPORTER_HPP
#define TSFGHOSTIMPORTER_HPP

#include "TSFVectorDecl.hpp"
#include "TSFVectorSpaceDecl.hpp"
#include "TSFGhostView.hpp"

namespace TSFExtended
{
  using namespace Teuchos;

  /**
   * In many applications it is necessary to view some subset of 
   * off-processor, or "ghost", elements of a vector. In 
   * matrix-vector multiplications,
   * access to off-processor elements is assumed to be handled internally
   * by the apply() method of LinearOp subtypes, so the TSF Vector type
   * does not need explicit accessors for ghost elements. However, in 
   * application interfaces such as finite element assembly engines, 
   * read-only access to ghost elements is sometimes required. The abstract
   * classes GhostImporter and GhostView define flexible interfaces
   * through which a set of required ghosts can be defined, ghost values
   * can be imported, and element values can be accessed through
   * global indices.
   *
   * Class GhostImporter is used to specify the set of ghost elements
   * that must be imported to this processor, and then to carry out the import.
   * It will often be the case that we do many imports with the same
   * set of ghost indices; for example, in a nonlinear problem the
   * import of the same set of ghost indices
   * will be repeated at each function evaluation. Therefore, it makes sense
   * to do the definition of the ghost index set and the import as
   * distinct methods. The definition of the ghost index set should be done
   * in the constructors of GhostImporter subclasses.
   */
  template <class Scalar>
  class GhostImporter
  {
  public:
    /** virtual dtor */
    virtual ~GhostImporter(){;}

    /** 
     * Import the ghost elements of the given vector
     * as specified during construction of this object. 
     */
    virtual void importView(const Vector<Scalar>& x,
                            RCP<GhostView<Scalar> >& ghostView) const = 0 ;

  };

}

#endif
