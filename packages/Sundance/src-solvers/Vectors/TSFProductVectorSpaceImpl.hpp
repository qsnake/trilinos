/* @HEADER@ */
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
 /* @HEADER@ */


#ifndef TSFPRODUCTVECTORSPACEIMPL_HPP
#define TSFPRODUCTVECTORSPACEIMPL_HPP

#include "Thyra_DefaultProductVectorSpace.hpp"
#include "TSFVectorSpaceDecl.hpp"
 
using namespace TSFExtended;
using namespace Teuchos;



namespace TSFExtended
{

/** */
template <class Scalar> inline
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > 
productSpace(const Array<VectorSpace<Scalar> >& spaces)
{
  Array<RCP<const Thyra::VectorSpaceBase<Scalar> > > data(spaces.size());
  for (int i=0; i<spaces.size(); i++)
  {
    data[i] = spaces[i].ptr();
  }
  return rcp(new Thyra::DefaultProductVectorSpace<Scalar>(data.size(), &(data[0])));
}


/** */
template <class Scalar> inline
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > 
productSpace(VectorSpace<Scalar>& s1)
{
  Array<VectorSpace<Scalar> > s;
  s.append(s1);
  return productSpace(s);
}

/** */
template <class Scalar> inline
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > 
productSpace(VectorSpace<Scalar>& s1, 
  VectorSpace<Scalar>& s2)
{
  Array<VectorSpace<Scalar> > s;
  s.append(s1);
  s.append(s2);
  return productSpace(s);
}

/** */
template <class Scalar> inline
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > 
productSpace(VectorSpace<Scalar>& s1,VectorSpace<Scalar>& s2,
  VectorSpace<Scalar>& s3)
{
  Array<VectorSpace<Scalar> > s;
  s.append(s1);
  s.append(s2);
  s.append(s3);
  return productSpace(s);
}

  
  
} 






#endif
