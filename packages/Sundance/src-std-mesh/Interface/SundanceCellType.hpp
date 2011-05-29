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

#ifndef CELLTOPOLOGYCODE_H
#define CELLTOPOLOGYCODE_H


#include "SundanceDefs.hpp"

namespace Sundance 
{

/** \defgroup Sundance_CellType_grp Cell Type Description
 */

/** \brief Enumeration of specific 1D, 2D, and 3D cell types.
 *
 * See the implementation of the nonmember functions <tt>dimension()</tt>,
 * <tt>numFacets()</tt>, and <tt>facetType()</tt> a full description of what
 * these cell types are.
 *
 * ToDo: Consider making CellType an abstract base class so that any type of
 * cell can be implemented?
 *
 * \ingroup Sundance_CellType_grp
 */
enum CellType {
  NullCell        ///< No cell specified
  ,PointCell      ///< 0D vertex cell
  ,LineCell       ///< 1D line, or edge, cell
  ,TriangleCell   ///< 2D triangle
  ,TetCell        ///< 3D tetrahedral cell
  ,QuadCell       ///< 2D quadrilateral cell
  ,BrickCell      ///< 3D "brick" cell
  ,PrismCell      ///< 3D prism cell
};

/** \brief Return a std::string representation of the cell type.
 *
 * \ingroup Sundance_CellType_grp
 */
std::string toString(const CellType& c) ;

/** \brief Return the dimension of the cell type.
 *
 * \ingroup Sundance_CellType_grp
 */
int dimension(const CellType& c) ;

/** \brief Return the number of faces of a given facet dimension for a cell type.
 *
 * <b>Preconditions:</b><ul>
 * <li><tt>facetDim <= dimension(c)</tt>
 * </ul>
 *
 * \returns -1 for null cells and point cells
 *
 * \ingroup Sundance_CellType_grp
 */
int numFacets(const CellType& c, int facetDim);

/** \return Return the type of facet of a given facet dimension and a
 * particualar facet.
 *
 * \param  c
 *           [in] Cell type
 * \param  facetDim
 *           [in] The dimension of the facet type requested
 * \param  facetIndex
 *           [in] The particualar index of the facet as defined
 *           by the convension of the cell type (which is actually
 *           partially defined in this function).
 *
 * <b>Preconditions:</b><ul>
 * <li><tt>facetDim <= dimension(c)</tt>
 * <li><tt>0 <= facetIndex < numFacets(c,facetDim)</tt>
 * </ul>
 *
 * \ingroup Sundance_CellType_grp
 */
CellType facetType(const CellType& c, int facetDim, int facetIndex);

/** \brief output stream operator for <tt>CellType</tt>.
 *
 * \ingroup Sundance_CellType_grp
 */
inline std::ostream& operator<<(std::ostream& os, const CellType& c)
{
  os << toString(c);
  return os;
}
  
} // namespace Sundance


#endif


