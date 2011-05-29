// @HEADER
// ***********************************************************************
// 
// Moocho: Multi-functional Object-Oriented arCHitecture for Optimization
//                  Copyright (2003) Sandia Corporation
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
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef DELETE_ROW_COL_H
#define DELETE_ROW_COL_H

#include "DenseLinAlgPack_Types.hpp"

namespace DenseLinAlgPack {

/** \brief . */
/* * Delete a symmetric row and a column form a triangular matrix.
 *
 * If #M# is a lower triangular matrix then we partition it
 * as:
 \verbatim

   1 |\
     |  \
     |    \
     | M11  \
     |________\ _
  kd |_________|_|
     |         | |\
     |         | |  \
     |         | |    \
     |   M31   | | M33  \
   n |         | |        \
     ----------------------
     1         kd         n

 \endverbatim
 * In order to delete row #kd# and column #kd# the rectangular
 * matrix #M31# is moved up one row and the triangular matrix
 * #M33# is moved up one row and to the left one column.
 *
 * If #M# is an upper triangular matrix then we partition it
 * as:
 \verbatim

  1         kd      n
  -------------------- 1
  \        | |       |
    \  M11 | |  M13  |
      \    | |       |
        \  | |       |
          \|_|_______|
           |_|_______| kd
             \       |
               \ M33 |
                 \   |
                   \ | n
 
 \endverbatim
 *
 * In order to delete row #kd# and column #kd# the matrix
 * #M13# is moved one column to the left and the upper
 * triangular matrix #M33# is moved one row up and
 * on column to the left.
 *
 * Preconditions:<ul>
 * <li> #M != NULL#
 * <li> #M->rows() >= 1#
 * <li> #1 <= kd <= M->rows()#
 * </ul>
 */
void delete_row_col( size_type kd, DMatrixSliceTriEle* M );

} // end namespace DenseLinAlgPack

#endif  // DELETE_ROW_COL_H
