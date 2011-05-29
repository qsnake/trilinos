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

#ifndef MATRIX_SYM_SECANT_UPDATEABLE_H
#define MATRIX_SYM_SECANT_UPDATEABLE_H

#include <stdexcept>

#include "AbstractLinAlgPack_MatrixSymInitDiag.hpp"

namespace AbstractLinAlgPack {

/** \brief Mix-in interface for all polymorphic symmetric matrices that support secant updating.
 *
 * This interface is ment to be incorrporated in with a concrete <tt>AbstractLinAlgPack::MatrixOp</tt>
 * object that can implement some secant updating method.  Note that this is purely abstract interface\
 * and can be used in any application.
 *
 * Note that the methods <tt>AbstractLinAlgPack::MatrixSymInitDiag::init_identity()</tt> and
 * <tt>AbstractLinAlgPack::MatrixSymInitDiag::init_diagonal()</tt> do not state any postconditions
 * on the state of \c this after they are performed.  Subclasses of this interface also do not have
 * adhere to the obvious strick postconditions that these methods suggest but they should do the
 * "right thing" for the application.  If the client needs the obvious strict postconditions
 * for correct behavior, then the client is wise to test to see if \c this really is the
 * identity matrix or is a diagonal matrix (these can be cheap tests).
 */
class MatrixSymSecant
  : virtual public AbstractLinAlgPack::MatrixSymInitDiag // doxygen needs the full name
{
public:

  /** \brief . */
  class UpdateSkippedException : public std::runtime_error
  {public: UpdateSkippedException(const std::string& what_arg) : std::runtime_error(what_arg) {}};

  /** \brief . */
  class UpdateFailedException : public std::runtime_error
  {public: UpdateFailedException(const std::string& what_arg) : std::runtime_error(what_arg) {}};

  /** \brief Perform a secant update of the matrix.
   *
   * The update is a secant update:
   *
   * <tt>B_hat * s = y</tt>
   *
   * It is assumed that <tt>||B_hat - B||</tt> (here \c B is the matrix before the
   * update and \c B_hat is the matrix after the update) is not too large.
   *
   * The update vectors \c s and \c y may be used for workspace and are therefore
   * not gaurented to be preserved.
   *
   * The vector \c Bs may be set by the client to <tt>B*s</tt>.  This may help the
   * implementing subclass from having to compute it also.  Again, \c Bs may
   * be used as workspace by subclass so it may change.
   *
   * If the update is not performed then an <tt>UpdateSkippedException</tt>
   * will be thrown with an informative error message enclosed.
   *
   * If the update failed in catestrophic way and the state of <tt>this</tt>
   * is invalid, then a <tt>UpdateFailedException</tt> will be thrown with an
   * informative error message enclosed.
   *
   * Subclasses may also throw other unspecified exceptions but they should all
   * be derived from <tt>std::exception</tt>.
   * 
   * @param  s   [in/work] On input must be set to the \c s vector (see above).
   *             May be used as workspace.
   * @param  y   [in/work] On input must be set to the \c y vector (see above)
   *             May be used as workspace.
   * @param  Bs  [in/work] On input (if not \c NULL), \c Bs may be set to the
   *             <tt>B*s</tt> vector (see above).  If <tt>Bs == NULL</tt> on
   *             input, then the subclass implementation will do without.
   *             May be used as workspace.
   *
   * Preconditions:<ul>
   * <li> <tt>s != NULL</tt> (throw <tt>???</tt>)
   * <li> \c s must be compatible with \c space_rows() of the underlying matrix.(throw <tt>???</tt>)
   * <li> <tt>y != NULL</tt> (throw <tt>???</tt>)
   * <li> \c y must be compatible with \c space_cols() of the underlying matrix.(throw <tt>???</tt>)
   * <li> [<tt>Bs != NULL</tt>] \c Bs must be compatible with \c space_cols() of the underlying matrix.(throw <tt>???</tt>)
   * </ul>
   *
   * Postconidtons:<ul>
   * <li> <tt>(*this) * x \approx y</tt>.  In other words, we expect the secant condition will hold.
   *      In almost every application, this is required for correct behavior so clients should
   *      meet this condition.
   * </ul>
   */
  virtual void secant_update(
    VectorMutable     *s
    ,VectorMutable    *y
    ,VectorMutable    *Bs = NULL
    ) = 0;
  
};	// end class MatrixSymSecant 

}	// end namespace AbstractLinAlgPack 

#endif	// MATRIX_SYM_SECANT_UPDATEABLE_H
