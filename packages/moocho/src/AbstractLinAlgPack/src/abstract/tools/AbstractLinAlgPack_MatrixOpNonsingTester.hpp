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

#ifndef MATRIX_WITH_OP_NONSINGULAR_TESTER_H
#define MATRIX_WITH_OP_NONSINGULAR_TESTER_H

#include <iosfwd>

#include "AbstractLinAlgPack_Types.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

namespace AbstractLinAlgPack {

/** \brief Testing class for \c MatrixOpNonsing interface.
 *
 * This testing class is basically a unit tester for \c MatrixOpNonsing.  The method \c test_matrix()
 * runs several different tests to check that \f$ M M^{-1} \approx I \f$ and \f$ M^{-T} M^T \approx I \f$ using
 * randomly generated vectors \a v and the methods \c MatrixNonsing::V_InvMtV() and \c MatrixOp::Vp_StMtV().
 * These test should only be performed, of course, on a fully initialized <tt>%MatrixOpNonsing</tt> object.
 *
 * The tests performed by this testing class are designed to allow some validation for even the larges systems
 * and will produce various levels of output so as to be usefull in debugging.
 *
 * ToDo:  Finish documentation!
 */
class MatrixOpNonsingTester {
public:

  /** @name Public types */
  //@{

  /** \brief . */
  enum ETestLevel {
    TEST_LEVEL_2_BLAS  = 1  ///< Test Level-2 BLAS operations only
    ,TEST_LEVEL_3_BLAS = 2  ///< Test Level-2 and Level-3 BLAS operations
  };
  /** \brief . */
  enum EPrintTestLevel {
    PRINT_NONE   = 0  ///< Don't print anything
    ,PRINT_BASIC = 1  ///< Print only very basic info
    ,PRINT_MORE  = 2  ///< Print greater detail about the tests.
    ,PRINT_ALL   = 3  ///< Print everything all the tests in great detail but output is independent of problem size.
  };

  //@}

  /** @name Set and access options */
  //@{

  /// Set the level of testing
  STANDARD_MEMBER_COMPOSITION_MEMBERS( ETestLevel, test_level );
  /// Set the level of output produced durring tests.
  STANDARD_MEMBER_COMPOSITION_MEMBERS( EPrintTestLevel, print_tests );
  /// Set whether vectors etc. are printed (warning, this may be a lot of output for larger systems).
  STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, dump_all );
  /// Set whether an exception that is thrown is thrown clear out of the testing function or not.
  STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, throw_exception );
  /// Set the number of random test cases created.
  STANDARD_MEMBER_COMPOSITION_MEMBERS( size_type, num_random_tests );
  /// Set the relative tolerance for numerical tests above which to print a warning.
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, warning_tol );
  /// Set the relative tolerance for numerical tests above which to return false from the testing function.
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, error_tol );

  //@}

  /** @name Constructors / initializers */
  //@{

  ///	Constructor (default options)
  MatrixOpNonsingTester(
    ETestLevel       test_level       = TEST_LEVEL_2_BLAS
    ,EPrintTestLevel print_tests      = PRINT_NONE
    ,bool            dump_all         = false
    ,bool            throw_exception  = true
    ,size_type       num_random_tests = 1
    ,value_type      warning_tol      = 1e-14
    ,value_type      error_tol        = 1e-8
    );

  //@}

  /** @name Test MatrixOpNonsing */
  //@{
 
  /** \brief Test a \c MatrixOpNonsing object.
   *
   * @param  M    [in] Matrix object being tested.
   * @param  M_name
   *              [in] Name given to the matrix object to be used in outputting and error
   *              reporting (i.e. throwing an exception).
   * @param  out  [in/out] If <tt>out != NULL</tt> any and all output will be sent here.  If
   *              <tt>out == NULL</tt> then no output will be produced.
   *
   * @return Returns \c true if all of the tests checked out and no unexpected exceptions were
   * thrown.
   *
   * The behavior of this method depends on a set of options and the input arguments.
   * <ul>
   * <li> <b><tt>throw_exception(bool)</tt></b>:
   *      If <tt>throw_exception()</tt> == true</tt>, then if any of the objects within
   *      this function throw exceptions, these exceptions will be be thrown clean
   *      out of this function for the caller to handle.  If <tt>throw_exception()</tt> == false</tt>,
   *      then if any object throws an exception, the exception is caught and this this function will
   *      return <tt>false</tt>.  In any case an error message will be printed
   *      to <tt>*out</tt> (if <tt>out != NULL</tt) before leaving the function (by \c return or \c throw).
   * <li> <b><tt>dump_all(bool)</tt></b>:
   *      If <tt>dump_all() == true</tt> then all of the computed quantities will but dumped to \c out.
   *      Note that this is a useful option for initial debugging of small systems but not a good idea for
   *      larger systems as it will result in an excessive amount of output.
   * <li> ToDo: Add rest of options!
   * </ul>
   */
  bool test_matrix(
    const MatrixOpNonsing   &M
    ,const char                     M_name[]
    ,std::ostream                   *out
    );

  //@}

}; // end class MatrixOpNonsingTester

} // end namespace AbstractLinAlgPack

#endif // MATRIX_WITH_OP_NONSINGULAR_TESTER_H
