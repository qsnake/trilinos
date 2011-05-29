// ///////////////////////////////////////////////////////////
// AbstractLinAlgPack_TestMatrixSymSecant.hpp

#include <iosfwd>

#include "AbstractLinAlgPack_Types.hpp"

namespace AbstractLinAlgPack {

/** \brief Checks the secant condition <tt>B*s = y</tt>.
 *
 * Call this function after calling <tt>B.secant_update(s,y,...)</tt>.
 *
 * Returns \c true if the secant condition satsifies \c error_tol
 * and \c false otherwise.
 *
 * @param  B  [in] The matrix object we are testing.
 * @param  s  [in] First secant vector
 * @param  y  [in] Second secant vector
 * @param  warning_tol
 *            [in] Any relative error above \c warning_tol will
 *            be noted and possibly printed.
 * @param  error_tol
 *            [in] Any relative error above \c error_tol will
 *            cause the tests to stop and false to
 *            be returned.
 * @param  print_all_warnings
 *            [in] If true then all relative errors greater than
 *            warning_tol will be printed.
 * @param  out
 *            [in/out] Stream that output or error messages are sent
 *            to follow the tests.
 * @param  trase
 *            [in] If \c trase==true then tests will be trased and
 *            if \c trase==false then only error messages will be
 *            output.
 */
bool TestMatrixSymSecant(
  const MatrixOp        &B
  ,const Vector       &s
  ,const Vector       &y
  ,value_type               warning_tol
  ,value_type               error_tol
  ,bool                     print_all_warnings
  ,std::ostream             *out
  ,bool                     trase                 = true
  );

} // end namespace AbstractLinAlgPack
