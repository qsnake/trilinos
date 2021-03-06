/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2010 Sandia National Laboratories.  Developed at the
    University of Wisconsin--Madison under SNL contract number
    624796.  The U.S. Government and the University of Wisconsin
    retain certain rights to this software.

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License 
    (lgpl.txt) along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

    (2009) kraftche@cae.wisc.edu    
    (2010) jwfrank@sandia.gov   

  ***************************************************************** */


/** \file Target3DShapeSizeUntangle.hpp
 *  \brief 
 *  \author Jason Franks
 */

#ifndef MSQ_TARGET_3D_SHAPE_SIZE_UNTANGLE_HPP
#define MSQ_TARGET_3D_SHAPE_SIZE_UNTANGLE_HPP

#include "Mesquite.hpp"
#include "TargetMetric3D.hpp"

namespace MESQUITE_NS {

/**\brief Untangle metric
 *
 * \f$ \mu(T) = |T|^3 - 3 \sqrt{3} \tau + \gamma (\tau - 1)^2 \f$
 *
 * \f$ \f(mu(T)) = \{|[1 - \epsilon] - \mu(T)| - ([1 - \epsilon] - \mu(T))\}^2
 *
 */
class Target3DShapeSizeUntangle : public TargetMetric3D
{
private:
  double mEps;
  double mGamma;

public:

  Target3DShapeSizeUntangle( double eps = 0.01, double gamma = 2.0 ) : mEps(eps), mGamma(gamma) {}

  MESQUITE_EXPORT virtual
  ~Target3DShapeSizeUntangle();

  MESQUITE_EXPORT virtual
  std::string get_name() const;

  MESQUITE_EXPORT virtual
  bool evaluate( const MsqMatrix<3,3>& A, 
                 const MsqMatrix<3,3>& W, 
                 double& result, 
                 MsqError& err );
   
  MESQUITE_EXPORT virtual
  bool evaluate_with_grad( const MsqMatrix<3,3>& A, 
                           const MsqMatrix<3,3>& W, 
                           double& result, 
                           MsqMatrix<3,3>& wrt_A, 
                           MsqError& err );

  /* MESQUITE_EXPORT virtual
  bool evaluate_with_hess( const MsqMatrix<3,3>& A,
                           const MsqMatrix<3,3>& W,
                           double& result,
                           MsqMatrix<3,3>& deriv_wrt_A,
                           MsqMatrix<3,3> second_wrt_A[3],
                           MsqError& err );*/
};



} // namespace MESQUITE_NS

#endif
