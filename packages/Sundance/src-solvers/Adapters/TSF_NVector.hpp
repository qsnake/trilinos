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

#ifndef TSF_NVECTOR_HPP
#define TSF_NVECTOR_HPP

#include "SundanceDefs.hpp"
#include "TSFVectorDecl.hpp"
#include "TSFVectorSpaceDecl.hpp"

#ifdef HAVE_SUNDIALS
#include "nvector.h"
#include "sundialstypes.h"

 /*
  * TSF implementation of SUNDIALS NVector interface
  */

class _TrilinosVectorContent
{
public:
  _TrilinosVectorContent() : data() {;}
  TSFExtended::Vector<realtype> data;
};

typedef struct _TrilinosVectorContent* TrilinosVectorContent;


extern "C"
{

#define NV_CONTENT_Trilinos(v)  ( (TrilinosVectorContent)(v->content) );

  /*
   *
   */
  inline Vector<realtype>& toTrilinos(N_Vector x)
  {
    TrilinosVectorContent xPtr = reinterpret_cast<TrilinosVectorContent>(x->content);

    return xPtr->data;
  }

  /*
   * -----------------------------------------------------------------
   * Function : N_VNew_Serial
   * -----------------------------------------------------------------
   * This function creates and allocates memory for a serial vector.
   * -----------------------------------------------------------------
   */

  N_Vector N_VNew_Trilinos(const VectorSpace<realtype>& space);

  /*
   * -----------------------------------------------------------------
   * Function : N_VPrint_Trilinos
   * -----------------------------------------------------------------
   * This function prints the content of a Trilinos vector to stdout.
   * -----------------------------------------------------------------
   */

  void N_VPrint_Trilinos(N_Vector v);

  /*
   * -----------------------------------------------------------------
   * Trilinos implementations of various useful vector operations
   * -----------------------------------------------------------------
   */

  N_Vector N_VCloneEmpty_Trilinos(N_Vector w);
  N_Vector N_VClone_Trilinos(N_Vector w);
  void N_VDestroy_Trilinos(N_Vector v);
  void N_VSpace_Trilinos(N_Vector v, long int *lrw, long int *liw);
  realtype *N_VGetArrayPointer_Trilinos(N_Vector v);
  void N_VSetArrayPointer_Trilinos(realtype *v_data, N_Vector v);
  void N_VLinearSum_Trilinos(realtype a, N_Vector x, realtype b, 
                             N_Vector y, N_Vector z);
  void N_VConst_Trilinos(realtype c, N_Vector z);
  void N_VProd_Trilinos(N_Vector x, N_Vector y, N_Vector z);
  void N_VDiv_Trilinos(N_Vector x, N_Vector y, N_Vector z);
  void N_VScale_Trilinos(realtype c, N_Vector x, N_Vector z);
  void N_VAbs_Trilinos(N_Vector x, N_Vector z);
  void N_VInv_Trilinos(N_Vector x, N_Vector z);
  void N_VAddConst_Trilinos(N_Vector x, realtype b, N_Vector z);
  realtype N_VDotProd_Trilinos(N_Vector x, N_Vector y);
  realtype N_VMaxNorm_Trilinos(N_Vector x);
  realtype N_VWrmsNorm_Trilinos(N_Vector x, N_Vector w);
  realtype N_VWrmsNormMask_Trilinos(N_Vector x, N_Vector w, N_Vector id);
  realtype N_VMin_Trilinos(N_Vector x);
  realtype N_VWL2Norm_Trilinos(N_Vector x, N_Vector w);
  realtype N_VL1Norm_Trilinos(N_Vector x);
  void N_VCompare_Trilinos(realtype c, N_Vector x, N_Vector z);
  booleantype N_VInvTest_Trilinos(N_Vector x, N_Vector z);
  booleantype N_VConstrMask_Trilinos(N_Vector c, N_Vector x, N_Vector m);
  realtype N_VMinQuotient_Trilinos(N_Vector num, N_Vector denom);
}
  
#endif // HAVE_SUNDIALS

#endif
