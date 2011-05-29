#include "TSF_NVector.hpp"
#include "TSFLinearCombinationDecl.hpp"
#include "Thyra_SUNDIALS_Ops.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "TSFVectorImpl.hpp"
#include "TSFLinearCombinationImpl.hpp"
#endif

using namespace TSFExtended;
using namespace TSFExtendedOps;

#ifdef TRILINOS_6

#ifdef HAVE_SUNDIALS
#include "sundialsmath.h"
#include "sundialstypes.h"

extern "C"
{

  /*
   * -----------------------------------------------------------------
   * exported functions
   * -----------------------------------------------------------------
   */

  /* ----------------------------------------------------------------------------
   * Function to create a new Trilinos vector
   */

  N_Vector N_VNew_Trilinos(const VectorSpace<realtype>& space)
  {
    /* Create vector */
    N_Vector v = new _generic_N_Vector(); 
    if (v == NULL) return(NULL);

    /* Create vector operation structure */
    N_Vector_Ops ops = new _generic_N_Vector_Ops();
    if (ops == NULL) 
      {
        delete v; 
        return(NULL);
      }


    ops->nvclone           = N_VClone_Trilinos;
    ops->nvcloneempty      = 0;//N_VCloneEmpty_Trilinos;
    ops->nvdestroy         = N_VDestroy_Trilinos;
    ops->nvspace           = N_VSpace_Trilinos;
    ops->nvgetarraypointer = 0;//N_VGetArrayPointer_Trilinos;
    ops->nvsetarraypointer = 0;//N_VSetArrayPointer_Trilinos;
    ops->nvlinearsum       = N_VLinearSum_Trilinos;
    ops->nvconst           = N_VConst_Trilinos;
    ops->nvprod            = N_VProd_Trilinos;
    ops->nvdiv             = N_VDiv_Trilinos;
    ops->nvscale           = N_VScale_Trilinos;
    ops->nvabs             = N_VAbs_Trilinos;
    ops->nvinv             = N_VInv_Trilinos;
    ops->nvaddconst        = N_VAddConst_Trilinos;
    ops->nvdotprod         = N_VDotProd_Trilinos;
    ops->nvmaxnorm         = N_VMaxNorm_Trilinos;
    ops->nvwrmsnormmask    = N_VWrmsNormMask_Trilinos;
    ops->nvwrmsnorm        = N_VWrmsNorm_Trilinos;
    ops->nvmin             = N_VMin_Trilinos;
    ops->nvwl2norm         = N_VWL2Norm_Trilinos;
    ops->nvl1norm          = N_VL1Norm_Trilinos;
    ops->nvcompare         = N_VCompare_Trilinos;
    ops->nvinvtest         = N_VInvTest_Trilinos;
    ops->nvconstrmask      = N_VConstrMask_Trilinos;
    ops->nvminquotient     = N_VMinQuotient_Trilinos;

    /* Create content */
    TrilinosVectorContent content = new _TrilinosVectorContent();
    if (content == NULL) 
      {
        delete ops; 
        delete v; 
        return(NULL);
      }

    content->data = space.createMember();

    /* Attach content and ops */
    v->content = content;
    v->ops = ops;

    return(v);
  }



  void N_VPrint_Trilinos(N_Vector x)
  {
    cout << toTrilinos(x) << std::endl;
  }

  /*
   * -----------------------------------------------------------------
   * implementation of vector operations
   * -----------------------------------------------------------------
   */


  N_Vector N_VClone_Trilinos(N_Vector w)
  {
    N_Vector v;

    Vector<realtype>& original = toTrilinos(w);
    v = N_VNew_Trilinos(original.space());
    toTrilinos(v).acceptCopyOf(original);

    return(v);
  }

  void N_VDestroy_Trilinos(N_Vector v)
  {
    TrilinosVectorContent xPtr = reinterpret_cast<TrilinosVectorContent>(v->content);
    delete xPtr;
    delete v->ops;
    delete v;
  }

  void N_VSpace_Trilinos(N_Vector v, long int *lrw, long int *liw)
  {
    *lrw = toTrilinos(v).space().dim();
    *liw = 1;
  }


  void N_VLinearSum_Trilinos(realtype a, N_Vector x, 
                             realtype b, N_Vector y, N_Vector z)
  {
    if (z==x)
      {
        toTrilinos(z).update(b, toTrilinos(y), a);
      }
    else if (z==y)
      {
        toTrilinos(z).update(a, toTrilinos(x), b);
      }
    else
      {
        toTrilinos(z).update(a, toTrilinos(x), 
                             b, toTrilinos(y), 0.0);
      }
  }

  void N_VConst_Trilinos(realtype c, N_Vector z)
  {
    toTrilinos(z).setToConstant(c);
  }

  void N_VProd_Trilinos(N_Vector x, N_Vector y, N_Vector z)
  {
    Thyra::VProd(*(toTrilinos(x).ptr()),
                 *(toTrilinos(y).ptr()),
                 toTrilinos(z).ptr().get());
  }

  void N_VDiv_Trilinos(N_Vector x, N_Vector y, N_Vector z)
  {
    Thyra::VDiv(*(toTrilinos(x).ptr()),
                *(toTrilinos(y).ptr()),
                toTrilinos(z).ptr().get());
  }

  void N_VScale_Trilinos(realtype c, N_Vector x, N_Vector z)
  {
    Thyra::VScale(c, *(toTrilinos(x).ptr()), toTrilinos(z).ptr().get());
  }

  void N_VAbs_Trilinos(N_Vector x, N_Vector z)
  {
    /* notice that the argument order is reversed. This is not a typo. */
    Thyra::abs(toTrilinos(z).ptr().get(), *(toTrilinos(x).ptr()));
  }

  void N_VInv_Trilinos(N_Vector x, N_Vector z)
  {
    /* notice that the argument order is reversed. This is not a typo. */
    Thyra::reciprocal(toTrilinos(z).ptr().get(), *(toTrilinos(x).ptr()));
  }

  void N_VAddConst_Trilinos(N_Vector x, realtype b, N_Vector z)
  {
    Thyra::VAddConst(b, *(toTrilinos(x).ptr()), toTrilinos(z).ptr().get());
  }

  realtype N_VDotProd_Trilinos(N_Vector x, N_Vector y)
  {
    realtype rtn = toTrilinos(x) * toTrilinos(y);
    return rtn;
  }

  realtype N_VMaxNorm_Trilinos(N_Vector x)
  {
    realtype rtn = toTrilinos(x).normInf();
    return rtn;
  }

  realtype N_VWrmsNorm_Trilinos(N_Vector x, N_Vector w)
  {
    realtype rtn = Thyra::VWrmsNorm(*(toTrilinos(x).ptr()),
                                    *(toTrilinos(w).ptr()));
    return rtn;
  }

  realtype N_VWrmsNormMask_Trilinos(N_Vector x, N_Vector w, N_Vector id)
  {
    realtype rtn = Thyra::VWrmsMaskNorm(*(toTrilinos(x).ptr()),
                                        *(toTrilinos(w).ptr()),
                                        *(toTrilinos(id).ptr()));
    return rtn;
  }

  realtype N_VMin_Trilinos(N_Vector x)
  {
    realtype rtn = toTrilinos(x).min();
    return rtn;
  }

  realtype N_VWL2Norm_Trilinos(N_Vector x, N_Vector w)
  {
    realtype rtn = Thyra::VWL2Norm(*(toTrilinos(x).ptr()),
                                   *(toTrilinos(w).ptr()));
    return rtn;
  }

  realtype N_VL1Norm_Trilinos(N_Vector x)
  {
    realtype rtn = toTrilinos(x).norm1();
    return rtn;
  }

  void N_VCompare_Trilinos(realtype c, N_Vector x, N_Vector z)
  {
    Thyra::VCompare(c, *(toTrilinos(x).ptr()),
                    toTrilinos(z).ptr().get());
  }

  booleantype N_VInvTest_Trilinos(N_Vector x, N_Vector z)
  {
    booleantype rtn = Thyra::VInvTest(*(toTrilinos(x).ptr()),
                                      toTrilinos(z).ptr().get());
    return rtn;
  }

  booleantype N_VConstrMask_Trilinos(N_Vector c, N_Vector x, N_Vector m)
  {
    booleantype rtn = Thyra::VConstrMask(*(toTrilinos(x).ptr()),
                                         *(toTrilinos(c).ptr()),
                                         toTrilinos(m).ptr().get());
    return rtn;
  }

  realtype N_VMinQuotient_Trilinos(N_Vector x, N_Vector y)
  {
    realtype rtn = Thyra::VMinQuotient(*(toTrilinos(x).ptr()),
                                       *(toTrilinos(y).ptr()));
    return rtn;
  }

}

#endif
#endif
