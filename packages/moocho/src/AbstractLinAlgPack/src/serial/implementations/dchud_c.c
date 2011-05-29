/* translated by f2c from dchud.f and hand modified. */

#include "Moocho_Config.h"


#if !defined(HAVE_MOOCHO_FORTRAN)


#include "Teuchos_BLAS_wrappers.hpp"
#include <math.h>


void dchud_c(double r__[], int *ldr, int *p,
  double x[], double z__[], int *ldz, int *nz,
  double y[], double *rho, double c__[], double s[]
  )
/*
double *r__;
int *ldr, *p;
double *x, *z__;
int *ldz, *nz;
double *y, *rho, *c__, *s;
*/
{
  /* System generated locals */
  int r_dim1, r_offset, z_dim1, z_offset, i__1, i__2;
  double d__1, d__2;

  /* Local variables */
  double zeta;
  int i__, j;
  double t, scale, azeta;
  double xj;
  int jm1;

/* ***FIRST EXECUTABLE STATEMENT  DCHUD */
  /* Parameter adjustments */
  r_dim1 = *ldr;
  r_offset = 1 + r_dim1 * 1;
  r__ -= r_offset;
  --x;
  z_dim1 = *ldz;
  z_offset = 1 + z_dim1 * 1;
  z__ -= z_offset;
  --y;
  --rho;
  --c__;
  --s;

  /* Function Body */
  i__1 = *p;
  for (j = 1; j <= i__1; ++j) {
    xj = x[j];

/*    APPLY THE PREVIOUS ROTATIONS. */

    jm1 = j - 1;
    if (jm1 < 1) {
      goto L20;
    }
    i__2 = jm1;
    for (i__ = 1; i__ <= i__2; ++i__) {
      t = c__[i__] * r__[i__ + j * r_dim1] + s[i__] * xj;
      xj = c__[i__] * xj - s[i__] * r__[i__ + j * r_dim1];
      r__[i__ + j * r_dim1] = t;
/* L10: */
    }
L20:

/*    COMPUTE THE NEXT ROTATION. */

    DROTG_F77(&r__[j + j * r_dim1], &xj, &c__[j], &s[j]);
/* L30: */
  }

/*   IF REQUIRED, UPDATE Z AND RHO. */

  if (*nz < 1) {
    goto L70;
  }
  i__1 = *nz;
  for (j = 1; j <= i__1; ++j) {
    zeta = y[j];
    i__2 = *p;
    for (i__ = 1; i__ <= i__2; ++i__) {
      t = c__[i__] * z__[i__ + j * z_dim1] + s[i__] * zeta;
      zeta = c__[i__] * zeta - s[i__] * z__[i__ + j * z_dim1];
      z__[i__ + j * z_dim1] = t;
/* L40: */
    }
    azeta = fabs(zeta);
    if (azeta == 0. || rho[j] < 0.) {
      goto L50;
    }
    scale = azeta + rho[j];
/* Computing 2nd power */
    d__1 = azeta / scale;
/* Computing 2nd power */
    d__2 = rho[j] / scale;
    rho[j] = scale * sqrt(d__1 * d__1 + d__2 * d__2);
L50:
/* L60: */
    ;
  }
L70:

  return;

} /* dchud_ */


#endif // !defined(HAVE_MOOCHO_FORTRAN)
