/*@HEADER
// ***********************************************************************
// 
//                Komplex: Complex Linear Solver Package
//                 Copyright (2002) Sandia Corporation
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
// ***********************************************************************
//@HEADER
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "az_aztec.h"
#include "azk_komplex.h"
/*! \file
\brief Destruction routine for deleting Komplex systems.

KOMPLEX is an add-on module to AZTEC that allows users to solve
complex-valued
linear systems.

KOMPLEX solves a complex-valued linear system Ax = b by solving
an equivalent real-valued system of twice the dimension.
Specifically,
writing in terms of real and imaginary parts, we have

 \f[ (A_r + i*A_i)*(x_r + i*x_i) = (b_r + i*b_i) \f]

  or by separating into real and imaginary equations we have

\f[
  \left( \begin{array}{rr}
                                    A_r & -A_i\\
                                    A_i &  A_r
                             \end{array}
   \right)
   \left( \begin{array}{r}
                                    x_r\\
                                    x_i
                             \end{array}
   \right)
   =
   \left( \begin{array}{r}
                                    b_r\\
                                    b_i
                             \end{array}
   \right)
\f]
  which is a real-valued system of twice the size.  If we find xr and
xi, we
  can form the solution to the original system as x = xr +i*xi.


KOMPLEX accept user linear systems in three forms with either global
or local index values.

1) The first form is true complex.  The user passes in an MSR or VBR
format matrix where the values are stored like Fortran complex
numbers.
Thus, the values array is of type double that is twice as long as the
number of complex values.  Each complex entry is stored with real part
followed by imaginary part (as in Fortran).

2) The second form stores real and imaginary parts separately, but the
pattern for each is identical.  Thus only the values of the imaginary
part are passed to the creation routines.

3) The third form accepts two real-valued matrices with no assumption
about the structure of the matrices.  Each matrix is multiplied by a
user-supplied complex constant.  This is the most general form.

Each of the above forms supports a global or local index set.  By this
we mean that the index values (stored in bindx) refer to the global
problem indices, or the local indices (for example after calling
AZ_transform).

*/


/*! \fn void AZK_destroy_linsys( int *options, double *params,
		       int *proc_config,
		       double **x, double **b,
		       AZ_MATRIX **Amat_komplex)
\brief Destroy a Komplex System.

Destroys a komplex system created by any of the AZK_create_linsys
functions.  Deletes any memory allocated by creation routine.

\param options (In)
       Determines specific solution method and other parameters.
\param params (In)
       Drop tolerance and convergence tolerance info.
\param proc_config (In)
       Machine configuration.  proc_config[AZ_node] is the node
       number.  proc_config[AZ_N_procs] is the number of processors.

\param x (Out)
       Deleted komplex version of solution.  Remember to call
       AZK_extract_solution_[k2c,g2k,ri2k] before calling this routine.
\param b (Out)
       Deleted komplex version of RHS.
\param Amat_komplex (Out)
       Deleted komplex version of matrix stored as an AZ_MATRIX structure.

*/

void AZK_destroy_linsys( int *options, double *params,
		       int *proc_config,
		       double **x, double **b,
		       AZ_MATRIX **Amat_komplex)
{
  AZ_KOMPLEX *linsys_pass_data;
  int *komplex_to_real, *komplex_to_imag;

  linsys_pass_data = (AZ_KOMPLEX *) (*Amat_komplex)->aux_ptr;

  if (linsys_pass_data->Form_of_Equations != AZK_Komplex_No_Copy)
    {
      /* Destroy RHS, initial guess and matrix */
      
      AZK_destroy_vector( options, params, proc_config, (*Amat_komplex), x);
      AZK_destroy_vector( options, params, proc_config, (*Amat_komplex), b);
      AZK_destroy_matrix( options, params, proc_config, Amat_komplex);
    }
  else
    {
      komplex_to_real = linsys_pass_data->komplex_to_real;
      komplex_to_imag = linsys_pass_data->komplex_to_imag;
    

      /* Free allocated memory */
      
      AZ_free((void *) komplex_to_real);
      AZ_free((void *) komplex_to_imag);
      AZ_free ((void **) x); 
      AZ_free ((void **) b);
      AZ_free((void *) linsys_pass_data);

      /* Free data_org if Aztec doesn't do it */
      if (!(*Amat_komplex)->must_free_data_org) 
	AZ_free((void *) (*Amat_komplex)->data_org);
      
      AZ_matrix_destroy (Amat_komplex);
    }

}
