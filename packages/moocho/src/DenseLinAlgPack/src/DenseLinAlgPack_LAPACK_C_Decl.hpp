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

#ifndef LAPACK_C_DECL_H
#define LAPACK_C_DECL_H

#include "Teuchos_F77_wrappers.h"

// C Declarations for calling LAPACK functions.

namespace LAPACK_C_Decl {

typedef FortranTypes::f_int			f_int;
typedef FortranTypes::f_real		f_real;
typedef FortranTypes::f_dbl_prec	f_dbl_prec;
typedef FortranTypes::f_logical		f_logical;
typedef FortranTypes::f_char		f_char;

// DPOTRF

void dpotrf(  const f_char& UPLO
  , const f_int& N, f_dbl_prec* A, const f_int& LDA
  , f_int* INFO  );

// DGEQRF

void dgeqrf( const f_int& M
  , const f_int& N, f_dbl_prec* A, const f_int& LDA
  , f_dbl_prec* TAU, f_dbl_prec* WORK
  , const f_int& LWORK, f_int* INFO  );

// DORMRQ

void dormqr( const f_char& SIDE
  , const f_char& TRANS, const f_int& M, const f_int& N
  , const f_int& K, const f_dbl_prec* A, const f_int& LDA
  , const f_dbl_prec* TAU, f_dbl_prec* C, const f_int& LDC
  , f_dbl_prec* WORK, const f_int& LWORK, f_int* INFO );

// DSYTRF

void dsytrf( const f_char& UPLO
  , const f_int& N, f_dbl_prec A[], const f_int& LDA
  , f_int IPIV[], f_dbl_prec WORK[], const f_int& LWORK
  , f_int* INFO );

// DSYTRS

void dsytrs( const f_char& UPLO
  , const f_int& N, const f_int& NRHS, const f_dbl_prec A[]
  , const f_int& LDA, const f_int IPIV[], f_dbl_prec B[]
  , const f_int& LDB, f_int* INFO );

// DGETRF

void dgetrf(
  const f_int& M, const f_int& N, f_dbl_prec A[], const f_int& LDA
  ,f_int IPIV[], f_int* INFO
  );

// DGETRS

void dgetrs(
  const f_char& TRANS
  ,const f_int& N, const f_int& NRHS, const f_dbl_prec A[]
  ,const f_int& LDA, const f_int IPIV[], f_dbl_prec B[]
  ,const f_int& LDB, f_int* INFO
  );


}	// end namespace LAPACK_C_Decl

#endif	// LAPACK_C_DECL_H
