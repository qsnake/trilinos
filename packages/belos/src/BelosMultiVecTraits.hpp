//@HEADER
// ************************************************************************
//
//                 Belos: Block Linear Solvers Package
//                  Copyright 2004 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER
//
#ifndef BELOS_MULTI_VEC_TRAITS_HPP
#define BELOS_MULTI_VEC_TRAITS_HPP

/*! \file BelosMultiVecTraits.hpp
    \brief Virtual base class which defines basic traits for the multivector type
*/

#include "Teuchos_RCP.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "BelosTypes.hpp"

namespace Belos {

  /*! \brief This is the default struct used by MultiVecTraits<ScalarType, MV> class to produce a
      compile time error when the specialization does not exist for multivector type <tt>MV</tt>.
  */
  template< class ScalarType, class MV >
  struct UndefinedMultiVecTraits
  {
    //! This function should not compile if there is an attempt to instantiate!
    /*! \note Any attempt to compile this function results in a compile time error.  This means
      that the template specialization of Belos::MultiVecTraits class for type <tt>MV</tt> does
      not exist, or is not complete.
    */
    static inline ScalarType notDefined() { return MV::this_type_is_missing_a_specialization(); };
  };

  /*! \brief Virtual base class which defines basic traits for the multi-vector type.

      An adapter for this traits class must exist for the <tt>MV</tt> type.
      If not, this class will produce a compile-time error.

      \ingroup belos_opvec_interfaces
  */
  template<class ScalarType, class MV>
  class MultiVecTraits 
  {
  public:

    //! @name Creation methods
    //@{

    /*! \brief Creates a new empty \c MV containing \c numvecs columns.

    \return Reference-counted pointer to the new multivector of type \c MV.
    */
    static Teuchos::RCP<MV> Clone( const MV& mv, const int numvecs )
    { UndefinedMultiVecTraits<ScalarType, MV>::notDefined(); return Teuchos::null; }     

    /*! \brief Creates a new \c MV and copies contents of \c mv into the new vector (deep copy).

      \return Reference-counted pointer to the new multivector of type \c MV.
    */
    static Teuchos::RCP<MV> CloneCopy( const MV& mv )
    { UndefinedMultiVecTraits<ScalarType, MV>::notDefined(); return Teuchos::null; }     

    /*! \brief Creates a new \c MV and copies the selected contents of \c mv into the new vector (deep copy).

      The copied vectors from \c mv are indicated by the \c index.size() indices in \c index.
      \return Reference-counted pointer to the new multivector of type \c MV.
    */
    static Teuchos::RCP<MV> CloneCopy( const MV& mv, const std::vector<int>& index )
    { UndefinedMultiVecTraits<ScalarType, MV>::notDefined(); return Teuchos::null; }     

    /*! \brief Creates a new \c MV that shares the selected contents of \c mv (shallow copy).

    The index of the \c numvecs vectors shallow copied from \c mv are indicated by the indices given in \c index.
    \return Reference-counted pointer to the new multivector of type \c MV.
    */
    static Teuchos::RCP<MV> CloneViewNonConst( MV& mv, const std::vector<int>& index )
    { UndefinedMultiVecTraits<ScalarType, MV>::notDefined(); return Teuchos::null; }     

    /*! \brief Creates a new const \c MV that shares the selected contents of \c mv (shallow copy).

    The index of the \c numvecs vectors shallow copied from \c mv are indicated by the indices given in \c index.
    \return Reference-counted pointer to the new const multivector of type \c MV.
    */
    static Teuchos::RCP<const MV> CloneView( const MV& mv, const std::vector<int>& index )
    { UndefinedMultiVecTraits<ScalarType, MV>::notDefined(); return Teuchos::null; }     

    //@}

    //! @name Attribute methods
    //@{

    //! Obtain the vector length of \c mv.
    static int GetVecLength( const MV& mv )
    { UndefinedMultiVecTraits<ScalarType, MV>::notDefined(); return 0; }     

    //! Obtain the number of vectors in \c mv
    static int GetNumberVecs( const MV& mv )
    { UndefinedMultiVecTraits<ScalarType, MV>::notDefined(); return 0; }     

    //@}

    //! @name Update methods
    //@{

    /*! \brief Update \c mv with \f$ \alpha AB + \beta mv \f$.
     */
    static void MvTimesMatAddMv( const ScalarType alpha, const MV& A, 
                                 const Teuchos::SerialDenseMatrix<int,ScalarType>& B, 
                                 const ScalarType beta, MV& mv )
    { UndefinedMultiVecTraits<ScalarType, MV>::notDefined(); }     

    /*! \brief Replace \c mv with \f$\alpha A + \beta B\f$.
     */
    static void MvAddMv( const ScalarType alpha, const MV& A, const ScalarType beta, const MV& B, MV& mv )
    { UndefinedMultiVecTraits<ScalarType, MV>::notDefined(); }     

    /*! \brief Scale each element of the vectors in \c mv with \c alpha.
     */
    static void MvScale ( MV& mv, const ScalarType alpha )
    { UndefinedMultiVecTraits<ScalarType, MV>::notDefined(); }
    
    /*! \brief Scale each element of the \c i-th vector in \c mv with \c alpha[i].
     */
    static void MvScale ( MV& mv, const std::vector<ScalarType>& alpha )
    { UndefinedMultiVecTraits<ScalarType, MV>::notDefined(); }

    /*! \brief Compute a dense matrix \c B through the matrix-matrix multiply \f$ \alpha A^Hmv \f$.
    */
    static void MvTransMv( const ScalarType alpha, const MV& A, const MV& mv, Teuchos::SerialDenseMatrix<int,ScalarType>& B)
    { UndefinedMultiVecTraits<ScalarType, MV>::notDefined(); }     

    /*! \brief Compute a vector \c b where the components are the individual dot-products of the \c i-th columns of \c A and \c mv, i.e.\f$b[i] = A[i]^Hmv[i]\f$.
     */
    static void MvDot ( const MV& mv, const MV& A, std::vector<ScalarType> &b) 
    { UndefinedMultiVecTraits<ScalarType, MV>::notDefined(); }     

    //@}
    //! @name Norm method
    //@{

    /*! \brief Compute the 2-norm of each individual vector of \c mv.
      Upon return, \c normvec[i] holds the value of \f$||mv_i||_2\f$, the \c i-th column of \c mv.
    */
    static void MvNorm( const MV& mv, std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType>& normvec, NormType type = TwoNorm )
    { UndefinedMultiVecTraits<ScalarType, MV>::notDefined(); }     

    //@}

    //! @name Initialization methods
    //@{
    /*! \brief Copy the vectors in \c A to a set of vectors in \c mv indicated by the indices given in \c index.

    The \c numvecs vectors in \c A are copied to a subset of vectors in \c mv indicated by the indices given in \c index,
    i.e.<tt> mv[index[i]] = A[i]</tt>.
    */
    static void SetBlock( const MV& A, const std::vector<int>& index, MV& mv )
    { UndefinedMultiVecTraits<ScalarType, MV>::notDefined(); }     

    /*! \brief Replace the vectors in \c mv with random vectors.
     */
    static void MvRandom( MV& mv )
    { UndefinedMultiVecTraits<ScalarType, MV>::notDefined(); }     

    /*! \brief Replace each element of the vectors in \c mv with \c alpha.
     */
    static void MvInit( MV& mv, const ScalarType alpha = Teuchos::ScalarTraits<ScalarType>::zero() )
    { UndefinedMultiVecTraits<ScalarType, MV>::notDefined(); }     

    //@}

    //! @name Print method
    //@{

    /*! \brief Print the \c mv multi-vector to the \c os output stream.
     */
    static void MvPrint( const MV& mv, std::ostream& os )
    { UndefinedMultiVecTraits<ScalarType, MV>::notDefined(); }     

    //@}
  };
  
} // namespace Belos

#endif // BELOS_MULTI_VEC_TRAITS_HPP
