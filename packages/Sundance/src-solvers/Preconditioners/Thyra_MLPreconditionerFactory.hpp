/*@HEADER
// ***********************************************************************
// 
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
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

#ifndef THYRA_ML_PRECONDITIONER_FACTORY_DECL_HPP
#define THYRA_ML_PRECONDITIONER_FACTORY_DECL_HPP

#include "SundanceDefs.hpp"

#ifndef TRILINOS_6

#include "Thyra_PreconditionerFactoryBase.hpp"
#include "Thyra_EpetraOperatorViewExtractorBase.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"

namespace Thyra {

  using Teuchos::ParameterList;
  using Teuchos::RefCountPtr;

  /** These can be used to choose one of ML's predefined default parameter sets */
  enum EMLProblemType {ML_SmoothedAggregation, 
                       ML_DomainDecomposition,
                       ML_DomainDecompositionML,
                       ML_Maxwell};

/** \brief Concrete preconditioner factory subclass based on ML.
 *
 * ToDo: Finish documentation!
 */
class MLPreconditionerFactory : public PreconditionerFactoryBase<double> {
public:

  /** @name Constructors/initializers/accessors */
  //@{

  /** \brief . */
  MLPreconditionerFactory(const RCP<ParameterList>& params);
  /** \brief . */
  MLPreconditionerFactory(const EMLProblemType& probType,
                          const ParameterList& revisions=ParameterList());

  /** */
  MLPreconditionerFactory(const std::string& probType,
                          const ParameterList& revisions=ParameterList());
    
  /** \brief Set the strategy object used to extract an
   * <tt>Epetra_Operator</tt> view of an input forward operator.
   *
   * This view will then be dynamically casted to <tt>Epetra_RowMatrix</tt>
   * before it is used.
   *
   * The default implementation used is <tt>EpetraOperatorViewExtractorBase</tt>.
   */
  STANDARD_COMPOSITION_MEMBERS( EpetraOperatorViewExtractorBase, epetraFwdOpViewExtractor )

  //@}

  /** @name Overridden from PreconditionerFactoryBase */
  //@{

  /** \brief . */
  bool isCompatible( const LinearOpBase<double> &fwdOp ) const;
  /** \brief . */
  bool applySupportsConj(EConj conj) const;
  /** \brief . */
  bool applyTransposeSupportsConj(EConj conj) const;
  /** \brief . */
  RCP<PreconditionerBase<double> > createPrec() const;
  /** \brief . */
  void initializePrec(
    const RCP<const LinearOpBase<double> >    &fwdOp
    ,PreconditionerBase<double>                                *prec
    ,const ESupportSolveUse                                    supportSolveUse
    ) const;
  /** \brief . */
  void uninitializePrec(
    PreconditionerBase<double>                          *prec
    ,RCP<const LinearOpBase<double> >  *fwdOp
    ,ESupportSolveUse                                   *supportSolveUse
    ) const;

  //@}

  /** @name Overridden from ParameterListAcceptor */
  //@{

  /** \brief . */
  void setParameterList(RCP<ParameterList> const& paramList);
  /** \brief . */
  RCP<ParameterList> getNonconstParameterList();
  /** \brief . */
  RCP<ParameterList> unsetParameterList();
  /** \brief . */
  RCP<const ParameterList> getParameterList() const;
  //@}

  /** \name Public functions overridden from Describable. */
  //@{

  /** \brief . */
  std::string description() const;

  //@}

private:

  RCP<ParameterList> reviseDefaultList(const ParameterList& defaults, 
                                               const ParameterList& revisions) const;

  std::string probToString(const EMLProblemType& probType) const ;

  RCP<ParameterList> defaultParameters(const EMLProblemType& probType) const ;

  RCP<ParameterList> defaultParameters(const std::string& probType) const ;

  // ////////////////////////////////
  // Private data members

  RCP<ParameterList>       paramList_;


  // ////////////////////////////////
  // Private member functions

  static RCP<const ParameterList> generateAndGetValidParameters();

};

} // namespace Thyra

#endif
#endif // THYRA_ML_PRECONDITIONER_FACTORY_DECL_HPP
