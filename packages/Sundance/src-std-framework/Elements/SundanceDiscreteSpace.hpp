/* @HEADER@ */
// ************************************************************************
// 
//                              Sundance
//                 Copyright (2005) Sandia Corporation
// 
// Copyright (year first published) Sandia Corporation.  Under the terms 
// of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government 
// retains certain rights in this software.
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
// Questions? Contact Kevin Long (krlong@sandia.gov), 
// Sandia National Laboratories, Livermore, California, USA
// 
// ************************************************************************
/* @HEADER@ */

#ifndef SUNDANCE_DISCRETESPACE_H
#define SUNDANCE_DISCRETESPACE_H

#include "SundanceDefs.hpp"
#include "SundanceMesh.hpp"
#include "SundanceBasisFamily.hpp"
#include "SundanceSpectralBasis.hpp"
#include "SundanceDOFMapBase.hpp"
#include "TSFVectorType.hpp"
#include "TSFVectorSpaceDecl.hpp"
#include "SundanceObjectWithVerbosity.hpp"
#include "SundanceDiscreteSpaceTransfBuilder.hpp"


namespace Sundance
{
class FunctionSupportResolver;
}

namespace Sundance
{
  using namespace Teuchos;
  using namespace TSFExtended;

  /** 
   * DiscreteSpace represents a discrete finite-element space (i.e., 
   * a mesh and a basis).
   */
  class DiscreteSpace
  {
  public:
    /** */
    DiscreteSpace(){;}
    /** */
    DiscreteSpace(const Mesh& mesh, const BasisFamily& basis,
      const VectorType<double>& vecType,
      int setupVerb = 0);
    /** */
    DiscreteSpace(const Mesh& mesh, const BasisArray& basis,
                  const VectorType<double>& vecType,
      int setupVerb = 0);

    /** */
    DiscreteSpace(const Mesh& mesh, const BasisArray& basis,
                  const Array<CellFilter>& regions,
                  const VectorType<double>& vecType,
      int setupVerb = 0);


    /** */
    DiscreteSpace(const Mesh& mesh, const BasisFamily& basis,
                  const CellFilter& regions,
                  const VectorType<double>& vecType,
      int setupVerb = 0);


    /** */
    DiscreteSpace(const Mesh& mesh, const BasisArray& basis,
                  const CellFilter& regions,
                  const VectorType<double>& vecType,
      int setupVerb = 0);

    /** */
    DiscreteSpace(const Mesh& mesh, const BasisArray& basis,
                  const RCP<DOFMapBase>& map,
                  const VectorType<double>& vecType,
      int setupVerb = 0);

    /** */
    DiscreteSpace(const Mesh& mesh, const BasisArray& basis,
      const RCP<DOFMapBase>& map,
      const RCP<Array<int> >& bcIndices,
      const VectorType<double>& vecType,
      int setupVerb = 0);

    /** */
    DiscreteSpace(const Mesh& mesh, const BasisFamily& basis,
                  const SpectralBasis& spBasis,
                  const VectorType<double>& vecType,
      int setupVerb = 0);
    /** */
    DiscreteSpace(const Mesh& mesh, const BasisArray& basis,
                  const SpectralBasis& spBasis,
                  const VectorType<double>& vecType,
      int setupVerb = 0);

    /** */
    DiscreteSpace(const Mesh& mesh, const BasisArray& basis,
      const RCP<FunctionSupportResolver>& fsr,
      const VectorType<double>& vecType,
      int setupVerb = 0);

    /** */
    const RCP<DOFMapBase>& map() const {return map_;}

    /** return the number of functions */
    int nFunc() const {return basis_.size();}

    /** */
    const BasisArray& basis() const {return basis_;}

    /** */
    Array<std::pair<int,int> > dimStructure() const {return vectorDimStructure(basis());}

    /** */
    Vector<double> createVector() const {return vecSpace_.createMember();}

    /** */
    VectorSpace<double> vecSpace() const {return vecSpace_;}

    /** */
    const Mesh& mesh() const {return mesh_;}

    /** */
    const VectorType<double>& vecType() const {return vecType_;}

    /** */
    void importGhosts(const Vector<double>& x,
                      RCP<GhostView<double> >& ghostView) const ;

    /** */
    void getAllowedFuncs(const CellFilter& cf, Set<int>& funcs) const ;

    /** */
    const CellFilter& cellFilters(int i) const {return subdomains_[i];}

    /** Return the transformation builder */
    const RCP<DiscreteSpaceTransfBuilder>& getTransformation() const
    		{ return transformationBuilder_; }

  private:

    /** */
    void init(const Array<CellFilter>& regions,
              const BasisArray& basis);

    /** */
    void init(const Array<CellFilter>& regions,
      const BasisArray& basis,
      const RCP<Array<int> >& isBCIndex, 
      bool partitionBCs);

    /** */
    Array<CellFilter> maximalRegions(int n) const ;

    /** */
    void initVectorSpace(
      const RCP<Array<int> >& isBCIndex, 
      bool partitionBCs);
    
    /** */
    void initImporter();

    /** */
    int setupVerb_;

    /** */
    RCP<DOFMapBase> map_;

    /** */
    Mesh mesh_;

    /** */
    Array<CellFilter> subdomains_;

    /** */
    BasisArray basis_;

    /** */
    VectorSpace<double> vecSpace_;

    /** */
    VectorType<double> vecType_;

    /** */
    RCP<GhostImporter<double> > ghostImporter_;

    /** Transformation builder in case when it is needed*/
    RCP<DiscreteSpaceTransfBuilder> transformationBuilder_;

  };

}



#endif
