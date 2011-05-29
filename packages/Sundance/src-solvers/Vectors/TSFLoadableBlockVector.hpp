/* @HEADER@ */
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
/* @HEADER@ */

#ifndef TSFLOADABLEBLOCKVECTOR_HPP
#define TSFLOADABLEBLOCKVECTOR_HPP

#include "SundanceDefs.hpp"
#include "Thyra_VectorBase.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY

namespace TSFExtended
{
  
  /**
   * LoadableBlockVector provides a LoadableVector interface to a 
   * physically-partitioned block 1x2 vector, making it appear to the
   * fill routine as if the block vector is a single vector. This 
   * is intended for filling systems where the internal and BC equations
   * and unknowns are stored in physically separate blocks.
   *
   * @author Kevin Long (krlong@sandia.gov)
   */
  class LoadableBlockVector : public LoadableVector<double>
    {
    public:
      /** */
      LoadableBlockVector(
        const Vector<double>& vec,
        int lowestLocalRow,
        int highestLocalRow,
        const RCP<Array<int> >& isBCRow
        ) : lowestLocalRow_(lowestLocalRow), 
            highestLocalRow_(highestLocalRow),
            bcVec_(),
            internalVec_(),
            isBCRow_(isBCRow)
        {
          TEST_FOR_EXCEPTION(vec.space().numBlocks() != 2,
            std::runtime_error, "LoadableBlockVector expected numBlocks=2, "
            "found " << vec.space().numBlocks());

          internalVec_ = vec.getBlock(0);
          bcVec_ = vec.getBlock(1);
        }


      /** virtual dtor */
      virtual ~LoadableBlockVector() {;}

      /** set a single element at the given global index */
      void setElement(OrdType globalIndex, const double& value) 
        {
          if (globalIndex < lowestLocalRow_ || globalIndex >= highestLocalRow_)
          {
            int localRow = globalIndex-lowestLocalRow_;
            if ((*isBCRow_)[localRow])
            {
              bcVec_.setElement(globalIndex, value);
            }
            else
            {
              internalVec_.setElement(globalIndex, value);
            }
          }
        }

      /** add to the existing value of 
       * a single element at the given global index */
      void addToElement(OrdType globalIndex, const double& value) 
        {
          if (globalIndex >= lowestLocalRow_ && globalIndex < highestLocalRow_)
          {
            int localRow = globalIndex-lowestLocalRow_;
            if ((*isBCRow_)[localRow])
            {
              bcVec_.addToElement(globalIndex, value);
            }
            else
            {
              internalVec_.addToElement(globalIndex, value);
            }
          }
        }

      const Vector<double>& bcBlock() const {return bcVec_;}
      const Vector<double>& internalBlock() const {return internalVec_;}
      
    private:
      int lowestLocalRow_;
      int highestLocalRow_;
      Vector<double> bcVec_;
      Vector<double> internalVec_;
      RCP<Array<int> > isBCRow_;
  };
  
}

#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
