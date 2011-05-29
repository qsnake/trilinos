#ifndef PYSUNDANCE_FIATSCALARADAPTER_H
#define PYSUNDANCE_FIATSCALARADAPTER_H
#include "Python.h"
#include "SundanceBasisFamilyBase.hpp"
#include "SundanceDefs.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include <stack>

namespace Sundance
{
  using namespace Sundance;
  using namespace Sundance;
  using namespace Sundance;
  
  using namespace Sundance;
  using namespace Sundance;
  using namespace Teuchos;


  class FIATScalarAdapter : public ScalarBasis
  {
  public:
    FIATScalarAdapter( PyObject *pyfamilyclass , int order );
    
    ~FIATScalarAdapter(); 

    bool supportsCellTypePair(
      const CellType& maximalCellType,
      const CellType& cellType
      ) const ;
      
    void getReferenceDOFs(
      const CellType& maximalCellType,
      const CellType& cellType,
      Array<Array<Array<int> > >& dofs) const;
    
    int nReferenceDOFs(
      const CellType& maximalCellType,
      const CellType& cellType
      ) const ;
    
    void refEval(
      const CellType& maximalCellType,
      const CellType& cellType,
      const Array<Point>& pts,
      const MultiIndex& deriv,
      Array<Array<Array<double> > >& result) const ;
    
    int order() const { return order_; }

    /** Needed for Printable interface */
    void print(std::ostream& os) const {os << "FIATScalarAdapter";}

    /* Needed for Handleable interface */
    GET_RCP(BasisFamilyBase);
    
  private:
    // we instantiate the basis for each order
    int order_;
    Array<PyObject *> bases_;
    Array<Array<Array<Array<int> > > > dof_;
  };
  
  
}
#endif
