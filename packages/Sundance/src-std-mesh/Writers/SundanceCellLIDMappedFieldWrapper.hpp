#ifndef SUNDANCE_CELLLID_MAPPED_FIELDWRAPPER_H
#define SUNDANCE_CELLLID_MAPPED_FIELDWRAPPER_H

#include "SundanceFieldBase.hpp"

namespace Sundance
{

using Teuchos::Array;
using Teuchos::RCP;

class CellLIDMappedFieldWrapper : public FieldBase
{
public:
  /** */
  CellLIDMappedFieldWrapper(int cellDim,
    int nFuncs, const RCP<Array<double> >& data)
    : cellDim_(cellDim), nFuncs_(nFuncs), data_(data) {}

  /** virtual dtor */
  virtual ~CellLIDMappedFieldWrapper(){}

  /** */
  double getData(int cellDim, int cellID, int elem) const 
    {
      TEST_FOR_EXCEPT(cellDim != cellDim_);
      return (*data_)[cellID*nFuncs_ + elem];
    }
  
  /** */
  bool isDefined(int cellDim, int cellID, int elem) const 
    {return true;}
  
  /** */
  int numElems() const {return nFuncs_;}
  
  /** */
  bool isPointData() const {return cellDim_==0;}
  
  /* */
  GET_RCP(FieldBase);

private:
  int cellDim_;
  int nFuncs_;
  RCP<Array<double> > data_;
};



}

#endif
