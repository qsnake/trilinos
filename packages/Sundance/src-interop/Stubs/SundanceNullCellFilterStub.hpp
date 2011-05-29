#ifndef SUNDANCE_NULLCELLFILTER_STUB_H
#define SUNDANCE_NULLCELLFILTER_STUB_H


#include "SundanceDefs.hpp"
#include "SundanceCellFilterStub.hpp"



namespace Sundance
{
using namespace Teuchos;
  
  
/** 
 *
 * <h4> Notes for framework interface implementors </h4>
 *
 *  
 */
class NullCellFilterStub : public CellFilterStub
{
public:
  /** Empty ctor */
  NullCellFilterStub();

  /** virtual dtor */
  virtual ~NullCellFilterStub(){;}

  /** Write to XML */
  virtual XMLObject toXML() const ;

  /** Ordering for storage in STL maps */
  virtual bool lessThan(const CellFilterStub* other) const ;

  /** \name Printable interface */
  //@{
  /** Print to a stream */
  virtual void print(std::ostream& os) const {os << toXML();}
  //@}

  /** \name Describable interface */
  //@{
  /** Print to a stream */
  virtual std::string description() const 
    {return "NullCellFilterStub";}
  //@}

  /** */
  virtual RCP<CellFilterStub> getRcp() {return rcp(this);}

};

}


#endif
