#ifndef SUNDANCERIVARAELEMENTITERATOR_H
#define SUNDANCERIVARAELEMENTITERATOR_H

#include "SundanceRivaraElement.hpp"

namespace Sundance
{
  namespace Rivara
  {
    class RivaraMesh;

    class ElementIterator
    {
    public:
      ElementIterator(const RivaraMesh* mesh);

      bool hasMoreElements() const ;

      const Element* getNextElement() const ;

    private:
      RivaraMesh* mesh_;

      mutable int rootIndex_;

      mutable Element* current_;

      mutable bool startingNewTree_;
    };
  }
}

#endif
