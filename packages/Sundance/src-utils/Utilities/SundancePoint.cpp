#include "SundancePoint.hpp"

using namespace Teuchos;
using namespace Sundance;



void Point::boundsCheck(int i) const
{
  TEST_FOR_EXCEPTION(i < 0 || i>=dim_, RuntimeError,
                     "Point::boundsCheck() dim:" << dim_ << " , i:" << i);
}

namespace std
{
  ostream& operator<<(std::ostream& os, const Point& point)
  {
    os << "(";
    os.precision(12);
    for (int i=0; i<point.dim(); i++)
      {
        
        os << point[i];
        if (i<(point.dim()-1)) os << ", ";
      }
    os << ")";
    return os;
  }
}

