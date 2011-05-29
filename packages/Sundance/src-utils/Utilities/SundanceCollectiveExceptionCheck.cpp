#include "SundanceCollectiveExceptionCheck.hpp"
#include "Teuchos_MPIComm.hpp"

using namespace Sundance;
using namespace Teuchos;

namespace Sundance
{
  void reportFailure(const MPIComm& comm)
  {
    int myBad = 1;
    int anyBad = 0;
    comm.allReduce((void*) &myBad, (void*) &anyBad, 1, MPIComm::INT,
                     MPIComm::SUM);
  }

  bool checkForFailures(const MPIComm& comm)
  {
    int myBad = 0;
    int anyBad = 0;
    comm.allReduce((void*) &myBad, (void*) &anyBad, 1, MPIComm::INT,
                     MPIComm::SUM);
    return anyBad > 0;
  }


}



	





