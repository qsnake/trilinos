#include "SundanceIntHashSet.hpp"
#include "SundanceSet.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_Time.hpp"

using namespace Sundance;
using namespace Teuchos;

int main(int argc, char** argv)
{
  try
    {
      GlobalMPISession session(&argc, &argv);
      
      int nRow = 40000;
      int nReps = 4;
      
      Time tSet("STL set time");
      Time tHash("hash set time");
      
      for (int nData=20; nData<360; nData+=20)
      {
        {
        tSet.start();
        Array<Set<int> > sets(nRow);


        for (int r=0; r<nRow; r++)
          {
            Set<int>& s = sets[r];
            for (int rep=0; rep<nReps; rep++)
              {
                for (int x=0; x<nData; x++)
                  {
                    int y = rand() % nData;
                    s.put(y);
                  }
              }
          }
        tSet.stop();
        }
        {
          tHash.start();
          Array<IntHashSet> sets(nRow+1);
          for (int r=0; r<nRow; r++)
            {
              sets[r].setCapacity(nData+1);
            }


          for (int r=0; r<nRow; r++)
            {
            IntHashSet& s = sets[r];
            for (int rep=0; rep<nReps; rep++)
              {
                for (int x=0; x<nData; x++)
                  {
                    int y = rand() % nData;
                    s.put(y);
                  }
              }
          }
        tHash.stop();
        }

        
        std::cerr << nData << "\t set=" << tSet.totalElapsedTime()
             << "\t hash=" << tHash.totalElapsedTime() 
             << "\t ratio=" 
             << tHash.totalElapsedTime()/tSet.totalElapsedTime()
             << std::endl;
      }

      
    }
  catch(std::exception& e)
    {
      std::cerr << "caught exception " << e.what() << std::endl;
    }
  
}
