#include "SundanceCombinatorialUtils.hpp"
#include "Teuchos_GlobalMPISession.hpp"


using namespace Sundance;
using namespace Teuchos;



#define TEST_MS(x) \
  {\
    Set<MultiSet<int> > subs = multisetSubsets(x);\
    write(x, subs);                                                \
    Set<MultiSet<MultiSet<int> > > parts = multisetPartitions(x);\
    write(x, parts);                                                 \
    Array<Array<MultiSet<int> > > comps = multisetCompositions(x);\
    write(x, comps);                                               \
  }



void write(const MultiSet<int>& x, 
           const Set<MultiSet<MultiSet<int> > >& y)
{
  std::cout << "---- Partitions of " << x << " ------------------"
       << std::endl;
  for (Set<MultiSet<MultiSet<int> > >::const_iterator 
         i=y.begin(); i!=y.end(); i++)
    {
      std::cout << *i << std::endl;
    }
}

void write(const MultiSet<int>& x, 
           const Array<Array<MultiSet<int> > >& y)
{
  std::cout << "---- Compositions of " << x << " ------------------"
       << std::endl;
  for (int i=0; i<y.size(); i++)
    {
      std::cout << y[i] << std::endl;
    }
}

void write(const MultiSet<int>& x, 
           const Set<MultiSet<int> >& y)
{
  std::cout << "---- Subsets of " << x << " ------------------"
       << std::endl;
  for (Set<MultiSet<int> >::const_iterator 
         i=y.begin(); i!=y.end(); i++)
    {
      std::cout << *i << std::endl;
    }
}



int main(int argc, char** argv)
{
  int stat = 0;
  try
		{
      GlobalMPISession session(&argc, &argv);


      bool bad = false;

      for (int n=1; n<=4; n++)
        {
          Array<Array<Array<int> > > c = compositions(n);
          std::cout << "N=" << n << " compositions=" << c << std::endl;

          MultiSet<int> mu;
          for (int m=1; m<=n; m++)
            {
              mu.put(m);

            }
          for (int m=1; m<=n; m++)
            {
              Array<Array<Array<int> > > b = binnings(mu, m);
              std::cout << "binnings = " << b << std::endl;
            }

          std::cout << "--------- non-neg compositions" << std::endl;
          for (int m=1; m<=n; m++)
            {
              for (int k=1; k<=n; k++)
                {
                  Array<Array<int> > a = nonNegCompositions(m, k);
                  std::cout << m << " " << k << " " << std::endl;
                  for (int l=0; l<a.size(); l++)
                    {
                      std::cout << "         " << a[l] << std::endl;
                    }
                }
            }
          
          std::cout << "-------- index combs ---- " << std::endl;
          Array<int> s = tuple(2,3,2);
          Array<Array<int> > C = indexCombinations(s);
          for (int m=0; m<C.size(); m++)
            {
              std::cout << C[m] << std::endl;
            }
        }

      std::cout << "--------- index tuples ----------------" << std::endl;

      Array<Array<int> > x = distinctIndexTuples(2, 6);

      std::cout << "num choices = " << x.size() << std::endl;

      for (int i=0; i<x.size(); i++) 
        {
          if ((i % 5)==0) std::cout << std::endl;
          std::cout << x[i] << std::endl;
        }
      
#ifdef BLAH
      TEST_MS(makeMultiSet(1));
      TEST_MS(makeMultiSet(1, 1));
      TEST_MS(makeMultiSet(1, 2));
      TEST_MS(makeMultiSet(1, 1, 2));
      TEST_MS(makeMultiSet(1, 1, 2, 2));
      TEST_MS(makeMultiSet(1, 2, 3));
      TEST_MS(makeMultiSet(1, 2, 2, 3));
      TEST_MS(makeMultiSet(1, 2, 2, 3, 3));
      TEST_MS(makeMultiSet(1, 2, 2, 3, 3, 3));
#endif // BLAH

      if (!bad) 
        {
          std::cerr << "all tests PASSED" << std::endl;
        }
      else
        {
          stat = -1;
          std::cerr << "a test has FAILED" << std::endl;
        }
    }
	catch(std::exception& e)
		{
      stat = -1;
      std::cerr << "detected exception " << e.what() << std::endl;
		}

  return stat;
}
