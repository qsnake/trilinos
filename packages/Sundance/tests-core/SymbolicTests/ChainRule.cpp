#include "SundanceChainRuleEvaluator.hpp"
#include "SundanceUnknownFunctionStub.hpp"
#include "SundanceCombinatorialUtils.hpp"

using namespace Teuchos;
using std::cout;
using std::exception;
using Sundance::List;


void write(const MultipleDeriv& md,
           const Array<MultiSet<int> >& K,
           const Array<MultipleDeriv>& L);

int main(int argc, char** argv)
{
  typedef Array<OrderedPair<Array<MultiSet<int> >, Array<MultipleDeriv> > > CR;
  try
		{
      GlobalMPISession session(&argc, &argv);

			Expr u = new UnknownFunctionStub("u");
			Expr v = new UnknownFunctionStub("v");
			Expr w = new UnknownFunctionStub("w");

      int nArgs = 2;
      MultipleDeriv md = makeDeriv(u,u);
      int order = md.order();
      

      for (int l=1; l<=order; l++)
        {
          Array<int> s(l, nArgs);
          Array<Array<int> > distinctQ = indexCombinations(s);
          Set<MultiSet<int> > q;
          for (int p=0; p<distinctQ.size(); p++)
            {
              q.put(makeMultiSet(distinctQ[p]));
            }
          if (l > 1) cout << " + " << std::endl;
          for (Set<MultiSet<int> >::const_iterator 
                 i=q.begin(); i!=q.end(); i++)
            {
              const MultiSet<int>& lambda = *i;
              if (lambda != *(q.begin())) cout << " + " << std::endl;
              cout << "f_" << lambda << " * [";
              for (int s=1; s<=md.order(); s++)
                {
                  CR p = chainRuleTerms(s, lambda, md);
                  bool firstTerm = true;
                  for (CR::const_iterator j=p.begin(); j!=p.end(); j++)
                    {
                      if (!firstTerm) cout << "+";
                      firstTerm = false;
                      Array<MultiSet<int> > K = j->first();
                      Array<MultipleDeriv> L = j->second();
                      write(md, K, L);
                    }
                }
              cout << "]" << std::endl;
            }
        }
    }
	catch(std::exception& e)
		{
			Out::println(e.what());
		}
}


void write(const MultipleDeriv& md,
           const Array<MultiSet<int> >& K,
           const Array<MultipleDeriv>& L)
{
  int factor = chainRuleMultiplicity(md, K, L);
  if (factor != 1) cout << factor << "*";
  bool firstTerm = true;
  
  for (int j=0; j<K.size(); j++)
    {
      if (!firstTerm) cout << "*";
      firstTerm = false;
      bool firstFactor = true;
      for (MultiSet<int>::const_iterator k=K[j].begin(); k!=K[j].end(); k++)
        {
          if (!firstFactor) cout << "*";
          firstFactor = false;
          int q = *k;
          cout << "D[q_" << q << ", " << L[j] << "]";
        }
    }
}
