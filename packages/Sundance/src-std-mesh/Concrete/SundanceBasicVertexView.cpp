#include "SundanceBasicVertexView.hpp"
#include "Teuchos_Utils.hpp"

using namespace Teuchos;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;

string VertexView::toString() const
{
  int* ptr = *base_ +  offset_*length_;
	string rtn="{";
	for (int i=0; i<length_; i++) 
		{
			rtn += Teuchos::toString(ptr[i]);
			if (i < length_-1) rtn += ", ";
		}
	rtn += "}";
	return rtn;
}

/*
 * Return a hash code for the vertex set. 
 */
int VertexView::hashCode() const
{
  int rtn = 0;
      int* p = *base_ + offset_*length_;

      for (int i=0; i<length_; i++)
        {
          rtn += p[i];
        }

      return rtn;

#ifdef BLARF
      // fails with sign int
  int* p = *base_ + offset_*length_;

  unsigned char* key = reinterpret_cast<unsigned char*>(p);
  unsigned int key_len = length_ * sizeof(int);

  /* Jenkins hash */

  unsigned int hash = 0;
  int i;
  
  for (i = 0; i < key_len; i++) {
    hash += key[i];
    hash += (hash << 10);
    hash ^= (hash >> 6);
  }
  hash += (hash << 3);
  hash ^= (hash >> 11);
  hash += (hash << 15);
  return hash;
#endif
}
