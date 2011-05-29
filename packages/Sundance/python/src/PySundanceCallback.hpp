#ifndef PYSUNDANCE_Callback_h
#define PYSUNDANCE_Callback_h

#include "Python.h"


class PySundanceCallback
{
public:
  PySundanceCallback();
  ~PySundanceCallback();

  PyObject* setFunction(PyObject* pyMethod);
  PyObject* getFunction();
  const PyObject* getFunction() const;

private:
  // Private and not implemented so never can be called
  PySundanceCallback(const PySundanceCallback & a_ref);
  const PySundanceCallback & operator = (const PySundanceCallback & a_ref); 

  PyObject* callback_;
};

#endif //Callback_h
