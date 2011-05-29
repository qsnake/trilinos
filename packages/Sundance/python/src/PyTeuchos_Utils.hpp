#ifndef PYTEUCHOS_UTILS_H
#define PYTEUCHOS_UTILS_H

#include "Teuchos_ParameterList.hpp"
#include "Python.h"

// Creates a newly allocated Teuchos parameter list from the input
// object, which must be a Python dictionary.
//
// "bool", "int", "double" and "string" are automatically recognized.
// Other types can be defined here as tuples. For example, the Python
// dictionary can be something like:
// List = {
//   "double parameter": 12.0,
//   "int parameter"   : 12,
//   "string parameter": "12"
// }
//
// \author Marzio Sala, SNL 9215
//
// \date Last modified on 08-Aug-05

Teuchos::ParameterList dict2ParameterList(PyObject* obj);



#endif
