#include "PySundanceCallback.hpp"

#include <iostream>

using std::cout;
using std::endl;

PySundanceCallback::PySundanceCallback()
  : callback_(NULL)
{
  // Nothing here yet
}

PySundanceCallback::~PySundanceCallback()
{
  Py_XDECREF(callback_);   /* Dispose of callback */
  callback_ = 0;
}

PyObject * PySundanceCallback::setFunction( PyObject * pyMethod)
{
    PyObject *p_result = NULL;
    PyObject *p_temp   = NULL;
    PyObject *p_func   = NULL;


    assert(0 != pyMethod && "Null argument passed to setFunction()");
    
    if (PyArg_ParseTuple(pyMethod, "O", &p_temp)) {
      // Assume that if this is a tuple, the item in the tuple is a
      // PyObject that is a pointer to a Python function. This is the
      // case when this function is called from Python.
      p_func = p_temp;
    } else {
      // Otherwise we assume that this function is directly passed a
      // PyObject that is a pointer to a Python function.  This is the
      // case when this function is called from C++.
      p_func = pyMethod;
    }

    if (!PyCallable_Check(p_func)) {
      PyErr_SetString(PyExc_TypeError,
		      "Function parameter must be callable");
      cout << "PyObject passed to function is not callable" << std::endl ;
      return NULL;
    }
    Py_XINCREF(p_func);          /* Add a reference to new callback */
    Py_XDECREF(callback_);     /* Dispose of previous callback    */
    callback_ = p_func;        /* Remember new callback           */
    /* Boilerplate to return "None" */
    Py_INCREF(Py_None);
    p_result = Py_None;

    assert(0 != callback_ && "Pointer to callback not set");
    return p_result;
}
 
PyObject * PySundanceCallback::getFunction()
{
  assert (0 != callback_ && "PySundanceCallback function not yet assigned");
  return callback_;
}

const PyObject * PySundanceCallback::getFunction() const
{
  assert (0 != callback_ && "PySundanceCallback function not yet assigned");
  return callback_;
}
