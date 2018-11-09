#include <Python.h>

#include <iostream>

static PyObject*
intersection_c(PyObject* self, PyObject* args)
{
  // Parse the iput arguments.
  float r1, r2, u1, u2;
  PyObject* surface;
  if (!PyArg_ParseTuple(args, "ffffO", &r1, &r2, &u1, &u2, &surface)) return nullptr;

  // Size of list
  // Py_ssize_t n_regions = PyList_Size(arg_list);
  std::cout << r1 << "\n";
  std::cout << r2 << "\n";
  std::cout << u1 << "\n";
  std::cout << u2 << "\n";

  PyObject* py_surf_type = PyObject_GetAttrString(surface, "typenum");
  long surf_type = PyLong_AsLong(py_surf_type);

  std::cout << surf_type << "\n";
  

//  for (auto i = 0; i < n_regions; i++) {
//    // Grab the region from the list.  Note that this returns a borrowed
//    // reference.  That means the refcount was not changed.
//    PyObject* region = PyList_GetItem(arg_list, i);
//
//    // This reterns a new reference so we have to decrease the refcount when
//    // we are done.
//    PyObject* material = PyObject_GetAttrString(region, "material");
//
//    // Do sommething with the material.
//    long mat = PyLong_AsLong(material);
//    std::cout << "region_list[" << i << "].material = " << mat << "\n";
//
//    // We're now done with the material so decrease its refcount.
//    Py_DECREF(material);
//  }
//
//  // All done.  In this case we want to return None.  None is a special object
//  // in that there is really only one None object and any Python session just
//  // has multiple copies of it floating around.  So to return None, we increase
//  // the refcount of that one special object and return a pointer to it.
//  Py_INCREF(Py_None);
  return Py_None;
} //function

//==============================================================================
// Boilerplate
//==============================================================================

static PyMethodDef MyExtensionMethods[] = {
  {"intersection_c", intersection_c, METH_VARARGS,
    "Calculate the intersection given a coordinate angle and surface"},
  // More methods go here.  They have the form
  // {"name_that_python_sees", name_of_the_C_function, METH_VARARGS, "comment"},
  {NULL, NULL, 0, NULL}
};

static struct PyModuleDef intersection_ext_module = {
  PyModuleDef_HEAD_INIT,
  "intersection_ext",
  NULL,
  -1,
  MyExtensionMethods
};

PyMODINIT_FUNC
PyInit_my_extension(void)
{
  PyObject* mod = PyModule_Create(&intersection_ext_module);
  return mod;
};
