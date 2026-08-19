/// @file BindPost.cpp
/// @brief Bindings for post processing (Post/)

#include "Bindings.h"

#include "TPZVTKGenerator.h"
#include "pzcmesh.h"

void InitPost(py::module_ &m) {
  /// only the simplest of the four constructors is exposed
  py::class_<TPZVTKGenerator>(m, "TPZVTKGenerator", "Writes the solution of a mesh as .vtk files.")
      .def(py::init<TPZCompMesh *, const TPZVec<std::string> &, std::string, int>(), py::arg("cmesh"), py::arg("fields"), py::arg("filename"), py::arg("vtkres") = 0, py::keep_alive<1, 2>(), "Field names are resolved by the material at runtime; the filename "
                                                                                                                                                                                              "carries no extension.")
      .def("Do", &TPZVTKGenerator::Do, "Writes one output step.")
      .def("SetNThreads", &TPZVTKGenerator::SetNThreads, py::arg("nthreads"));
}
