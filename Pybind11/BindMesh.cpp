/// @file BindMesh.cpp
/// @brief Bindings for the computational meshes (Mesh/)

#include "Bindings.h"

#include "pzcmesh.h"
#include "TPZMultiphysicsCompMesh.h"

void InitMesh(py::module_ &m)
{
     /// the default constructor is protected, so the geometric mesh is required
     py::class_<TPZCompMesh>(m, "TPZCompMesh",
                             "Computational mesh built on a geometric mesh.")
         .def(py::init<TPZGeoMesh *>(), py::arg("gmesh"), py::keep_alive<1, 2>())
         .def("NElements", &TPZCompMesh::NElements)
         .def("NConnects", &TPZCompMesh::NConnects)
         .def("NEquations", &TPZCompMesh::NEquations)
         .def("Dimension", &TPZCompMesh::Dimension)
         .def("Reference", &TPZCompMesh::Reference,
              py::return_value_policy::reference,
              "Geometric mesh this mesh refers to but does not own.")
         .def("SetDefaultOrder", &TPZCompMesh::SetDefaultOrder, py::arg("order"))
         .def("__str__", [](const TPZCompMesh &cmesh)
              { return PrintToString(cmesh); });

     py::class_<TPZMultiphysicsCompMesh, TPZCompMesh>(
         m, "TPZMultiphysicsCompMesh",
         "Mesh combining several atomic approximation spaces.")
         .def(py::init<>())
         .def(py::init<TPZGeoMesh *>(), py::arg("gmesh"), py::keep_alive<1, 2>());
}
