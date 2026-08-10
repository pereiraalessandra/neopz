/// @file BindGeom.cpp
/// @brief Bindings for the geometric mesh and its generators (Geom/, Pre/)

#include "Bindings.h"

#include "pzgmesh.h"
#include "TPZGenGrid2D.h"
#include "MMeshType.h"

void InitGeom(py::module_ &m)
{
     py::class_<TPZGeoMesh>(m, "TPZGeoMesh",
                            "Nodes, elements and their connectivity.")
         .def(py::init<>())
         .def("NElements", &TPZGeoMesh::NElements, "Includes boundary elements.")
         .def("NNodes", &TPZGeoMesh::NNodes)
         .def("Dimension", &TPZGeoMesh::Dimension)
         .def("SetDimension", &TPZGeoMesh::SetDimension, py::arg("dim"))
         .def("SetName", &TPZGeoMesh::SetName, py::arg("name"))
         .def("BuildConnectivity", &TPZGeoMesh::BuildConnectivity,
              "Call after the elements are created.")
         .def("__str__", [](const TPZGeoMesh &gmesh)
              { return PrintToString(gmesh); });

     py::class_<TPZGenGrid2D>(m, "TPZGenGrid2D",
                              "Structured grid generator for 2D domains.")
         .def(py::init<const TPZVec<int> &, const TPZVec<REAL> &, const TPZVec<REAL> &>(),
              py::arg("nx"), py::arg("x0"), py::arg("x1"),
              "nx divisions per direction; x0 and x1 are opposite corners in 3D.")
         .def("SetElementType", &TPZGenGrid2D::SetElementType, py::arg("type"))
         .def("Read", py::overload_cast<TPZGeoMesh *, int>(&TPZGenGrid2D::Read),
              py::arg("gmesh"), py::arg("matid") = 1,
              "Fills the geometric mesh with the generated elements.")
         .def("SetBC", py::overload_cast<TPZGeoMesh *, int, int>(&TPZGenGrid2D::SetBC),
              py::arg("gmesh"), py::arg("side"), py::arg("bcid"),
              "Sides 4 to 7 are the edges of the global rectangle.");
}
