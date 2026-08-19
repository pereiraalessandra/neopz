/// @file BindGeom.cpp
/// @brief Bindings for the geometric mesh and its generators (Geom/, Pre/)

#include "Bindings.h"

#include "MMeshType.h"
#include "TPZGenGrid2D.h"
#include "TPZGenGrid3D.h"
#include "TPZGmshReader.h"
#include "pzcheckgeom.h"
#include "pzgmesh.h"

void InitGeom(py::module_ &m) {
  py::class_<TPZGeoMesh>(m, "TPZGeoMesh")
      .def(py::init<>())
      .def("NElements", &TPZGeoMesh::NElements, "Includes boundary elements.")
      .def("NNodes", &TPZGeoMesh::NNodes)
      .def("Dimension", &TPZGeoMesh::Dimension)
      .def("SetDimension", &TPZGeoMesh::SetDimension, py::arg("dim"))
      .def("SetName", &TPZGeoMesh::SetName, py::arg("name"))
      .def("BuildConnectivity", &TPZGeoMesh::BuildConnectivity, "Call after the elements are created.")
      .def("__str__", [](const TPZGeoMesh &gmesh) { return PrintToString(gmesh); });

  py::class_<TPZGenGrid2D>(m, "TPZGenGrid2D")
      .def(py::init<const TPZVec<int> &, const TPZVec<REAL> &, const TPZVec<REAL> &>(), py::arg("nx"), py::arg("x0"), py::arg("x1"), "nx divisions per direction; x0 and x1 are opposite corners in 3D.")
      .def("SetElementType", &TPZGenGrid2D::SetElementType, py::arg("type"))
      .def("Read", py::overload_cast<TPZGeoMesh *, int>(&TPZGenGrid2D::Read), py::arg("gmesh"), py::arg("matid") = 1)
      .def("SetBC", py::overload_cast<TPZGeoMesh *, int, int>(&TPZGenGrid2D::SetBC), py::arg("gmesh"), py::arg("side"), py::arg("bcid"), "Sides 4 to 7 are the edges of the global rectangle.");

  py::class_<TPZGenGrid3D>(m, "TPZGenGrid3D")
      .def(py::init<const TPZVec<REAL> &, const TPZVec<REAL> &, const TPZVec<int> &, MMeshType>(), py::arg("minX"), py::arg("maxX"), py::arg("nelDiv"), py::arg("elType"))
      // Both return the same internally-owned mesh. BuildVolumetricElements
      // uses reference to avoid a second owning wrapper; BuildBoundaryElements
      // keeps the default owning policy, the one a script actually keeps.
      .def("BuildVolumetricElements", &TPZGenGrid3D::BuildVolumetricElements, py::arg("matIdDomain"), py::return_value_policy::reference)
      .def("BuildBoundaryElements", &TPZGenGrid3D::BuildBoundaryElements, py::arg("matIdZmin"), py::arg("matIdXmin"), py::arg("matIdYmin"), py::arg("matIdXmax"), py::arg("matIdYmax"), py::arg("matIdZmax"));

  py::class_<TPZGmshReader>(m, "TPZGmshReader", "Reads geometric meshes from Gmsh .msh files "
                                                "(format 3 or 4, auto-detected).")
      .def(py::init<>())
      // Passing an existing mesh would leave the same pointer owned by two
      // Python wrappers, double-freeing it at collection; reference is used
      // only then. gmesh=None still returns an owning wrapper, unchanged.
      .def("GeometricGmshMesh", [](TPZGmshReader &self, const std::string &file_name, TPZGeoMesh *gmesh, bool addNonAssignedEls) {
                  bool owning = (gmesh == nullptr);
                  TPZGeoMesh *result = self.GeometricGmshMesh(file_name, gmesh, addNonAssignedEls);
                  return py::cast(result, owning ? py::return_value_policy::take_ownership
                                                  : py::return_value_policy::reference); }, py::arg("file_name"), py::arg("gmesh") = static_cast<TPZGeoMesh *>(nullptr), py::arg("addNonAssignedEls") = true, "Reads the file into a newly created mesh (or gmesh, if given); "
                                                                                                                                                                                                                                                        "BuildConnectivity is called internally, unlike TPZGenGrid2D/3D.")
      .def("SetCharacteristiclength", &TPZGmshReader::SetCharacteristiclength, py::arg("length"))
      .def("SetVerbose", &TPZGmshReader::SetVerbose, py::arg("verbose"))
      .def("GetDimPhysicalTagName", &TPZGmshReader::GetDimPhysicalTagName, "Physical group names, indexed by dimension then material id, "
                                                                           "as declared in the file.");

  // Mostly a mesh-consistency checker; only the refinement entry point is
  // bound. Does not own the mesh (no destructor deletes fMesh); keep_alive
  // is used instead of transferring ownership.
  py::class_<TPZCheckGeom>(m, "TPZCheckGeom")
      .def(py::init<TPZGeoMesh *>(), py::arg("gmesh"), py::keep_alive<1, 2>())
      .def("UniformRefine", &TPZCheckGeom::UniformRefine, py::arg("nDiv"), "Divides every element, domain and boundary alike, nDiv times.");
}
