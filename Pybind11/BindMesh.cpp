/// @file BindMesh.cpp
/// @brief Bindings for the computational meshes (Mesh/)

#include "Bindings.h"

#include "TPZMultiphysicsCompMesh.h"
#include "pzbuildmultiphysicsmesh.h"
#include "pzcmesh.h"

void InitMesh(py::module_ &m) {
  /// the default constructor is protected; the geometric mesh is required
  py::class_<TPZCompMesh>(m, "TPZCompMesh")
      .def(py::init<TPZGeoMesh *>(), py::arg("gmesh"), py::keep_alive<1, 2>())
      .def("NElements", &TPZCompMesh::NElements)
      .def("NConnects", &TPZCompMesh::NConnects)
      .def("NEquations", &TPZCompMesh::NEquations)
      .def("Dimension", &TPZCompMesh::Dimension)
      .def("Reference", &TPZCompMesh::Reference, py::return_value_policy::reference, "Geometric mesh this mesh refers to but does not own.")
      .def("SetDefaultOrder", &TPZCompMesh::SetDefaultOrder, py::arg("order"))
      .def("LoadReferences", &TPZCompMesh::LoadReferences, "Rebuilds the CompEl<->GeoEl mapping; call before EvaluateError.")
      // C++ is EvaluateError(bool store_error, TPZVec<REAL>& errorSum), an
      // out-param Python can't mutate the same way; nErrors stands in for
      // errorSum's size, and the filled result comes back as the return
      // value. storeError keeps its C++ position and has no default, matching C++.
      .def("EvaluateError", [](TPZCompMesh &self, bool storeError, int nErrors) {
                  TPZManVector<REAL> errorSum(nErrors, 0.0);
                  self.EvaluateError(storeError, errorSum);
                  return errorSum; }, py::arg("storeError"), py::arg("nErrors"), "One norm per component reported by the domain material's "
                                                                                                                                "NEvalErrors().")
      .def("__str__", [](const TPZCompMesh &cmesh) { return PrintToString(cmesh); });

  py::class_<TPZMultiphysicsCompMesh, TPZCompMesh>(
      m, "TPZMultiphysicsCompMesh", "Mesh combining several atomic approximation spaces.")
      .def(py::init<>())
      .def(py::init<TPZGeoMesh *>(), py::arg("gmesh"), py::keep_alive<1, 2>())
      .def("MeshVector", &TPZMultiphysicsCompMesh::MeshVector, py::return_value_policy::reference, "Atomic meshes this mesh was built from.");

  // Free function, not a method: TPZBuildMultiphysicsMesh has no state,
  // only this static utility. Copies the multiphysics solution back into
  // the atomic meshes; without it, EvaluateError reads a stale (empty)
  // solution on the atomic meshes and reports a spurious error.
  m.def("TransferFromMultiPhysics", py::overload_cast<const TPZVec<TPZCompMesh *> &, TPZCompMesh *>(&TPZBuildMultiphysicsMesh::TransferFromMultiPhysics), py::arg("meshvec"), py::arg("cmesh"), "Call after solving and before EvaluateError or reading atomic "
                                                                                                                                                                                                "mesh fields.");
}
