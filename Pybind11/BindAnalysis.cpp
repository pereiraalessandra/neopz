/// @file BindAnalysis.cpp
/// @brief Bindings for the assembly and solution of the linear system
/// (Analysis/, StrMatrix/, Solvers/)

#include "Bindings.h"
#include <pybind11/numpy.h>

#include "TPZLinearAnalysis.h"
#include "TPZMatrixSolver.h"
#include "TPZSolver.h"
#include "TPZStructMatrix.h"
#include "pzcmesh.h"
#include "pzmatrix.h"
#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"

namespace {
/// Clones the assembled matrix into a direct solver to use as a
/// preconditioner. Valid for the linear problems this binding solves:
/// assembled once, never reassembled; the frozen copy never goes stale.
/// NeoPZ's SetReferenceMatrix + UpdateFrom mechanism syncs a
/// preconditioner across repeated assemblies instead (e.g. Newton), but
/// has no test coverage anywhere in NeoPZ. TPZAutoPointer stays inside
/// this function; only the TPZStepSolver is returned.
TPZStepSolver<STATE> *MakePreconditionerFrom(TPZLinearAnalysis &an, DecomposeType precondDecompose) {
  auto *solver = dynamic_cast<TPZMatrixSolver<STATE> *>(an.Solver());
  if (!solver || !solver->Matrix()) {
    std::cout << "MakePreconditionerFrom: call Assemble() first\n";
    DebugStop();
  }
  TPZAutoPointer<TPZMatrix<STATE>> matClone = solver->Matrix()->Clone();
  auto *precond = new TPZStepSolver<STATE>(matClone);
  precond->SetDirect(precondDecompose);
  return precond;
}
} // namespace

void InitAnalysis(py::module_ &m) {
  /// bare registrations: they only need to exist for the setters below
  py::class_<TPZStructMatrix>(m, "TPZStructMatrix");
  py::class_<TPZSolver>(m, "TPZSolver");
  // Registered so a TPZStepSolver can be passed as the TPZMatrixSolver&
  // SetCG/SetGMRES expect for precond — without this intermediate level,
  // pybind only knows TPZSolver and rejects the call.
  py::class_<TPZMatrixSolver<STATE>, TPZSolver>(m, "TPZMatrixSolver");

  py::class_<TPZSkylineStructMatrix<STATE>, TPZStructMatrix>(
      m, "TPZSkylineStructMatrix", "Structural matrix in skyline format.")
      .def(py::init<TPZCompMesh *>(), py::arg("cmesh"))
      // TPZStrMatParInterface (where SetNumThreads lives) is a virtual
      // base reached through TPar=TPZStructMatrixOR<STATE>: no pointer to
      // member, called directly instead.
      .def("SetNumThreads", [](TPZSkylineStructMatrix<STATE> &self, int n) { self.SetNumThreads(n); }, py::arg("n"));

  py::class_<TPZStepSolver<STATE>, TPZMatrixSolver<STATE>>(m, "TPZStepSolver")
      .def(py::init<>())
      .def("SetDirect", &TPZStepSolver<STATE>::SetDirect, py::arg("decomposeType"), "Mixed formulations are symmetric indefinite: use ELDLt.")
      .def("SetCG", &TPZStepSolver<STATE>::SetCG, py::arg("numIterations"), py::arg("precond"), py::arg("tol"), py::arg("fromCurrent") = 0, "precond is cloned internally, not referenced; a TPZStepSolver "
                                                                                                                                            "configured with SetJacobi or SetDirect works as one. Only "
                                                                                                                                            "guaranteed to converge for a symmetric positive definite system.")
      .def("SetGMRES", &TPZStepSolver<STATE>::SetGMRES, py::arg("numIterations"), py::arg("numVectors"), py::arg("precond"), py::arg("tol"), py::arg("fromCurrent") = 0, "Unlike SetCG, works for symmetric indefinite systems such as "
                                                                                                                                                                         "the mixed formulation's saddle point.")
      .def("SetJacobi", &TPZStepSolver<STATE>::SetJacobi, py::arg("numIterations"), py::arg("tol"), py::arg("fromCurrent") = 0);

  py::class_<TPZLinearAnalysis>(m, "TPZLinearAnalysis")
      .def(py::init<TPZCompMesh *>(), py::arg("cmesh"), py::keep_alive<1, 2>())
      /// no keep_alive: the argument is cloned
      .def("SetStructuralMatrix", py::overload_cast<TPZStructMatrix &>(&TPZLinearAnalysis::SetStructuralMatrix), py::arg("strmatrix"))
      .def("SetSolver", &TPZLinearAnalysis::SetSolver, py::arg("solver"))
      .def("Solver", &TPZLinearAnalysis::Solver, py::return_value_policy::reference, "The current solver, cloned when SetSolver was called; needed "
                                                                                     "to reach its assembled matrix (e.g. via MakePreconditionerFrom).")
      .def("Assemble", &TPZLinearAnalysis::Assemble)
      .def("Solve", &TPZLinearAnalysis::Solve, "Loads the solution into the mesh.")
      .def("SolutionNorm", [](TPZLinearAnalysis &self) { return Norm(self.Solution()); })
      // Zero-copy view: TPZFMatrix stores column-major (fElem[col*rows+row]),
      // so the strides below are the transpose of what a C-order array would
      // use. base=self keeps the analysis (and the buffer it owns) alive for
      // as long as the returned array is.
      .def("Solution", [](py::object self) {
                  TPZFMatrix<STATE> &sol = self.cast<TPZLinearAnalysis &>().Solution();
                  return py::array_t<STATE>(
                      {sol.Rows(), sol.Cols()},
                      {sizeof(STATE), sizeof(STATE) * sol.Rows()},
                      sol.Elem(), self); }, "Zero-copy view of the solution vector/matrix as a numpy array.");

  m.def("MakePreconditionerFrom", &MakePreconditionerFrom, py::arg("an"), py::arg("precondDecompose"), "Builds a preconditioner from a snapshot of an already-assembled "
                                                                                                       "analysis's matrix. Call an.Assemble() first.");
}
