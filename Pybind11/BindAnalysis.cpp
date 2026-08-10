/// @file BindAnalysis.cpp
/// @brief Bindings for the assembly and solution of the linear system
/// (Analysis/, StrMatrix/, Solvers/)

#include "Bindings.h"

#include "pzcmesh.h"
#include "TPZLinearAnalysis.h"
#include "TPZStructMatrix.h"
#include "pzskylstrmatrix.h"
#include "TPZSolver.h"
#include "pzstepsolver.h"

void InitAnalysis(py::module_ &m)
{
     /// bare registrations: they only need to exist for the setters below
     py::class_<TPZStructMatrix>(m, "TPZStructMatrix", "Base of the structural matrices.");
     py::class_<TPZSolver>(m, "TPZSolver", "Base of the solvers.");

     py::class_<TPZSkylineStructMatrix<STATE>, TPZStructMatrix>(
         m, "TPZSkylineStructMatrix", "Structural matrix in skyline format.")
         .def(py::init<TPZCompMesh *>(), py::arg("cmesh"));

     py::class_<TPZStepSolver<STATE>, TPZSolver>(m, "TPZStepSolver", "Solver driver.")
         .def(py::init<>())
         .def("SetDirect", &TPZStepSolver<STATE>::SetDirect, py::arg("decomposeType"),
              "Mixed formulations are symmetric indefinite: use ELDLt.");

     py::class_<TPZLinearAnalysis>(m, "TPZLinearAnalysis",
                                   "Assembles and solves a linear problem.")
         .def(py::init<TPZCompMesh *>(), py::arg("cmesh"), py::keep_alive<1, 2>())
         /// no keep_alive: the argument is cloned
         .def("SetStructuralMatrix",
              py::overload_cast<TPZStructMatrix &>(&TPZLinearAnalysis::SetStructuralMatrix),
              py::arg("strmatrix"))
         .def("SetSolver", &TPZLinearAnalysis::SetSolver, py::arg("solver"))
         .def("Assemble", &TPZLinearAnalysis::Assemble)
         .def("Solve", &TPZLinearAnalysis::Solve,
              "Loads the solution into the mesh.")
         .def("SolutionNorm", [](TPZLinearAnalysis &self)
              { return Norm(self.Solution()); });
}
