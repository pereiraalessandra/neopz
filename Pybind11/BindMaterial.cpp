/// @file BindMaterial.cpp
/// @brief Bindings for the material interfaces and concrete materials (Material/)
/// @note materials inserted into a mesh are deleted by TPZCompMesh::CleanUp,
/// hence the py::nodelete holder. TPZBndCond is a sibling of TPZMaterial, not
/// derived from it, which is why InsertMaterialObject is overloaded on both

#include "Bindings.h"

#include "TPZMaterial.h"
#include "TPZBndCond.h"
#include "TPZBndCondT.h"
#include "pzfmatrix.h"
#include "DarcyFlow/TPZDarcyFlow.h"
#include "DarcyFlow/TPZMixedDarcyFlow.h"

namespace
{
    template <class TMaterial>
    TPZBndCond *CreateBoundaryCondition(TMaterial &material, int id, int type,
                                        const TPZFMatrix<STATE> &val1,
                                        const TPZVec<STATE> &val2)
    {
        return material.CreateBC(&material, id, type, val1, val2);
    }

    /// TPZIsotropicPermeability is a virtual base: no pointer to member
    template <class TMaterial>
    void SetConstantPermeability(TMaterial &material, STATE permeability)
    {
        material.SetConstantPermeability(permeability);
    }
} // namespace

void InitMaterial(py::module_ &m)
{
    py::class_<TPZMaterial>(m, "TPZMaterial", "Interface of every NeoPZ material.")
        .def("Id", &TPZMaterial::Id)
        .def("NStateVariables", &TPZMaterial::NStateVariables)
        .def("Dimension", &TPZMaterial::Dimension)
        .def("__str__", [](const TPZMaterial &material)
             { return PrintToString(material); });

    py::class_<TPZBndCond>(m, "TPZBndCond",
                           "Type agnostic interface of a boundary condition.")
        .def("Id", &TPZBndCond::Id)
        .def("Type", &TPZBndCond::Type, "0 is Dirichlet, 1 is Neumann.")
        .def("__str__", [](const TPZBndCond &bndcond)
             { return PrintToString(bndcond); });

    py::class_<TPZDarcyFlow, TPZMaterial,
               std::unique_ptr<TPZDarcyFlow, py::nodelete>>(
        m, "TPZDarcyFlow", "Darcy material for a single H1 space.")
        .def(py::init<int, int>(), py::arg("id"), py::arg("dim"))
        .def("SetConstantPermeability", &SetConstantPermeability<TPZDarcyFlow>,
             py::arg("permeability"))
        .def("CreateBC", &CreateBoundaryCondition<TPZDarcyFlow>,
             py::arg("id"), py::arg("type"), py::arg("val1"), py::arg("val2"),
             py::return_value_policy::reference, py::keep_alive<0, 1>());

    py::class_<TPZMixedDarcyFlow, TPZMaterial,
               std::unique_ptr<TPZMixedDarcyFlow, py::nodelete>>(
        m, "TPZMixedDarcyFlow",
        "Darcy material for combined spaces (HDiv flux, L2 pressure).")
        .def(py::init<int, int>(), py::arg("id"), py::arg("dim"))
        .def("SetConstantPermeability", &SetConstantPermeability<TPZMixedDarcyFlow>,
             py::arg("permeability"))
        .def("CreateBC", &CreateBoundaryCondition<TPZMixedDarcyFlow>,
             py::arg("id"), py::arg("type"), py::arg("val1"), py::arg("val2"),
             py::return_value_policy::reference, py::keep_alive<0, 1>());
}
