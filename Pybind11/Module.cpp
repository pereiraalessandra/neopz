/// @file Module.cpp
/// @brief Entry point of the pyneopz extension module
/// @note the only file carrying PYBIND11_MODULE; the others expose init functions

#include "Bindings.h"

PYBIND11_MODULE(pyneopz, m)
{
    m.doc() = "Python bindings for NeoPZ";

    /// dependency order: bases and enums before whatever mentions them
    InitEnums(m);
    InitGeom(m);
    InitMaterial(m);
    InitMesh(m);
    InitApproxCreator(m);
    InitAnalysis(m);
    InitPost(m);
}