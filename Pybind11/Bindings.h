/// @file Bindings.h
/// @brief Shared includes and the declaration of every init function

#ifndef PZ_PYBIND_BINDINGS_H
#define PZ_PYBIND_BINDINGS_H

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "PybindCasters.h"
#include "PybindUtils.h"

namespace py = pybind11;

/// one init function per binding file, called from Module.cpp in dependency order
void InitEnums(py::module_ &m);
void InitGeom(py::module_ &m);
void InitMaterial(py::module_ &m);
void InitMesh(py::module_ &m);
void InitApproxCreator(py::module_ &m);
void InitAnalysis(py::module_ &m);
void InitPost(py::module_ &m);

#endif // PZ_PYBIND_BINDINGS_H