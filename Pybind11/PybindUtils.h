/// @file PybindUtils.h
/// @brief Helpers shared by the binding files

#ifndef PZ_PYBIND_UTILS_H
#define PZ_PYBIND_UTILS_H

#include <sstream>
#include <string>

/// @brief Adapts Print(std::ostream&) to a Python string; constness deduced
template <class T>
std::string PrintToString(T &obj)
{
    std::ostringstream out;
    obj.Print(out);
    return out.str();
}

#endif // PZ_PYBIND_UTILS_H
