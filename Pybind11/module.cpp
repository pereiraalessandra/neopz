#include<pybind11/pybind11.h>

PYBIND11_MODULE(pyneopz, m) {
    m.doc() = "pybind11 example plugin"; // optional module docstring}
}