pybind11_add_module(pyneopz module.cpp)

target_link_libraries(pyneopz PRIVATE pz)