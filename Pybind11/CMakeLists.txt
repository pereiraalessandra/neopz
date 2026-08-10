pybind11_add_module(pyneopz
    Module.cpp
    BindEnums.cpp
    BindGeom.cpp
    BindMaterial.cpp
    BindMesh.cpp
    BindApproxCreator.cpp
    BindAnalysis.cpp
    BindPost.cpp
    )

target_link_libraries(pyneopz PRIVATE pz)

set_target_properties(pyneopz PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/python)
