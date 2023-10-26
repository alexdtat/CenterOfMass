add_library(MolecularDynamics
        MolecularDynamics.h
        MolecularDynamics.cpp
)

target_link_libraries("${PROJECT_NAME}" PUBLIC MolecularDynamics)
target_include_directories("${PROJECT_NAME}" PUBLIC
        "${PROJECT_BINARY_DIR}"
        "${PROJECT_SOURCE_DIR}/src")