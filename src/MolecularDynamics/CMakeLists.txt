add_library(MolecularDynamics
        MolecularDynamics.h
        MolecularDynamics.cpp
)

# include_directories("../matplotplusplus/source/")
target_link_libraries("${PROJECT_NAME}" PUBLIC MolecularDynamics)
# target_link_libraries("${PROJECT_NAME}" PUBLIC matplot)
target_include_directories("${PROJECT_NAME}" PUBLIC
        "${PROJECT_BINARY_DIR}"
        "${PROJECT_SOURCE_DIR}/src")