cmake_minimum_required(VERSION 3.26)
project(CenterOfMass)

set(CMAKE_CXX_STANDARD 17)


add_executable("${PROJECT_NAME}" main.cpp
        src/MolecularDynamics/MolecularDynamics.h)
add_subdirectory(src)

target_link_libraries("${PROJECT_NAME}" PUBLIC MolecularDynamics)
target_include_directories("${PROJECT_NAME}" PUBLIC
        "${PROJECT_BINARY_DIR}"
        "${PROJECT_SOURCE_DIR}/src")

install(TARGETS "${PROJECT_NAME}" DESTINATION bin)
