cmake_minimum_required (VERSION 2.6.0)
project (Cerberus CXX)

set(CMAKE_MODULE_PATH /usr/local/share/cmake/Modules)

find_package(Eigen3)
find_package (Ceres REQUIRED)

set (SOURCES
    Cerberus.cpp
  )

add_library (Cerberus ${SOURCES})
target_link_libraries (Cerberus ${CERES_LIBRARIES})
target_include_directories (Cerberus PUBLIC
        ${EIGEN3_INCLUDE_DIR}
        ${CERES_INCLUDE_DIRS}
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/Residuals)

add_executable (CerberusCLI EXCLUDE_FROM_ALL main.cpp ${SOURCES})
target_link_libraries (CerberusCLI ${CERES_LIBRARIES})
target_include_directories (CerberusCLI PUBLIC
        ${EIGEN3_INCLUDE_DIR}
        ${CERES_INCLUDE_DIRS}
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/Residuals)