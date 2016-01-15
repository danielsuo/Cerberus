cmake_minimum_required (VERSION 2.6.0)
project (Cerberus CXX)

find_package (Ceres REQUIRED)

set (SOURCES
    main.cpp
    Cerberus.cpp
  )

add_library (Cerberus ${SOURCES})
target_link_libraries (Cerberus ${CERES_LIBRARIES})
target_include_directories (Cerberus PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories (Cerberus PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/Residuals)
