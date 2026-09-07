cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

find_package(cetmodules 3.13.00 QUIET)

# Common version must be applied manually pending ups compatible solution
# set(LARCONTENT_VERSION 05.03.00)

if(cetmodules_FOUND)
  project(larpandoracontent VERSION 05.04.00 LANGUAGES CXX)
else()
  project(LArContent VERSION 05.05.00 LANGUAGES CXX)
endif()

if(cetmodules_FOUND)
    # Deal with MRB changing top-level CMake invocation location
    set(PANDORA_PROJECT_ROOT "${CMAKE_CURRENT_SOURCE_DIR}" CACHE PATH "Top-level Pandora project root")
    message(STATUS "Configuring with cetmodules build system ${PANDORA_PROJECT_ROOT}.")
    include(CetCMakeEnv)
    cet_cmake_env()
    include(cmake/cetmodules_build.cmake)
else()
    message(STATUS "Configuring with standalone build system.")
    include(cmake/standalone_build.cmake)
endif()
