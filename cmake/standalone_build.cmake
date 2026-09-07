# Enforce out-of-source builds.
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    message(FATAL_ERROR "${CMAKE_PROJECT_NAME} requires an out-of-source build.")
endif()

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/install" CACHE PATH "Default install path" FORCE)
endif()

# Build Options
option(PANDORA_LIBTORCH "Build the LArDLContent library against LibTorch" OFF)
option(LArContent_BUILD_DOCS "Build documentation for ${PROJECT_NAME}" OFF)

# Find Dependencies
if (NOT TARGET PandoraPFA::PandoraSDK)
    find_package(PandoraSDK 05.01.00 REQUIRED)
endif()
option(PANDORA_MONITORING "Enable Pandora Monitoring" ON)
if(PANDORA_MONITORING)
    if (NOT TARGET PandoraPFA::PandoraMonitoring)
        find_package(PandoraMonitoring 05.00.00 REQUIRED)
    endif()
endif()
find_package(Eigen3 3.3 REQUIRED)
if(PANDORA_LIBTORCH)
    find_package(TBB QUIET)
    find_package(Torch REQUIRED)
endif()

#include(PandoraCMakeSettings)

# --- Target: LArContent ---
include(cmake/LArContent_sources.cmake)
add_library(${PROJECT_NAME} SHARED ${LAR_CONTENT_SRCS})
add_library(PandoraPFA::${PROJECT_NAME} ALIAS ${PROJECT_NAME})

set_target_properties(${PROJECT_NAME} PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}
)

target_include_directories(LArContent PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
)

target_link_libraries(${PROJECT_NAME} PUBLIC
    PandoraPFA::PandoraSDK
    Eigen3::Eigen
)

if(PandoraMonitoring_FOUND)
    target_link_libraries(${PROJECT_NAME} PUBLIC PandoraPFA::PandoraMonitoring)
    target_compile_definitions(${PROJECT_NAME} PUBLIC MONITORING)
endif()

set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED ON)

#-------------------------------------------------------------------------------------------------------------------------------------------
# Target: LArDLContent (Optional)
if(PANDORA_LIBTORCH)
    set(DL_PROJECT_NAME "LArDLContent")

    include(cmake/LArDLContent_sources.cmake)
    add_library(${DL_PROJECT_NAME} SHARED ${LAR_DL_CONTENT_SRCS})
    add_library(PandoraPFA::${DL_PROJECT_NAME} ALIAS ${DL_PROJECT_NAME})

    set_target_properties(${DL_PROJECT_NAME} PROPERTIES
        VERSION ${PROJECT_VERSION}
        SOVERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}
    )

    target_include_directories(${DL_PROJECT_NAME} PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:include>
    )

######## Modern cmake version not currently support by LibTorch ###
#    target_link_libraries(${DL_PROJECT_NAME} PUBLIC
#        PandoraPFA::${PROJECT_NAME}
#        Torch::torch
#    )
####### End modern version

######## Backward compatible version
    target_link_libraries(${DL_PROJECT_NAME} PUBLIC PandoraPFA::LArContent)
    target_include_directories(${DL_PROJECT_NAME} PUBLIC ${TORCH_INCLUDE_DIRS})
    target_link_libraries(${DL_PROJECT_NAME} PUBLIC ${TORCH_LIBRARIES})
    target_compile_options(${DL_PROJECT_NAME} PRIVATE ${TORCH_CXX_FLAGS})
######## End backward compatible version

    set_target_properties(${DL_PROJECT_NAME} PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED ON)
endif()

# --- Installation ---
include(CMakePackageConfigHelpers)
include(GNUInstallDirs)

set(LIBRARIES_TO_INSTALL ${PROJECT_NAME})
if(PANDORA_LIBTORCH)
    list(APPEND LIBRARIES_TO_INSTALL ${DL_PROJECT_NAME})
endif()

foreach(LIB_NAME IN LISTS LIBRARIES_TO_INSTALL)
    install(TARGETS ${LIB_NAME}
        EXPORT ${LIB_NAME}Targets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
    )

    string(TOLOWER ${LIB_NAME} SOURCE_DIR_NAME) # LArContent -> larcontent
    if (LIB_NAME STREQUAL "LArContent")
        set(SOURCE_DIR_NAME "larpandoracontent")
    elseif(LIB_NAME STREQUAL "LArDLContent")
        set(SOURCE_DIR_NAME "larpandoradlcontent")
    endif()
    install(DIRECTORY ${SOURCE_DIR_NAME}/
        DESTINATION include/${SOURCE_DIR_NAME}
        COMPONENT Development
        FILES_MATCHING PATTERN "*.h"
    )

    write_basic_package_version_file(
        "${CMAKE_CURRENT_BINARY_DIR}/${LIB_NAME}ConfigVersion.cmake"
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY AnyNewerVersion
    )

    configure_package_config_file(
        "cmake/${LIB_NAME}Config.cmake.in"
        "${CMAKE_CURRENT_BINARY_DIR}/${LIB_NAME}Config.cmake"
        INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${LIB_NAME}
    )

    install(FILES
        "${CMAKE_CURRENT_BINARY_DIR}/${LIB_NAME}Config.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/${LIB_NAME}ConfigVersion.cmake"
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${LIB_NAME}
        COMPONENT Development
    )

    install(EXPORT ${LIB_NAME}Targets
        FILE ${LIB_NAME}Targets.cmake
        NAMESPACE PandoraPFA::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${LIB_NAME}
        COMPONENT Development
    )
endforeach()

if(LArContent_BUILD_DOCS)
    add_subdirectory(doc)
endif()

