
cet_set_compiler_flags(
  DIAGS CAUTIOUS
  WERROR
  NO_UNDEFINED
  EXTRA_FLAGS -pedantic
  $<$<COMPILE_LANGUAGE:CXX>:-Wno-deprecated-enum-float-conversion>
  $<$<COMPILE_LANGUAGE:CXX>:-Wno-deprecated-enum-enum-conversion>
)
cet_report_compiler_flags(REPORT_THRESHOLD VERBOSE)

find_package(PandoraSDK 05.01.00 REQUIRED EXPORT)

option(PANDORA_MONITORING "Enable Pandora Monitoring" TRUE)
if(PANDORA_MONITORING)
  find_package(PandoraMonitoring 05.00.00 REQUIRED EXPORT)
endif()

find_package(Eigen3 3.3 REQUIRED EXPORT)

find_package(TBB QUIET EXPORT)
find_package(Torch QUIET EXPORT)
if(Torch_FOUND)
  set(PANDORA_LIBTORCH ON)
  message(STATUS "LibTorch found — building DL content (LArPandoraDLContent)")
else()
  set(PANDORA_LIBTORCH OFF)
  message(STATUS "LibTorch not found — skipping DL content build")
endif()

include(${PANDORA_PROJECT_ROOT}/cmake/LArContent_sources.cmake)

if(NOT DEFINED LAR_CONTENT_SRCS OR LAR_CONTENT_SRCS STREQUAL "")
  message(FATAL_ERROR "LAR_CONTENT_SRCS not defined or empty — check LArContent_sources.cmake")
endif()

cet_make_library(
  LIBRARY_NAME LArPandoraContent
  VERSION ${PROJECT_VERSION}
  SOVERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}
  SOURCE ${LAR_CONTENT_SRCS}
  LIBRARIES
    PUBLIC
      PandoraPFA::PandoraSDK
      PandoraPFA::PandoraMonitoring
    PRIVATE
      Eigen3::Eigen
)

if(PANDORA_MONITORING)
  target_compile_definitions(LArPandoraContent PUBLIC MONITORING)
endif()

if(PANDORA_LIBTORCH)
  include(${PANDORA_PROJECT_ROOT}/cmake/LArDLContent_sources.cmake)

  if(NOT DEFINED LAR_DL_CONTENT_SRCS OR LAR_DL_CONTENT_SRCS STREQUAL "")
    message(FATAL_ERROR "LAR_DL_CONTENT_SRCS not defined or empty — check LArDLContent_sources.cmake")
  endif()

  cet_make_library(
    LIBRARY_NAME LArPandoraDLContent
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}
    SOURCE ${LAR_DL_CONTENT_SRCS}
    LIBRARIES
      PUBLIC
        LArPandoraContent
        PandoraPFA::PandoraSDK
        torch
        Eigen3::Eigen
  )

  target_compile_definitions(LArPandoraDLContent PUBLIC PANDORA_LIBTORCH)
endif()

set(LAR_CONTENT_DIRS
  larpandoracontent
  larpandoracontent/LArCheating
  larpandoracontent/LArControlFlow
  larpandoracontent/LArCustomParticles
  larpandoracontent/LArHelpers
  larpandoracontent/LArMetrics
  larpandoracontent/LArMonitoring
  larpandoracontent/LArObjects
  larpandoracontent/LArPersistency
  larpandoracontent/LArPlugins
  larpandoracontent/LArReclustering
  larpandoracontent/LArReclustering/LArExample
  larpandoracontent/LArShowerRefinement
  larpandoracontent/LArThreeDReco
  larpandoracontent/LArThreeDReco/LArCosmicRay
  larpandoracontent/LArThreeDReco/LArEventBuilding
  larpandoracontent/LArThreeDReco/LArHitCreation
  larpandoracontent/LArThreeDReco/LArLongitudinalTrackMatching
  larpandoracontent/LArThreeDReco/LArPfoMopUp
  larpandoracontent/LArThreeDReco/LArPfoRecovery
  larpandoracontent/LArThreeDReco/LArPfoStitching
  larpandoracontent/LArThreeDReco/LArShowerFragments
  larpandoracontent/LArThreeDReco/LArShowerMatching
  larpandoracontent/LArThreeDReco/LArThreeDBase
  larpandoracontent/LArThreeDReco/LArTrackFragments
  larpandoracontent/LArThreeDReco/LArTransverseTrackMatching
  larpandoracontent/LArThreeDReco/LArTwoViewMatching
  larpandoracontent/LArTrackShowerId
  larpandoracontent/LArTwoDReco
  larpandoracontent/LArTwoDReco/LArClusterAssociation
  larpandoracontent/LArTwoDReco/LArClusterCreation
  larpandoracontent/LArTwoDReco/LArClusterMopUp
  larpandoracontent/LArTwoDReco/LArClusterSplitting
  larpandoracontent/LArTwoDReco/LArCosmicRay
  larpandoracontent/LArUtility
  larpandoracontent/LArVertex
)
install_source(SUBDIRS ${LAR_CONTENT_DIRS})
install_headers(SUBDIRS ${LAR_CONTENT_DIRS})

if(PANDORA_LIBTORCH)
  set(LAR_DL_DIRS
    larpandoradlcontent
    larpandoradlcontent/LArCheating
    larpandoradlcontent/LArControlFlow
    larpandoradlcontent/LArEventClassification
    larpandoradlcontent/LArHelpers
    larpandoradlcontent/LArMonitoring
    larpandoradlcontent/LArObjects
    larpandoradlcontent/LArShowerGrowing
    larpandoradlcontent/LArSignalId
    larpandoradlcontent/LArThreeDReco
    larpandoradlcontent/LArThreeDReco/LArEventBuilding
    larpandoradlcontent/LArTrackShowerId
    larpandoradlcontent/LArTwoDReco
    larpandoradlcontent/LArVertex
  )
  install_source(SUBDIRS ${LAR_DL_DIRS})
  install_headers(SUBDIRS ${LAR_DL_DIRS})
endif()

cet_cmake_config()
