cmake_minimum_required(VERSION 3.16 FATAL_ERROR)

project(Collimator_Temperature_Rise
        VERSION 1.0
        LANGUAGES CXX)

# Geant4 is required.
option(WITH_GEANT4_UIVIS "Build with Geant4 UI and visualization drivers" ON)

if(WITH_GEANT4_UIVIS)
    find_package(Geant4 REQUIRED ui_all vis_all)
else()
    find_package(Geant4 REQUIRED)
endif()

# ROOT is used for the analysis output and analysis macros.
find_package(ROOT REQUIRED)

include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)

file(GLOB SOURCES CONFIGURE_DEPENDS ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB HEADERS CONFIGURE_DEPENDS ${PROJECT_SOURCE_DIR}/include/*.hh)

add_executable(collimator_temp_rise
    main.cc
    ${SOURCES}
    ${HEADERS}
)

target_compile_features(collimator_temp_rise PRIVATE cxx_std_17)
target_link_libraries(collimator_temp_rise
    PRIVATE
    ${Geant4_LIBRARIES}
    ROOT::Core
    ROOT::RIO
    ROOT::Tree
)

# Make the Geant4 command files available in the build directory.
configure_file(${PROJECT_SOURCE_DIR}/run.mac
               ${PROJECT_BINARY_DIR}/run.mac COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/vis.mac
               ${PROJECT_BINARY_DIR}/vis.mac COPYONLY)

install(TARGETS collimator_temp_rise DESTINATION bin)
