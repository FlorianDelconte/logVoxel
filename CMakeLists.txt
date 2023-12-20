
# -----------------------------------------------------------------------------
# Project definition
# -----------------------------------------------------------------------------
PROJECT(TLDDC)
cmake_minimum_required (VERSION 3.1)
cmake_policy(SET CMP0057 NEW)
# -----------------------------------------------------------------------------
# CPP11
# -----------------------------------------------------------------------------
set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_STANDARD_REQUIRED TRUE)
# -----------------------------------------------------------------------------
# DGtal dependencies
# -----------------------------------------------------------------------------
FIND_PACKAGE(DGtal 1.2 REQUIRED)
INCLUDE_DIRECTORIES(${DGTAL_INCLUDE_DIRS})
LINK_DIRECTORIES(${DGTAL_LIBRARY_DIRS})
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
message(STATUS "DGtal found.")
# -----------------------------------------------------------------------------
# Thread dependencies
# -----------------------------------------------------------------------------
FIND_PACKAGE(Threads)
# -----------------------------------------------------------------------------
# CLI11
# -----------------------------------------------------------------------------
include_directories( "${PROJECT_SOURCE_DIR}/ext/" )

ADD_EXECUTABLE(voxelizer Main_voxelizer Morphology3D)
TARGET_LINK_LIBRARIES(voxelizer  ${DGTAL_LIBRARIES})
