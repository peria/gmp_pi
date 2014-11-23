project(gmp_pi)

cmake_minimum_required(VERSION 2.8)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif()

enable_testing()

find_path(GFLAGS_PATH gflags/gflags.h)
find_library(GFLAGS_LIBRARY gflags)
find_path(GLOG_PATH glog/logging.h)
find_library(GLOG_LIBRARY glog)

if(GFLAGS_LIBRARY AND GLOG_LIBRARY)
  set(USE_GFLAGS_AND_GLOG 1)
  include_directories(${GFLAGS_PATH})
  include_directories(${GLOG_PATH})
endif()

# -------------------
# Set include directories, c++ options, etc.

function(pi_add_cxx_flags flags)
  if(${CMAKE_GENERATOR} MATCHES "Visual Studio")
  else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${flags}" PARENT_SCOPE)
  endif()
endfunction()

include_directories(.)
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

# ---------------
# Add subdirectories

function(pi_target_link_libraries target)
  target_link_libraries(${target} ${GFLAGS_LIBRARY})
  target_link_libraries(${target} ${GLOG_LIBRARY})
  target_link_libraries(${target} gmp)
  target_link_libraries(${target} ${ARGN})
endfunction()

add_subdirectory(base)
add_subdirectory(drm)
