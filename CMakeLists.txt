cmake_minimum_required(VERSION 3.15...3.29)
project(FindAFactor LANGUAGES CXX)

set(PYBIND11_FINDPYTHON ON)
find_package(pybind11 CONFIG REQUIRED)
find_package(Boost CONFIG REQUIRED)

# Check if the system is Ubuntu
if (CMAKE_SYSTEM_NAME MATCHES "Linux")
    execute_process(
        COMMAND lsb_release -a
        OUTPUT_VARIABLE LSB_RELEASE_OUTPUT
    )

    # Extract the version number
    string(REGEX MATCH "Release:[ \t]+([0-9.]+)" 
           UBUNTU_VERSION "${LSB_RELEASE_OUTPUT}")
    set(UBUNTU_VERSION "${CMAKE_MATCH_1}")

    message(STATUS "Ubuntu version: ${UBUNTU_VERSION}")

    if ("${UBUNTU_VERSION}" STREQUAL "22.04")
        option (USE_MT_RNG "Use mt19937 random number generator instead of taus88" ON)
    else ("${UBUNTU_VERSION}" STREQUAL "22.04")
        option (USE_MT_RNG "Use mt19937 random number generator instead of taus88" OFF)
    endif ("${UBUNTU_VERSION}" STREQUAL "22.04")
else (CMAKE_SYSTEM_NAME MATCHES "Linux")
    option (USE_MT_RNG "Use mt19937 random number generator instead of taus88" OFF)
endif (CMAKE_SYSTEM_NAME MATCHES "Linux")

message(STATUS "USE_MT_RNG: ${USE_MT_RNG}")

configure_file(FindAFactor/include/common/config.h.in ${CMAKE_CURRENT_BINARY_DIR}/include/common/config.h @ONLY)

pybind11_add_module(_find_a_factor FindAFactor/_find_a_factor.cpp FindAFactor/dispatchqueue.cpp)

if (DEFINED ENV{BOOST_ROOT})
  target_include_directories(_find_a_factor PUBLIC FindAFactor/include ${CMAKE_CURRENT_BINARY_DIR}/include $ENV{BOOST_ROOT}/include/boost)
else (DEFINED ENV{BOOST_ROOT})
  target_include_directories(_find_a_factor PUBLIC FindAFactor/include ${CMAKE_CURRENT_BINARY_DIR}/include)
endif (DEFINED ENV{BOOST_ROOT})

if (MSVC)
  target_compile_options(_find_a_factor PUBLIC /O2 /std:c++14)
else (MSVC)
  target_compile_options(_find_a_factor PUBLIC -O3 -std=c++14 -lpthread)
endif (MSVC)

install(TARGETS _find_a_factor DESTINATION .)