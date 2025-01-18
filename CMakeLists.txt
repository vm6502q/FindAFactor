cmake_minimum_required(VERSION 3.15...3.29)
project(FindAFactor LANGUAGES CXX)

set(PYBIND11_FINDPYTHON ON)
find_package(pybind11 CONFIG REQUIRED)
find_package(Boost CONFIG REQUIRED)

set(_USE_MT_RNG OFF)
if (CMAKE_SYSTEM_NAME MATCHES "Linux")
    if (CMAKE_SYSTEM_VERSION MATCHES "22\\.04")
        set(_USE_MT_RNG ON)
    endif()
endif()

option (USE_MT_RNG "Use mt19937 random number generator instead of taus88" _USE_MT_RNG)

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