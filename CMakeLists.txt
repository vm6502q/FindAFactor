cmake_minimum_required(VERSION 3.15...3.29)
project(FindAFactor LANGUAGES CXX)

set(PYBIND11_FINDPYTHON ON)
find_package(pybind11 CONFIG REQUIRED)
find_package(Boost CONFIG REQUIRED)

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