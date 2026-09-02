cmake_minimum_required(VERSION 3.16)
project(resonator_fem CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(FetchContent)
FetchContent_Declare(
  Eigen
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
  GIT_TAG 3.4.0
)
set(EIGEN_BUILD_DOC OFF CACHE BOOL "" FORCE)
set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(Eigen)

add_library(resonator_core src/beam.cpp src/environment.cpp)
target_include_directories(resonator_core PUBLIC include)
target_link_libraries(resonator_core PUBLIC Eigen3::Eigen)

foreach(ex modal_convergence prestress_sweep coupled_history)
  add_executable(${ex} examples/${ex}.cpp)
  target_link_libraries(${ex} PRIVATE resonator_core)
endforeach()

enable_testing()
foreach(t test_beam test_environment)
  add_executable(${t} tests/${t}.cpp)
  target_link_libraries(${t} PRIVATE resonator_core)
  add_test(NAME ${t} COMMAND ${t})
endforeach()
