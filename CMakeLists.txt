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

add_executable(resonator src/main.cpp src/beam.cpp src/environment.cpp)
target_include_directories(resonator PRIVATE include)
target_link_libraries(resonator PRIVATE Eigen3::Eigen)
