cmake_minimum_required(VERSION 3.10.0)

project(vector-field-decomposition)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
### Configure output locations
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Print the build type
# set(CMAKE_BUILD_TYPE "Debug")
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug Release" FORCE)
endif()
message(STATUS "cmake build type: ${CMAKE_BUILD_TYPE}")

### Configure the compiler
# This is a basic, decent setup that should do something sane on most compilers

if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")

  # using Clang (linux or apple) or GCC
  message("Using clang/gcc compiler flags")
  # -fopenmp
  # -O3 really make a difference! 
  SET(BASE_CXX_FLAGS "-std=c++11 -Wall -Wextra -O3 -fopenmp")
  SET(DISABLED_WARNINGS " -Wno-unused-parameter -Wno-unused-variable -Wno-unused-function -Wno-deprecated-declarations -Wno-missing-braces -Wno-unused-private-field")
  SET(TRACE_INCLUDES " -H -Wno-error=unused-command-line-argument")

  if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    message("Setting clang-specific options")
    SET(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} -ferror-limit=3 -fcolor-diagnostics")
    SET(CMAKE_CXX_FLAGS_DEBUG          "-g3 -fsanitize=address -fno-limit-debug-info")
  elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    message("Setting gcc-specific options")
    SET(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} -fmax-errors=5")
    SET(CMAKE_CXX_FLAGS_DEBUG          "-g3")
    SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} -Wno-maybe-uninitialized -Wno-format-zero-length -Wno-unused-but-set-parameter -Wno-unused-but-set-variable")
  endif()
  
  SET(CMAKE_CXX_FLAGS "${BASE_CXX_FLAGS} ${DISABLED_WARNINGS}")
  
  # https://gitlab.kitware.com/cmake/cmake/-/issues/22484
  include(CheckCXXCompilerFlag)
  unset(COMPILER_SUPPORTS_MARCH_NATIVE CACHE)
  CHECK_CXX_COMPILER_FLAG(-march=native COMPILER_SUPPORTS_MARCH_NATIVE)
  if(COMPILER_SUPPORTS_MARCH_NATIVE)
    SET(CMAKE_CXX_FLAGS_RELEASE        "-O3 -march=native -DNDEBUG")
  else()
    SET(CMAKE_CXX_FLAGS_RELEASE        "-O3 -mcpu=apple-m1 -DNDEBUG") # Apple M1
  endif()

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  # using Visual Studio C++
  message("Using Visual Studio compiler flags")
  set(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} /W4")
  set(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} /MP") # parallel build
  SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} /wd\"4267\"")  # ignore conversion to smaller type (fires more aggressively than the gcc version, which is annoying)
  SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} /wd\"4244\"")  # ignore conversion to smaller type (fires more aggressively than the gcc version, which is annoying)
  SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} /wd\"4305\"")  # ignore truncation on initialization
  SET(CMAKE_CXX_FLAGS "${BASE_CXX_FLAGS} ${DISABLED_WARNINGS}")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MD")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MDd")
  set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

  add_definitions(/D "_CRT_SECURE_NO_WARNINGS")
  add_definitions(-DNOMINMAX)
  add_definitions(-D_USE_MATH_DEFINES)
else()
  # unrecognized
  message( FATAL_ERROR "Unrecognized compiler [${CMAKE_CXX_COMPILER_ID}]" )
endif()

# Edit files to build here.
set(SRCS 
  src/main.cpp
  src/include/Curve.cpp
#   src/tree-cotree.cpp
#   src/harmonic-bases.cpp
#   src/hodge-decomposition.cpp
#   ../poisson-problem/src/scalar-poisson-problem.cpp
#   ../../utils/src/solvers.cpp
)

# set(TEST_SRCS
#   src/tree-cotree.cpp
#   src/harmonic-bases.cpp
#   src/hodge-decomposition.cpp
#   ../../utils/src/solvers.cpp
#   src/test-decomp.cpp
# )
list(APPEND CMAKE_PREFIX_PATH "/opt/homebrew/share/eigen3/cmake")

# 寻找 Eigen3
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# 检查是否找到 Eigen
if(NOT Eigen3_FOUND)
    message(FATAL_ERROR "Eigen not found")
endif()
# add_subdirectory(${CMAKE_SOURCE_DIR}/deps/geometry-central)
add_subdirectory(${CMAKE_SOURCE_DIR}/deps/polyscope)
# find_library(Eigen NAMES eigen PATHS /opt/homebrew/lib /usr/local/lib)
# add_compile_definitions()

find_path(SUITESPARSE_INCLUDE_DIR SuiteSparse_config.h
          PATHS /opt/homebrew/include /usr/local/include)

find_library(SUITESPARSE_LIBRARY suitesparseconfig
             PATHS /opt/homebrew/lib /usr/local/lib)

find_path(SUITESPARSE_INCLUDE_DIR NAMES SuiteSparse_config.h PATHS /opt/homebrew/include/suitesparse /usr/local/include)
include_directories(${SUITESPARSE_INCLUDE_DIR})

find_library(CHOLMOD_LIBRARY NAMES cholmod PATHS /opt/homebrew/lib /usr/local/lib)
# find_package(SuiteSparse QUIET)

# find_package(SuiteSparse CONFIG QUIET)
# if(SuiteSparse_FOUND)
  add_definitions(-DUSE_SUITESPARSE)
  message(STATUS "Using SuiteSparse: ${SuiteSparse_VERSION}")
  include_directories(${SuiteSparse_INCLUDE_DIRS})
  link_libraries(${SuiteSparse_LIBRARIES})
# else()
  # message(STATUS "SuiteSparse not found, using Eigen for solvers")
# endif()

# add_subdirectory(../../deps/googletest deps/googletest)

find_library(CHOLMOD_LIBRARY NAMES cholmod PATHS /opt/homebrew/lib /usr/local/lib REQUIRED)
find_library(AMD_LIBRARY NAMES amd PATHS /opt/homebrew/lib /usr/local/lib REQUIRED)
find_library(COLAMD_LIBRARY NAMES colamd PATHS /opt/homebrew/lib /usr/local/lib REQUIRED)


# Change the executable name by changing "DER" in the lines below.
add_executable(DER "${SRCS}")
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT DER)
# target_include_directories(DER PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include/")
# target_include_directories(DER PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../utils/include")
# target_include_directories(DER PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../poisson-problem/include")
# add the args.hxx project which we use for command line args
# target_include_directories(DER PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../deps/polyscope/deps/args")
target_link_libraries(DER  polyscope Eigen3::Eigen
${CHOLMOD_LIBRARY}
${AMD_LIBRARY}
${COLAMD_LIBRARY})

# target_link_libraries(DER PUBLIC 

#   # 其他依赖库
# )
# if(SuiteSparse_FOUND)
# target_link_directories(DER PUBLIC ${SuiteSparse_INCLUDE_DIRS})
# target_link_libraries(DER ${SuiteSparse_LIBRARIES})
# endif()
# add_executable(test-decomp "${TEST_SRCS}")
# target_include_directories(test-decomp PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include/")
# target_include_directories(test-decomp PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../utils/include")
# target_link_libraries(test-decomp geometry-central polyscope gtest)
