cmake_minimum_required(VERSION 3.16.0)
project(ReMatching LANGUAGES C CXX)


# Build and install libCUT
include(ExternalProject)
ExternalProject_Add(CUT
                    SOURCE_DIR "${CMAKE_SOURCE_DIR}/ext/cut"
                    BINARY_DIR "${CMAKE_SOURCE_DIR}/build-cut"
                    INSTALL_DIR "${CMAKE_SOURCE_DIR}/install"
                    CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX=${CMAKE_SOURCE_DIR}/install")


# Include dependencies
include_directories("${CMAKE_SOURCE_DIR}/ext/eigen")
include_directories("${CMAKE_SOURCE_DIR}/ext/libigl/include")
include_directories("${CMAKE_SOURCE_DIR}/ext/json/include")
include_directories("${CMAKE_SOURCE_DIR}/install/cut/include")
if(MSVC)
    set(libCUT "${CMAKE_SOURCE_DIR}/install/cut/lib/cut.lib")
elseif(UNIX)
    set(libCUT "${CMAKE_SOURCE_DIR}/install/cut/lib/libcut.a")
endif()
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)



# Include rematching library
include_directories("${CMAKE_SOURCE_DIR}/include")
add_library(RMT STATIC  "${CMAKE_SOURCE_DIR}/src/rmt/graph.cpp"
                        "${CMAKE_SOURCE_DIR}/src/rmt/mesh.cpp"
                        "${CMAKE_SOURCE_DIR}/src/rmt/voronoifps.cpp"
                        "${CMAKE_SOURCE_DIR}/src/rmt/region.cpp"
                        "${CMAKE_SOURCE_DIR}/src/rmt/regiondict.cpp"
                        "${CMAKE_SOURCE_DIR}/src/rmt/flatunion.cpp"
                        "${CMAKE_SOURCE_DIR}/src/rmt/reconstruction.cpp"
                        "${CMAKE_SOURCE_DIR}/src/rmt/preprocess.cpp"
                        "${CMAKE_SOURCE_DIR}/src/rmt/clean.cpp"
                        "${CMAKE_SOURCE_DIR}/src/rmt/weightmap.cpp"
                        "${CMAKE_SOURCE_DIR}/src/rmt/eval.cpp"
                        "${CMAKE_SOURCE_DIR}/src/rmt/io.cpp"
                        "${CMAKE_SOURCE_DIR}/src/rmt/rmt.cpp")
target_link_libraries(RMT ${libCUT} Threads::Threads)
set_target_properties(RMT PROPERTIES CXX_STANDARD 17)


set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/pyrmt")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/pyrmt")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/pyrmt")
foreach( OUTPUTCONFIG ${CMAKE_CONFIGURATION_TYPES} )
    string( TOUPPER ${OUTPUTCONFIG} OUTPUTCONFIG )
    set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY_${OUTPUTCONFIG} "${CMAKE_BINARY_DIR}/pyrmt")
    set( CMAKE_LIBRARY_OUTPUT_DIRECTORY_${OUTPUTCONFIG} "${CMAKE_BINARY_DIR}/pyrmt")
    set( CMAKE_RUNTIME_OUTPUT_DIRECTORY_${OUTPUTCONFIG} "${CMAKE_BINARY_DIR}/pyrmt")
endforeach( OUTPUTCONFIG CMAKE_CONFIGURATION_TYPES )


# Find PyBind
set(pybind11_DIR "path/to/pybind11/share/cmake/pybind11")
find_package(pybind11 REQUIRED)

# Add pybind module
pybind11_add_module(PyRMT "${CMAKE_SOURCE_DIR}/src/python/pyrmt.cpp")
target_link_libraries(PyRMT PRIVATE RMT)

# Copy example
configure_file("${CMAKE_SOURCE_DIR}/python/sample.py" "${CMAKE_BINARY_DIR}/pyrmt/sample.py" COPYONLY)