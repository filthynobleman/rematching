cmake_minimum_required(VERSION 3.16.0)
project(ReMatching LANGUAGES C CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)  # Add this line before any libraries are defined

# Build and install libCUT
message("Patching cmake file of cut library...") # this is necessary so that we can build it with the project
set(CUT_CMAKES_LIST ${CMAKE_SOURCE_DIR}/ext/cut/CMakeLists.txt)
file(READ ${CUT_CMAKES_LIST} FILE_CONTENTS)
string(REPLACE "CMAKE_SOURCE_DIR" "PROJECT_SOURCE_DIR" FILE_CONTENTS "${FILE_CONTENTS}")
file(WRITE ${CUT_CMAKES_LIST} "${FILE_CONTENTS}")
message("Done.")
add_subdirectory(${CMAKE_SOURCE_DIR}/ext/cut)
include_directories("${CMAKE_SOURCE_DIR}/ext/cut/include/")


# Include dependencies
include_directories("${CMAKE_SOURCE_DIR}/ext/eigen")
include_directories("${CMAKE_SOURCE_DIR}/ext/libigl/include")
include_directories("${CMAKE_SOURCE_DIR}/ext/json/include")
include_directories("${CMAKE_SOURCE_DIR}/install/cut/include")
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
target_link_libraries(RMT cut Threads::Threads)
set_target_properties(RMT PROPERTIES 
                      CXX_STANDARD 17
                      POSITION_INDEPENDENT_CODE ON)

# we don't want to reset all those files => otherwise our setup script doesn't find the python package and it will consequently not be copied correctly to our site-packages directory
#set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/pyrmt")
#set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/pyrmt")
#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/pyrmt")
#foreach( OUTPUTCONFIG ${CMAKE_CONFIGURATION_TYPES} )
#    string( TOUPPER ${OUTPUTCONFIG} OUTPUTCONFIG )
#    set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY_${OUTPUTCONFIG} "${CMAKE_BINARY_DIR}/pyrmt")
#    set( CMAKE_LIBRARY_OUTPUT_DIRECTORY_${OUTPUTCONFIG} "${CMAKE_BINARY_DIR}/pyrmt")
#    set( CMAKE_RUNTIME_OUTPUT_DIRECTORY_${OUTPUTCONFIG} "${CMAKE_BINARY_DIR}/pyrmt")
#endforeach( OUTPUTCONFIG CMAKE_CONFIGURATION_TYPES )


# Include PyBind
add_subdirectory(ext/pybind11)

# Add pybind module
pybind11_add_module(PyRMT "${CMAKE_SOURCE_DIR}/src/python/pyrmt.cpp")
target_link_libraries(PyRMT PRIVATE RMT)

# Copy example
#configure_file("${CMAKE_SOURCE_DIR}/python/sample.py" "${CMAKE_BINARY_DIR}/pyrmt/sample.py" COPYONLY)