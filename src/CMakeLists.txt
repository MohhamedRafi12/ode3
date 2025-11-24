# src/CMakeLists.txt
cmake_minimum_required(VERSION 3.20)

find_package(ROOT REQUIRED COMPONENTS Core Graf Hist Gpad Tree Gui)
find_package(GSL REQUIRED)

include_directories(${ROOT_INCLUDE_DIRS})
include_directories(${CMAKE_SOURCE_DIR}/odelib)

set(EXECUTABLES
    baseball1
    baseball2
    RKnTest
    vterm
    RKnDemo
)

foreach(exec ${EXECUTABLES})
    add_executable(${exec} ${exec}.cpp)
    target_link_libraries(${exec} PRIVATE ode ${ROOT_LIBRARIES})
    set_target_properties(${exec} PROPERTIES
        MACOSX_RPATH ON
        INSTALL_RPATH "@loader_path/../odelib"
    )
endforeach()

# projGSL needs GSL
add_executable(projGSL projGSL.cpp)
target_link_libraries(projGSL PRIVATE ode ${ROOT_LIBRARIES} GSL::gsl GSL::gslcblas)
set_target_properties(projGSL PROPERTIES
    MACOSX_RPATH ON
    INSTALL_RPATH "@loader_path/../odelib"
)
