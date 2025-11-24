# odelib/CMakeLists.txt
cmake_minimum_required(VERSION 3.20)

find_package(ROOT REQUIRED COMPONENTS Core Graf Hist Gpad Tree Gui)

include_directories(${ROOT_INCLUDE_DIRS})

add_library(ode SHARED RKn.cpp)
add_library(ode_static STATIC RKn.cpp)

target_link_libraries(ode PRIVATE ${ROOT_LIBRARIES})
target_link_libraries(ode_static PRIVATE ${ROOT_LIBRARIES})

# Proper macOS install_name
set_target_properties(ode PROPERTIES
    MACOSX_RPATH ON
    INSTALL_RPATH "@loader_path"
)

# Make headers available to src/
target_include_directories(ode PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(ode_static PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
