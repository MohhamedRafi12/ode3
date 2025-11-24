cmake_minimum_required(VERSION 3.20)
project(ODE3 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Enable RPATH on macOS + UNIX
set(CMAKE_MACOSX_RPATH ON)

# Add subdirectories
add_subdirectory(odelib)
add_subdirectory(src)
