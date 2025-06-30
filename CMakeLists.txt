# Minimum CMake version & project
cmake_minimum_required(VERSION 3.21)
project(AVBDViewer LANGUAGES CXX)

# --------------------------------------------------------------
# Fetch SDL2 automatically (header‑only helper + static library)
# --------------------------------------------------------------
include(FetchContent)

# Tell CMake we want SDL2, disable tests/examples so the build stays tiny
set(SDL2_DISABLE_INSTALL ON CACHE BOOL "" FORCE)
set(SDL2_DISABLE_TESTS   ON CACHE BOOL "" FORCE)
set(SDL2_DISABLE_SAMPLES ON CACHE BOOL "" FORCE)

FetchContent_Declare(
  SDL2
  GIT_REPOSITORY https://github.com/libsdl-org/SDL.git
  GIT_TAG        release-2.32.0  # ↳ pick a specific tag for CI repeatability
)
FetchContent_MakeAvailable(SDL2)

# --------------------------------------------------------------
# Project sources
# --------------------------------------------------------------
add_executable(avbd_viewer
    AvbdToy.cpp
    physics/World.cpp
    # ... add other *.cpp here ...
)

# C++17, warnings etc.
set_target_properties(avbd_viewer PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED YES
)

# Link against SDL2::SDL2 and platform libs that SDL needs
# FetchContent gives us the alias target automatically.
if (WIN32)
    target_link_libraries(avbd_viewer PRIVATE SDL2::SDL2 SDL2::SDL2main)
else()
    target_link_libraries(avbd_viewer PRIVATE SDL2::SDL2)
    # pthread / dl are usually transitively linked, but add if needed:
    # target_link_libraries(avbd_viewer PRIVATE pthread dl)
endif()

# On Windows copy SDL2.dll next to the executable after build.
if (WIN32)
    add_custom_command(TARGET avbd_viewer POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<TARGET_FILE:SDL2::SDL2>
            $<TARGET_FILE_DIR:avbd_viewer>)
endif()