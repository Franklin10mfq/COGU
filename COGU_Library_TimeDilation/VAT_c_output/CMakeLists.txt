cmake_minimum_required(VERSION 3.10)
project(cogu_scvx C)

# Cross compiling: pass -DCMAKE_TOOLCHAIN_FILE=<your-toolchain>.cmake
# e.g. cmake -B build -DCMAKE_TOOLCHAIN_FILE=arm-none-eabi.cmake

# CVXPYgen solver; exports cpg_include and cpg_src to the parent scope
add_subdirectory(solver/c)

add_library(cogu_scvx SHARED
    scvx_main.c
    dynamics.c
    ${cpg_src}
)
# COGU_API in cogu_scvx.h expands to __declspec(dllexport) only when this is
# defined, so a *consumer* of the .dll (Python via ctypes, a Simulink S-
# function) does not need to define COGU_USE_DLL to link against it.
target_compile_definitions(cogu_scvx PRIVATE COGU_BUILD_DLL)

target_include_directories(cogu_scvx PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${cpg_include}
)

if(NOT MSVC)
    target_link_libraries(cogu_scvx PRIVATE m)
endif()


# Static library, for Simulink's C Caller block (see COGU_SCVx_C_caller.h):
# a static .lib/.a is what the C Caller block links against, no separate
# .dll/.so needs to sit next to the built model at runtime. Not used by the
# Python callers -- those load the SHARED library above via ctypes, which
# needs a .dll/.so it can dlopen() at runtime. COGU_API expands to nothing
# here (COGU_BUILD_DLL is not defined for this target), which is correct: a
# static archive has no import/export boundary to decorate.
add_library(cogu_scvx_static STATIC
    scvx_main.c
    dynamics.c
    ${cpg_src}
)
set_target_properties(cogu_scvx_static PROPERTIES OUTPUT_NAME cogu_scvx_static)
target_include_directories(cogu_scvx_static PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${cpg_include}
)
if(NOT MSVC)
    target_link_libraries(cogu_scvx_static PRIVATE m)
endif()
