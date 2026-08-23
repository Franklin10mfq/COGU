cmake_minimum_required(VERSION 3.10)
project(cogu_scvx C)

# Compilacion cruzada: pasar -DCMAKE_TOOLCHAIN_FILE=arm.cmake para ARM
# Ejemplo: cmake -B build -DCMAKE_TOOLCHAIN_FILE=arm-none-eabi.cmake

# Solver generado por CVXPYgen — expone cpg_include y cpg_src al scope padre
add_subdirectory(solver/c)

# Ejecutable principal
add_executable(cogu_scvx
    scvx_main.c
    dynamics.c
    ${cpg_src}
)

# Includes: raiz del output (cpg_compat.h, dynamics.h) + headers del solver
target_include_directories(cogu_scvx PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${cpg_include}
)

# Matematicas
if(NOT MSVC)
    target_link_libraries(cogu_scvx PRIVATE m)
endif()
