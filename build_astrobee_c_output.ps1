# Compilar codigo generado por pipeline COGU en astrobee_c_output/
# Correr desde: c:\Users\luiss\COGU\
#
# Rutas correctas segun CVXPYgen: solver/c/src/, solver/c/include/, etc.

$env:PATH = "C:\msys64\ucrt64\bin;C:\msys64\usr\bin;" + $env:PATH

$OUT    = "astrobee_c_output"
$SOLVER = "$OUT\solver\c"

gcc -O0 -std=c99 -D_USE_MATH_DEFINES `
    "$OUT\scvx_main.c" `
    "$OUT\dynamics.c" `
    "$SOLVER\src\cpg_solve.c" `
    "$SOLVER\src\cpg_workspace.c" `
    "$SOLVER\solver_code\src\ecos.c" `
    "$SOLVER\solver_code\src\cone.c" `
    "$SOLVER\solver_code\src\ctrlc.c" `
    "$SOLVER\solver_code\src\equil.c" `
    "$SOLVER\solver_code\src\expcone.c" `
    "$SOLVER\solver_code\src\kkt.c" `
    "$SOLVER\solver_code\src\preproc.c" `
    "$SOLVER\solver_code\src\spla.c" `
    "$SOLVER\solver_code\src\splamm.c" `
    "$SOLVER\solver_code\src\timer.c" `
    "$SOLVER\solver_code\src\wright_omega.c" `
    "$SOLVER\solver_code\external\amd\src\amd_1.c" `
    "$SOLVER\solver_code\external\amd\src\amd_2.c" `
    "$SOLVER\solver_code\external\amd\src\amd_aat.c" `
    "$SOLVER\solver_code\external\amd\src\amd_control.c" `
    "$SOLVER\solver_code\external\amd\src\amd_defaults.c" `
    "$SOLVER\solver_code\external\amd\src\amd_dump.c" `
    "$SOLVER\solver_code\external\amd\src\amd_global.c" `
    "$SOLVER\solver_code\external\amd\src\amd_info.c" `
    "$SOLVER\solver_code\external\amd\src\amd_order.c" `
    "$SOLVER\solver_code\external\amd\src\amd_post_tree.c" `
    "$SOLVER\solver_code\external\amd\src\amd_postorder.c" `
    "$SOLVER\solver_code\external\amd\src\amd_preprocess.c" `
    "$SOLVER\solver_code\external\amd\src\amd_valid.c" `
    "$SOLVER\solver_code\external\ldl\src\ldl.c" `
    -I "$OUT" `
    -I "$SOLVER\include" `
    -I "$SOLVER\solver_code\include" `
    -I "$SOLVER\solver_code\external\amd\include" `
    -I "$SOLVER\solver_code\external\ldl\include" `
    -I "$SOLVER\solver_code\external\SuiteSparse_config" `
    -o "$OUT\astrobee_cogu.exe" `
    -lm

if ($LASTEXITCODE -eq 0) {
    Write-Host "`n=== COMPILACION EXITOSA ===" -ForegroundColor Green
    Write-Host "Ejecutable: $OUT\astrobee_cogu.exe"
} else {
    Write-Host "`n=== ERROR DE COMPILACION (codigo $LASTEXITCODE) ===" -ForegroundColor Red
}
