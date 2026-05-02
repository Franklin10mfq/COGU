# Compilar Astrobee T=30 desde COGU_Dev_Library
# Correr desde: c:\Users\luiss\COGU\

$env:PATH = "C:\msys64\ucrt64\bin;C:\msys64\usr\bin;" + $env:PATH

$SOLVER   = "COGU_Dev_Library\c_output_t30\solver\c"
$TEMPLATE = "COGU_Dev_Library\astrobee_t30"

gcc -O0 -std=c99 -D_USE_MATH_DEFINES `
    "$TEMPLATE\code_c_ZOH.c" `
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
    -I "$SOLVER\include" `
    -I "$SOLVER\solver_code\include" `
    -I "$SOLVER\solver_code\external\amd\include" `
    -I "$SOLVER\solver_code\external\ldl\include" `
    -I "$SOLVER\solver_code\external\SuiteSparse_config" `
    -I "$TEMPLATE" `
    -o "$TEMPLATE\astrobee_cogu_t30.exe" `
    -lm

if ($LASTEXITCODE -eq 0) {
    Write-Host "`n=== COMPILACION EXITOSA ===" -ForegroundColor Green
    Write-Host "Ejecutable: $TEMPLATE\astrobee_cogu_t30.exe"
} else {
    Write-Host "`n=== ERROR DE COMPILACION (codigo $LASTEXITCODE) ===" -ForegroundColor Red
}
