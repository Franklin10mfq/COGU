$env:PATH = "C:\msys64\ucrt64\bin;C:\msys64\usr\bin;" + $env:PATH
$gcc = "C:\msys64\ucrt64\bin\gcc.exe"

$base   = "c:\Users\luiss\COGU\COGU_algorithms\Astrobee_template"
$solver = "$base\Astrobee_SCVX_ECOS_T51_l2\c"
$amd    = "$solver\solver_code\external\amd\src"

$files = @(
  "$base\ASTROBEE_TEMPLATE_C\code_c_ZOH.c",
  "$solver\src\cpg_solve.c",
  "$solver\src\cpg_workspace.c",
  "$solver\solver_code\src\ecos.c",
  "$solver\solver_code\src\kkt.c",
  "$solver\solver_code\src\cone.c",
  "$solver\solver_code\src\spla.c",
  "$solver\solver_code\src\timer.c",
  "$solver\solver_code\src\preproc.c",
  "$solver\solver_code\src\splamm.c",
  "$solver\solver_code\src\ctrlc.c",
  "$solver\solver_code\src\equil.c",
  "$solver\solver_code\src\expcone.c",
  "$solver\solver_code\src\wright_omega.c",
  "$amd\amd_1.c", "$amd\amd_2.c", "$amd\amd_aat.c",
  "$amd\amd_control.c", "$amd\amd_defaults.c", "$amd\amd_dump.c",
  "$amd\amd_global.c", "$amd\amd_info.c", "$amd\amd_order.c",
  "$amd\amd_post_tree.c", "$amd\amd_postorder.c", "$amd\amd_preprocess.c",
  "$amd\amd_valid.c",
  "$solver\solver_code\external\ldl\src\ldl.c"
)

$iflags = @(
  "-I$solver\include",
  "-I$solver\solver_code\include",
  "-I$solver\solver_code\external\SuiteSparse_config",
  "-I$solver\solver_code\external\amd\include",
  "-I$solver\solver_code\external\ldl\include"
)

$out = "$base\ASTROBEE_TEMPLATE_C\astrobee_cogu.exe"

Write-Output "Compilando... (puede tardar ~30 segundos)"
$output = & $gcc @("-O0", "-std=c99", "-D_USE_MATH_DEFINES") @files @iflags -o $out -lm 2>&1
if ($output) { Write-Output "Mensajes del compilador:"; $output }
Write-Output "Exit: $LASTEXITCODE"

if (Test-Path $out) {
    Write-Output "`nCOMPILACION EXITOSA: $out"
} else {
    Write-Output "`nCOMPILACION FALLO - no se creo el ejecutable"
}
