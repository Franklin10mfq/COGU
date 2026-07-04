To generate the .lib:

$files = Get-ChildItem -Filter *.c | ForEach-Object { $_.Name }; gcc -shared -o cpg_example.dll "-Wl,--out-implib,cpg_example.lib" -O2 $files