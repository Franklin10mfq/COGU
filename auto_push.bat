@echo off
setlocal

REM ==============================
REM CONFIGURACION
REM ==============================

REM Ruta absoluta de tu repositorio
set REPO_PATH=C:\Users\franklin\Documents\Kyutech\Thesis\Papers\githubs\github_COGU

REM Mensaje de commit (puede cambiarse al ejecutar)
set COMMIT_MSG=%1

REM Si no se pasa mensaje como argumento
if "%COMMIT_MSG%"=="" (
    set COMMIT_MSG=Auto commit
)

REM ==============================
REM EJECUCION
REM ==============================

echo.
echo ==============================
echo Going to repository...
echo ==============================

cd /d %https://github.com/Franklin10mfq/COGU%

if not exist ".git" (
    echo ERROR: This is not a git repository.
    pause
    exit /b
)

echo.
echo ==============================
echo Git Add
echo ==============================
git add .

echo.
echo ==============================
echo Git Commit
echo ==============================
git commit -m "%COMMIT_MSG%"

echo.
echo ==============================
echo Git Push
echo ==============================
git push

echo.
echo Done 🚀
pause
