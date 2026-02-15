@echo off
setlocal EnableDelayedExpansion

REM ==============================
REM CONFIGURACION
REM ==============================

set REPO_PATH=C:\Users\franklin\Documents\Kyutech\Thesis\Papers\githubs\github_COGU

REM Captura TODOS los argumentos
set COMMIT_MSG=%*

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

cd /d "%https://github.com/Franklin10mfq/COGU%"

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