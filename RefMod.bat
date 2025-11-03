
@echo off
setlocal

echo Attempting to launch RefMod...
REM Get the current directory
set "BAT_DIR=%~dp0"
python "%BAT_DIR%src\RefMod_GUI.py"
if %errorlevel% equ 0 (
    goto :eof
)

REM echo Launch failed. Checking Python and pip...
REM python -m ensurepip >nul 2>&1
REM if %errorlevel% neq 0 (
REM     echo Python or pip is not installed. Please install Python first.
REM     pause
REM     exit /b
REM )

echo Launch failed. Checking Python installation...
where python >nul 2>&1
if %errorlevel% neq 0 (
    echo Python is not installed or not in PATH. Please install Python.
    pause
    exit /b
)

echo Checking pip installation...
python -m pip --version >nul 2>&1
if %errorlevel% neq 0 (
    echo Pip not found. Attempting to install pip...
    python -m ensurepip
)


echo Installing requirements...
pip install -r requirements.txt

echo Trying to launch RefMod again...
python src\RefMod_GUI.py
if %errorlevel% equ 0 (
    goto :eof
)

echo Failed to launch RefMod.
pause
exit /b
