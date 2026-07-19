@echo off
setlocal
cd /d "%~dp0"

where py >nul 2>nul
if %errorlevel%==0 (
  py -3 start_local.py
  exit /b %errorlevel%
)

where python >nul 2>nul
if %errorlevel%==0 (
  python start_local.py
  exit /b %errorlevel%
)

echo CallVCF requires Python 3.8 or newer.
echo Download Python from https://www.python.org/downloads/windows/
pause
exit /b 1
