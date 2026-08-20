@echo off
chcp 65001 >nul
title Spectral Splitting - Reproducibility Run
where python >nul 2>nul && (set "PY=python") || (set "PY=py")

echo [1/3] Checking/installing dependencies (first run only)...
%PY% -m pip install -r requirements.txt

echo [2/3] Running main.py - takes tens of minutes, please wait...
cd code
%PY% -u main.py
cd ..

echo [3/3] Finished. Checking outputs:
dir /b code\results
echo.
echo If you see 4 .xlsx and 2 .png files above, reproduction succeeded.
pause
