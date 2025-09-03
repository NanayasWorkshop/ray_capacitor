@echo off
echo Building Dynamic Capacitor with MinGW...
set EMBREE_PATH=C:\embree
set PATH=C:\mingw64\bin;%PATH%
g++ -std=c++11 -I"%EMBREE_PATH%\include" src\dynamic_capacitor.cpp -L"%EMBREE_PATH%\lib" -lembree4 -o bin\dynamic_capacitor.exe
echo Done!
pause