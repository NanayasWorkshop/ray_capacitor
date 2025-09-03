@echo off
echo Building with MinGW...
set EMBREE_PATH=C:\embree
set PATH=C:\mingw64\bin;%PATH%
g++ -std=c++11 -I"%EMBREE_PATH%\include" src\ray_distance.cpp -L"%EMBREE_PATH%\lib" -lembree4 -o bin\ray_distance.exe
echo Done!
pause