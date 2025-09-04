@echo off
echo Building Multi-Sensor Capacitor with MinGW...
set EMBREE_PATH=C:\embree
set PATH=C:\mingw64\bin;%PATH%
g++ -std=c++11 -I"%EMBREE_PATH%\include" src\multi_sensor_capacitor.cpp -L"%EMBREE_PATH%\lib" -lembree4 -o bin\multi_sensor_capacitor.exe
echo Done!
pause