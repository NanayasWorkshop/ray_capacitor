@echo off
echo Building Model Dimension Checker...
set PATH=C:\mingw64\bin;%PATH%
g++ -std=c++11 src\model_dimension_checker.cpp -o bin\model_dimension_checker.exe
echo Done!
pause