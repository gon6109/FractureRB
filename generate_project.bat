echo set current directory
cd /d %~dp0

if not exist ..\vcpkg\vcpkg.exe (
    cd ..
    git clone https://github.com/microsoft/vcpkg -b 2022.06.16.1
    cd FractureRB
    ..\vcpkg\bootstrap-vcpkg.bat
)

..\vcpkg\vcpkg.exe install eigen3:x64-windows tclap:x64-windows boost:x64-windows tbb:x64-windows zlib:x64-windows openexr:x64-windows

git submodule update --init

if not exist build mkdir build
cd build
cmake .. -DCMAKE_TOOLCHAIN_FILE=../../vcpkg/scripts/buildsystems/vcpkg.cmake