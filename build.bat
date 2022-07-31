echo set current directory
cd /d %~dp0

git clone https://github.com/bulletphysics/bullet3.git
cd bullet3
mkdir build
cd build
cmake .. -DBUILD_CPU_DEMOS=OFF -DBUILD_BULLET3=OFF -DBUILD_ENET=OFF -DBUILD_CLSOCKET=OFF -DBUILD_BULLET2_DEMOS=OFF -DBUILD_UNIT_TESTS=OFF -DINSTALL_CMAKE_FILES=OFF -DINSTALL_LIBS=ON
cmake --build . --config Release
cmake --install . --prefix .

cd ../..
git clone https://github.com/AcademySoftwareFoundation/openvdb.git -b v2.2.0
cd openvdb
mkdir build
cd build
