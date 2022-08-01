include(ExternalProject)

# eigen
list(APPEND THIRDPARTY_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}/eigen)

# tclap
list(APPEND THIRDPARTY_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}/tclap/include)

# boost

set(Boost_Bootstrap_Command ${CMAKE_CURRENT_SOURCE_DIR}/boost/bootstrap.bat)
set(Boost_b2_Command ${CMAKE_CURRENT_SOURCE_DIR}/boost/b2.exe install
        --prefix=${CMAKE_BINARY_DIR}/Install/boost
        --threading=multi
        --link=static
        runtime-link=static
)
ExternalProject_Add(EP_boost
        SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/boost
        BUILD_IN_SOURCE   1
        CONFIGURE_COMMAND ${Boost_Bootstrap_Command}
        BUILD_COMMAND     ${Boost_b2_Command}
        INSTALL_COMMAND   ""
        PATCH_COMMAND     ""
        UPDATE_COMMAND    ""
        INSTALL_DIR       ${CMAKE_BINARY_DIR}/Install/boost
)

list(APPEND THIRDPARTY_INCLUDES ${CMAKE_CURRENT_BINARY_DIR}/Install/boost/include)

# openvdb

ExternalProject_Add(
    EP_openvdb
    SOURCE_DIR      ${CMAKE_CURRENT_SOURCE_DIR}/openvdb/openvdb
    PREFIX          ${CMAKE_CURRENT_BINARY_DIR}/Build/openvdb
    INSTALL_DIR     ${CMAKE_CURRENT_BINARY_DIR}/Install/openvdb
    CMAKE_ARGS      -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}/Install/openvdb
                    -DCMAKE_USER_MAKE_RULES_OVERRIDE=${CMAKE_CURRENT_SOURCE_DIR}/CFlagOverrides.cmake
                    -DUSE_MSVC_RUNTIME_LIBRARY_DLL=${USE_MSVC_RUNTIME_LIBRARY_DLL}
                    -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE}
                    -DCMAKE_DEBUG_POSTFIX=d
)


list(APPEND THIRDPARTY_INCLUDES ${CMAKE_CURRENT_BINARY_DIR}/Install/openvdb/include)
list(APPEND THIRDPARTY_LIBRARY_DIRECTORIES ${CMAKE_CURRENT_BINARY_DIR}/Install/openvdb/lib)

# bullet3

ExternalProject_Add(
    EP_bullet3
    SOURCE_DIR      ${CMAKE_CURRENT_SOURCE_DIR}/bullet3
    PREFIX          ${CMAKE_CURRENT_BINARY_DIR}/Build/bullet3
    INSTALL_DIR     ${CMAKE_CURRENT_BINARY_DIR}/Install/bullet3
    CMAKE_ARGS      -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}/Install/bullet3
                    -DCMAKE_USER_MAKE_RULES_OVERRIDE=${CMAKE_CURRENT_SOURCE_DIR}/CFlagOverrides.cmake
                    -DCMAKE_DEBUG_POSTFIX=d
                    -DUSE_MSVC_RUNTIME_LIBRARY_DLL=${USE_MSVC_RUNTIME_LIBRARY_DLL}
                    -DBUILD_SHARED_LIBS=OFF
                    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                    -DCMAKE_POSITION_INDEPENDENT_CODE=ON
                    -DUSE_GRAPHICAL_BENCHMARK=OFF
                    -DUSE_SOFT_BODY_MULTI_BODY_DYNAMICS_WORLD=OFF
                    -DENABLE_VHACD=OFF
                    -DBULLET2_MULTITHREADING=${BULLET_MULTITHREADING}
                    -DBULLET2_USE_OPEN_MP_MULTITHREADING=${BULLET_MULTITHREADING}
                    -DBUILD_CPU_DEMOS=OFF
                    -DUSE_GLUT=OFF
                    -DBUILD_BULLET3=OFF
                    -DBUILD_ENET=OFF
                    -DBUILD_CLSOCKET=OFF
                    -DBUILD_EGL=OFF
                    -DBUILD_OPENGL3_DEMOS=OFF
                    -DBUILD_BULLET2_DEMOS=OFF
                    -DBUILD_EXTRAS=ON
                    -DBUILD_UNIT_TESTS=OFF
                    -DINSTALL_CMAKE_FILES=OFF
                    -DINSTALL_LIBS=ON
                    -DUSE_MSVC_DISABLE_RTTI=OFF
)

list(APPEND THIRDPARTY_INCLUDES ${CMAKE_CURRENT_BINARY_DIR}/Install/bullet3/include/bullet)
list(APPEND THIRDPARTY_LIBRARY_DIRECTORIES ${CMAKE_CURRENT_BINARY_DIR}/Install/bullet3/lib)

# openexr

ExternalProject_Add(
    EP_openexr
    SOURCE_DIR      ${CMAKE_CURRENT_SOURCE_DIR}/openexr
    PREFIX          ${CMAKE_CURRENT_BINARY_DIR}/Build/openexr
    INSTALL_DIR     ${CMAKE_CURRENT_BINARY_DIR}/Install/openexr
    CMAKE_ARGS      -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}/Install/openexr
                    -DCMAKE_USER_MAKE_RULES_OVERRIDE=${CMAKE_CURRENT_SOURCE_DIR}/CFlagOverrides.cmake
                    -DUSE_MSVC_RUNTIME_LIBRARY_DLL=${USE_MSVC_RUNTIME_LIBRARY_DLL}
                    -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE}
                    -DCMAKE_DEBUG_POSTFIX=d
                    -DPYILMBASE_ENABLE=OFF
)

list(APPEND THIRDPARTY_INCLUDES ${CMAKE_CURRENT_BINARY_DIR}/Install/openexr/include)
list(APPEND THIRDPARTY_LIBRARY_DIRECTORIES ${CMAKE_CURRENT_BINARY_DIR}/Install/openexr/lib)

# tbb

ExternalProject_Add(
    EP_tbb
    SOURCE_DIR      ${CMAKE_CURRENT_SOURCE_DIR}/tbb
    PREFIX          ${CMAKE_CURRENT_BINARY_DIR}/Build/tbb
    INSTALL_DIR     ${CMAKE_CURRENT_BINARY_DIR}/Install/tbb
    CMAKE_ARGS      -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}/Install/tbb
                    -DCMAKE_USER_MAKE_RULES_OVERRIDE=${CMAKE_CURRENT_SOURCE_DIR}/CFlagOverrides.cmake
                    -DUSE_MSVC_RUNTIME_LIBRARY_DLL=${USE_MSVC_RUNTIME_LIBRARY_DLL}
                    -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE}
                    -DCMAKE_DEBUG_POSTFIX=d
                    -DTBB_BUILD_TESTS=OFF
)

list(APPEND THIRDPARTY_INCLUDES ${CMAKE_CURRENT_BINARY_DIR}/Install/tbb/include)
list(APPEND THIRDPARTY_LIBRARY_DIRECTORIES ${CMAKE_CURRENT_BINARY_DIR}/Install/tbb/lib)

# zlib

ExternalProject_Add(
    EP_zlib
    SOURCE_DIR      ${CMAKE_CURRENT_SOURCE_DIR}/zlib
    PREFIX          ${CMAKE_CURRENT_BINARY_DIR}/Build/zlib
    INSTALL_DIR     ${CMAKE_CURRENT_BINARY_DIR}/Install/zlib
    CMAKE_ARGS      -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}/Install/zlib
                    -DCMAKE_USER_MAKE_RULES_OVERRIDE=${CMAKE_CURRENT_SOURCE_DIR}/CFlagOverrides.cmake
                    -DUSE_MSVC_RUNTIME_LIBRARY_DLL=${USE_MSVC_RUNTIME_LIBRARY_DLL}
                    -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE}
                    -DCMAKE_DEBUG_POSTFIX=d
)

list(APPEND THIRDPARTY_INCLUDES ${CMAKE_CURRENT_BINARY_DIR}/Install/zlib/include)
list(APPEND THIRDPARTY_LIBRARY_DIRECTORIES ${CMAKE_CURRENT_BINARY_DIR}/Install/zlib/lib)

# assign valiables into the parent scope
set(THIRDPARTY_INCLUDES ${THIRDPARTY_INCLUDES} PARENT_SCOPE)
set(THIRDPARTY_LIBRARY_DIRECTORIES ${THIRDPARTY_LIBRARY_DIRECTORIES} PARENT_SCOPE)