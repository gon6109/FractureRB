cmake_minimum_required(VERSION 3.21)

if(MSVC)
	cmake_policy(SET CMP0091 NEW)
endif()

project(FractureRB)

if(MSVC)
	option(USE_MSVC_RUNTIME_LIBRARY_DLL "Bulid as MultithreadedDLL" ON)
endif()

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
	set(CMAKE_INSTALL_PREFIX ".." CACHE PATH "comment" FORCE)
endif()

set(CMAKE_CXX_STANDARD 14)

# specify project type
set(CMAKE_CONFIGURATION_TYPES Debug Release)

# specify destionation directories
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG "${CMAKE_BINARY_DIR}/DEBUG")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_DEBUG "${CMAKE_BINARY_DIR}/DEBUG")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE "${CMAKE_BINARY_DIR}/RELEASE")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_RELEASE "${CMAKE_BINARY_DIR}/RELEASE")

# specify msvc runtime library
if(MSVC)
	if(USE_MSVC_RUNTIME_LIBRARY_DLL)
		cmake_policy(SET CMP0091 NEW)
		set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>DLL")
	else()
		cmake_policy(SET CMP0091 NEW)
		set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
	endif()

	message(${CMAKE_MSVC_RUNTIME_LIBRARY})
endif()

option(BUILD_HYENA "If on, HyENA library parts are built from source, otherwise link to a pre-built version" OFF)

if(BUILD_HYENA)
	set(LIBTYPE SHARED)
endif()

option(BULLET_DOUBLE_PRECISION "If on, we'll link to double-precision Bullet libs; recommended" OFF)

add_definitions(-DASSEMBLE_GALERKIN_OMP) # use OpenMP in HyENA to assemble SGBEM matrices

if(MSVC)
	# stuff for Visual Studio compiler
	add_definitions(-TP) # force MSVC to use C++ compiler, otherwise it falls back to ANSI C
	add_definitions(-D_USE_MATH_DEFINES) # make sure we have M_PI available
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /openmp") # enable OpenMP
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W0") # get rid of warnings

# add_definitions(-D_UNICODE -DUNICODE)
else()
	# stuff for GNU compiler
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp") # enable OpenMP
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w") # get rid of warnings

	if(BUILD_HYENA)
		SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -fopenmp")
	endif()
endif()

# my Header-only libraries are all in one directory, find that dir by searching for Eigen/Eigen
find_path(HLIB_INC Eigen/Eigen)
include_directories(${HLIB_INC}) # takes care of other header-only libs, need Eigen, TCLAP, and parts of Boost

# the HyENA library is in the project dir, because we've made a few adjustments to it
include_directories(${PROJECT_SOURCE_DIR})

# and here's where all the other headers are
include_directories(${PROJECT_SOURCE_DIR}/src)

add_subdirectory(thirdparty)

include_directories(${THIRDPARTY_INCLUDES})
link_directories(${THIRDPARTY_LIBRARY_DIRECTORIES})

# OpenVDB libraries and includes
add_definitions(-DOPENVDB_STATICLIB) # use static OpenVDB lib
find_library(tbblib tbb)
find_library(tbblibd tbb_debug)
find_package(IlmBase)
find_package(ZLIB)
set(OpenVDB_libs
	debug ${tbblibd}
	optimized ${tbblib}
	ZLIB::ZLIB
	IlmBase::Half
	debug openvdbd
	optimized openvdb)

# Bullet libraries and includes
if(BULLET_DOUBLE_PRECISION)
	add_definitions(-DBT_USE_DOUBLE_PRECISION)
endif()

# Bullet core libraries
set(BULLET_LIBS
	debug BulletDynamicsd
	optimized BulletDynamics
	debug BulletCollisiond
	optimized BulletCollision
	debug LinearMathd
	optimized LinearMath) # order is very important for GCC builds! (MSVC is ok either way)

# Bullet Extras libs
set(BULLET_EXT_LIBS
	debug BulletWorldImporterd
	optimized BulletWorldImporter
	debug BulletFileLoaderd
	optimized BulletFileLoader)

if(BUILD_HYENA)
	add_library(HyENAlib2 ${LIBTYPE}
		"hyena/wrap/HyENAWrapper.cpp"
		"hyena/wrap/QuadratureWrapper.cpp"
	)
	find_package(FFTW3)
	target_link_libraries(HyENAlib2 PRIVATE FFTW3::fftw3)
	install(TARGETS HyENAlib2 RUNTIME DESTINATION bin LIBRARY DESTINATION bin ARCHIVE DESTINATION bin) # throw this in the same dir as the binaries
endif()

if(MSVC)
	add_compile_options(/bigobj)
endif()

add_library(FractureLib STATIC
	"src/FractureBEM.cpp"
	"src/FractureModel.cpp"
	"src/Reader.cpp"
	"src/Reader_VCG.cpp"
	"src/PostProcessor.cpp"
	"src/VDBWrapper.cpp"
	"src/VDBWrapper_mesh.cpp"
	"src/SubsampledCrackTip.cpp"
	"src/MaterialModel.cpp"
	"vcg/wrap/ply/plylib.cpp"
)

if(BUILD_HYENA)
	target_link_libraries(FractureLib HyENAlib2)
else()
	add_library(HyENAlib2 SHARED IMPORTED)
	set_property(TARGET HyENAlib2 APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
	set_property(TARGET HyENAlib2 APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
	set_target_properties(HyENAlib2 PROPERTIES IMPORTED_LOCATION ${CMAKE_CURRENT_SOURCE_DIR}/hyena/HyENAlib2_x64_MSVC_2010/HyENAlib2.dll)
	set_target_properties(HyENAlib2 PROPERTIES IMPORTED_IMPLIB ${CMAKE_CURRENT_SOURCE_DIR}/hyena/HyENAlib2_x64_MSVC_2010/HyENAlib2.lib)
	target_link_libraries(FractureLib HyENAlib2)
endif()

add_dependencies(FractureLib EP_openvdb)
install(TARGETS FractureLib RUNTIME DESTINATION bin LIBRARY DESTINATION bin ARCHIVE DESTINATION bin) # throw this in the same dir as the binaries

FILE(GLOB SRC_HEADERS ${PROJECT_SOURCE_DIR}/src/*.h) # for convenience -> puts the headers into the MSVC project list
add_executable(FractureRB ${SRC_HEADERS}
	"src/main.cpp"
	"src/BulletWrapper.cpp"
	"src/FractureRB.cpp"
	"src/FractureRB_fragments.cpp"
	"src/ColliderData.cpp"
)
target_link_libraries(FractureRB FractureLib ${OpenVDB_libs} ${BULLET_EXT_LIBS} ${BULLET_LIBS}) # order is important for GCC
add_dependencies(FractureRB EP_openvdb EP_bullet3)
install(TARGETS FractureRB RUNTIME DESTINATION bin LIBRARY DESTINATION lib ARCHIVE DESTINATION lib)

add_executable(FractureRB_vis ${SRC_HEADERS}
	"src/mainVis.cpp"
	"src/VDBLoader.cpp"
)
target_link_libraries(FractureRB_vis FractureLib ${OpenVDB_libs})
add_dependencies(FractureRB_vis EP_openvdb)
add_compile_definitions(FracutreRB_vis _SILENCE_STDEXT_HASH_DEPRECATION_WARNINGS)
install(TARGETS FractureRB_vis RUNTIME DESTINATION bin LIBRARY DESTINATION lib ARCHIVE DESTINATION lib)

# old FractureBEM build targets -- mostly used for testing
#
add_executable(FractureBEM ${SRC_HEADERS}
	"src/mainFractureBEM.cpp"
)
target_link_libraries(FractureBEM FractureLib ${OpenVDB_libs})
add_dependencies(FractureBEM EP_openvdb)
install(TARGETS FractureBEM RUNTIME DESTINATION bin LIBRARY DESTINATION lib ARCHIVE DESTINATION lib)

#
add_executable(SegmentVDB
	"src/mainSegment.cpp"
	"src/mySegment.h"
)
target_link_libraries(SegmentVDB ${OpenVDB_libs})
add_dependencies(SegmentVDB EP_openvdb)
install(TARGETS SegmentVDB RUNTIME DESTINATION bin LIBRARY DESTINATION lib ARCHIVE DESTINATION lib)

#
add_executable(FractureBEM_vis ${SRC_HEADERS}
	"src/mainVisFractureBEM.cpp"
	"src/VDBLoader.cpp"
)
target_link_libraries(FractureBEM_vis FractureLib ${OpenVDB_libs})
add_dependencies(FractureBEM_vis EP_openvdb)
install(TARGETS FractureBEM_vis RUNTIME DESTINATION bin LIBRARY DESTINATION lib ARCHIVE DESTINATION lib)

add_executable(GenerateFractureData ${SRC_HEADERS}
	"src/mainGenerateFractureData.cpp"
	"src/FractureDataGenerator.cpp"
	"src/FractureRB.cpp"
	"src/FractureRB_fragments.cpp"
	"src/ColliderData.cpp"
)
target_link_libraries(GenerateFractureData FractureLib ${OpenVDB_libs} ${BULLET_EXT_LIBS} ${BULLET_LIBS})
add_dependencies(GenerateFractureData EP_bullet3 EP_openvdb)
install(TARGETS GenerateFractureData RUNTIME DESTINATION bin LIBRARY DESTINATION lib ARCHIVE DESTINATION lib)

if(MSVC)
	set_target_properties(FractureRB PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/$(Configuration)")
	set_target_properties(FractureRB_vis PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/$(Configuration)")
	set_target_properties(FractureBEM PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/$(Configuration)")
	set_target_properties(FractureBEM_vis PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/$(Configuration)")
	set_target_properties(GenerateFractureData PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/$(Configuration)")
endif()