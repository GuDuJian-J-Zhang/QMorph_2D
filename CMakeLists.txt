cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
# do not warn non-definition of WIN32 on cygwin
set(CMAKE_LEGANCY_CYGWIN_WIN32 0)

if(DEFINED CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Choose build type")
else(DEFINED CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose build type")
endif()

project(qmorph2D)

set(QMORPH_SRC
    #src/main.cpp
    src/predicates.cxx
	src/QuadMesher.cpp
	src/sldelaunay2d.cxx
	src/Vector3D.cpp
	src/TriMergeGenerator.cpp
	src/QMorphGenerator.cpp
	src/quadMeshAPI.cpp
	src/MeshOptimizer.cpp
)

set(QMORPH_INC
    inc/foundation.h
    inc/predicates.h
    inc/QuadMesher.h
    inc/sldelaunay2d.h
    inc/sldelaunay2d_std.h
    inc/sldelaunaymshtools.hpp
    inc/slharray1d.h
    inc/toolsdefine.h
    inc/USE_CGAL.h
    inc/Vector3D.h
    inc/QuadMeshGenerator.h
	inc/QMorphGenerator.h
	inc/TriMergeGenerator.h
	inc/MeshOptimizer.h
	quadMeshAPI2.h
)

macro(find_all_libraries VARNAME LISTNAME PATH SUFFIX)
  set(${VARNAME})
  list(LENGTH ${LISTNAME} NUM_LIST)
  foreach(LIB ${${LISTNAME}})
    if("${PATH}" STREQUAL "")
      find_library(FOUND_LIB ${LIB} PATH_SUFFIXES ${SUFFIX})
    else("${PATH}" STREQUAL "")
      find_library(FOUND_LIB ${LIB} PATHS ${PATH} NO_DEFAULT_PATH)
    endif("${PATH}" STREQUAL "")
    if(FOUND_LIB)
      list(APPEND ${VARNAME} ${FOUND_LIB})
    endif(FOUND_LIB)
    unset(FOUND_LIB CACHE)
  endforeach(LIB)
  list(LENGTH ${VARNAME} NUM_FOUND_LIBRARIES)
  if(NUM_FOUND_LIBRARIES LESS NUM_LIST)
    set(${VARNAME})
  endif(NUM_FOUND_LIBRARIES LESS NUM_LIST)
endmacro()

if(MSVC)
    set(DIR_3RD_PARTY "$ENV{DIR_3RDPARTY}")
    set(BOOST_DIR "${DIR_3RD_PARTY}/boost_1_59_0")
    set(CGAL_DIR "${DIR_3RD_PARTY}/CGAL-4.8")
    set(OPENGL_DIR "${DIR_3RD_PARTY}/OpenGL")
elseif(UNIX)
    set(DIR_3RD_PARTY "/home/gudujian/gudujian/software/3RD_PARTY")
    set(BOOST_DIR "${DIR_3RD_PARTY}/boost_1_62_0")
    set(CGAL_DIR "${DIR_3RD_PARTY}/CGAL")
    #set(OPENGL_DIR "${DIR_3RD_PARTY}/OpenGL")
endif()

if(MSVC)
    include_directories("${BOOST_DIR}")
    include_directories("${CGAL_DIR}/include")
    include_directories("${CGAL_DIR}/auxiliary/gmp/include")
    include_directories("${CGAL_DIR}/build/include")
elseif(UNIX)
    include_directories("${BOOST_DIR}")
    include_directories("${CGAL_DIR}/include")
endif()
	#include_directories("${OPENGL_DIR}/include")
include_directories(inc)

if(MSVC)
    link_directories("${BOOST_DIR}/lib64-msvc-11.0")
	link_directories("${CGAL_DIR}/build/lib")
	link_directories("${CGAL_DIR}/auxiliary/gmp/lib")
	#link_directories("${OPENGL_DIR}/lib")
elseif(UNIX)
    add_definitions(-DLINUX)
	set(BOOST_LIB_DIR "${BOOST_DIR}/lib")
	set(BOOST_LIBS_REQUIRED
        boost_filesystem boost_thread boost_system boost_timer boost_chrono boost_python)
    find_all_libraries(BOOST_LIBS BOOST_LIBS_REQUIRED "${BOOST_LIB_DIR}" "")
	set(CGAL_LIB_DIR "${CGAL_DIR}/lib")
	link_directories("${BOOST_LIB_DIR}")
	link_directories("${CGAL_LIB_DIR}/Debug")
endif()

if(MSVC)
    #set(OPENGL_LIBS_REQUIRED opengl32 glu32 glut32 glut64)

    #list(APPEND EXTERNAL_LIBRARIES_DEBUG "${OPENGL_DIR}/lib/freeglut.lib")
	#list(APPEND EXTERNAL_LIBRARIES_DEBUG "${OPENGL_DIR}/lib/glew32.lib")
	#list(APPEND EXTERNAL_LIBRARIES_OPT "${OPENGL_DIR}/lib/glut64.lib")
endif()

#if(IS_EXE)
#    add_executable(quadMesh ${QMORPH_SRC} ${QMORPH_INC})
#else()
#    add_library(quadMesh SHARED ${QMORPH_SRC} ${QMORPH_INC})
#endif()

add_library(quadMesh SHARED ${QMORPH_SRC} ${QMORPH_INC})
#add_executable(qmorph ${QMORPH_SRC} ${QMORPH_INC})

if(MSVC)
    set_target_properties(quadMesh PROPERTIES COMPILE_FLAGS "/bigobj")
elseif(UNIX)
    set_target_properties(quadMesh PROPERTIES COMPILE_FLAGS "-std=c++11")
endif()

target_link_libraries(quadMesh ${EXTERNAL_LIBRARIES} debug  "${EXTERNAL_LIBRARIES_DEBUG}" optimized  "${EXTERNAL_LIBRARIES_OPT}")
