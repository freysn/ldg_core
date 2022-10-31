# https://cliutils.gitlab.io/modern-cmake/chapters/basics.html

cmake_minimum_required(VERSION 3.14)

project(ldg_core VERSION 0.1
  DESCRIPTION "Level-of-Detail Grids"
  LANGUAGES CXX)

# ~/dev/tmp/modern-cmake/examples/fetch/CMakeLists.txt

include(FetchContent)

FetchContent_Declare(
  helper
  GIT_REPOSITORY git@github.com:freysn/helper.git
  GIT_TAG master
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
)


FetchContent_MakeAvailable(helper)



# There are several things to unpack here. one is both the name of the executable file generated, and the name of the CMake target created (you'll hear a lot more about targets soon, I promise). The source file list comes next, and you can list as many as you'd like. CMake is smart, and will only compile source file extensions. The headers will be, for most intents and purposes, ignored; the only reason to list them is to get them to show up in IDEs. Targets show up as folders in many IDEs. More about the general build system and targets is available at buildsystem.

add_executable(ldg_core supertiles_place.cpp)

# PUBLICdoesn't mean much for an executable; for a library it lets CMake know that any targets that link to this target must also need that include directory.
target_include_directories(ldg_core PUBLIC .)

# target_link_libraries(another PUBLIC one)
target_compile_features(ldg_core PUBLIC cxx_std_17)


find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(ldg_core PUBLIC OpenMP::OpenMP_CXX)
else()
    target_compile_definitions(ldg_core PUBLIC NO_OMP=1)
endif()


find_package(PNG REQUIRED)
include_directories(${PNG_INCLUDE_DIR})
target_link_libraries(ldg_core PUBLIC ${PNG_LIBRARY})


add_library(helper INTERFACE)
target_include_directories(helper INTERFACE ${helper_SOURCE_DIR})

target_link_libraries(ldg_core PUBLIC helper)


find_package(BZip2 REQUIRED)
include_directories(${BZIP2_INCLUDE_DIRS})
target_link_libraries(ldg_core PUBLIC ${BZIP2_LIBRARIES})

# find_package(CAIRO REQUIRED)
# include_directories(${CAIRO_INCLUDE_DIR})
# target_link_libraries(ldg_core PUBLIC ${CAIRO_LIBRARY})

#find_package(CIMG REQUIRED)
# include_directories(${PNG_INCLUDE_DIR})

# the `pkg_check_modules` function is created with this call
find_package(PkgConfig REQUIRED) 

# these calls create special `PkgConfig::<MODULE>` variables
pkg_check_modules(MY_PKG REQUIRED IMPORTED_TARGET cairo)

target_link_libraries(ldg_core PUBLIC PkgConfig::MY_PKG)


# Detect if the "port" command is valid on this system; if so, return full path
EXECUTE_PROCESS(COMMAND which port RESULT_VARIABLE DETECT_MACPORTS OUTPUT_VARIABLE MACPORTS_PREFIX ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)

IF (${DETECT_MACPORTS} EQUAL 0)
    # "/opt/local/bin/port" doesn't have libs, so we get the parent directory
    GET_FILENAME_COMPONENT(MACPORTS_PREFIX ${MACPORTS_PREFIX} DIRECTORY)

    # "/opt/local/bin" doesn't have libs, so we get the parent directory
    GET_FILENAME_COMPONENT(MACPORTS_PREFIX ${MACPORTS_PREFIX} DIRECTORY)

    # "/opt/local" is where MacPorts lives, add `/lib` suffix and link
    LINK_DIRECTORIES(${LINK_DIRECTORIES} ${MACPORTS_PREFIX}/lib)

    INCLUDE_DIRECTORIES(${INCLUDE_DIRECTORIES} ${MACPORTS_PREFIX}/include)

    MESSAGE("WINNING!: ${MACPORTS_PREFIX}/lib")
ENDIF()
