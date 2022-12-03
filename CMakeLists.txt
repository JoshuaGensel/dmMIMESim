cmake_minimum_required(VERSION 3.10)
# soll man eigentlich nicht, aber ist nötig für filesystem

project(MIMESim CXX)

#Option to run Unit Tests
option(MIMESimTest_WITH_COV "Enable profiling and coverage report analysis" OFF)
option(UNIT_TESTING "Enable Unit Testing" ON)
SET(CMAKE_CODEBLOCKS_EXCLUDE_EXTERNAL_FILES OFF)

set(MIMESim_VERSION_MAJOR 0)
set(MIMESim_VERSION_MINOR 1)
set(MIMESim_VERSION_PATCH 0)
set(MIMESim_VERSION "${MIMESim_VERSION_MAJOR}.${MIMESim_VERSION_MINOR}.${MIMESim_VERSION_PATCH}")

set(CMAKE_CXX_STANDARD 17)
#set(CMAKE_CXX_FLAGS_RELEASE "${MAKE_CXX_FLAGS_RELEASE} -O3 -DNDEBUG")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")

#include_directories(extern MIMESim TestMIMESim)
#TODO ggf checken ob beim MIME + DCA überhaupt benötigt wird, dann hier rein schieben
include_directories(extern)
#TODO ggc noch Egien hier mit aufnehmen, sonst muss es irgendwo vorhanden sein
#including static library for sam2counts
#link_directories(../libs)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/609281088cfefc76f9d0ce82e1ff6c30cc3591e5.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)


#TODO nötig? Oder setzte ich das am Anfang im Projekt?
#set_property(TARGET MIMESim PROPERTY VERSION ${MIMESimVERSION})

# Some other library that we import that was also built using CMake
# and has an exported target.
#find_package(sam2counts REQUIRED)

enable_testing()

# The compiled library code is here
add_subdirectory(src)
add_subdirectory(test)



################## Testing #########################################################

if(UNIT_TESTING)
    MESSAGE(${CMAKE_CXX_COMPILER_ID})
    #gcov is only available with gcc
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND MIMESimTest_WITH_COV)
        #Flags to analyse the code coverage of the unit tests
        # --coverage =  generate additional info for gcov, is a synonym for -ftest-coverage, -fprofile-arcs and -lgcov
        #-g = debugging info, -p = collect profiling stats during runtime for code coverage analysis
        #-o0 = compiler allocates distinct stack memory for variables without register storage
        #--coverage = generate additional info for gcov
        #-ftest-coverage = generate info in '.gcno' to reconstruct the basic block graphs and assign source line numbers to blocks.
        #-fprofile-arcs = a separate ‘.gcda’ file is created for eachobject file compiled, Add code so that program flowarcsare instrumented
        #TODO manual of gcov says either .gcda or .gcno ??
        SET(GCC_COVERAGE_COMPILE_FLAGS "-O0 --coverage -fprofile-arcs")
       # SET(GCC_COVERAGE_LINK_FLAGS    " -lgcov")
        SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}" )
        MESSAGE(${CMAKE_CXX_FLAGS})
       # SET( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${GCC_COVERAGE_LINK_FLAGS}" )
    endif()
    #add_subdirectory(test)

endif()