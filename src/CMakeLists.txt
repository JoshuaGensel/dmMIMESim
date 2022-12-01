

#TODO kann man das schöner mache mit den Pfaden?
set(MIMESIM_HEADER_LIST
        ../include/MIMESim/Constants.hpp
        ../include/MIMESim/Mutation.hpp
        ../include/MIMESim/Species.hpp
        ../include/MIMESim/Utils.hpp
        ../include/MIMESim/FunctionalSequence.hpp
        ../include/MIMESim/BindingCompetition.hpp
        ../include/MIMESim/Count.hpp
        ../include/MIMESim/Generator.hpp)

set(MIMESIM_SOURCE_LIST
        Constants.cpp
        Mutation.cpp
        Species.cpp
        Utils.cpp
        FunctionalSequence.cpp
        BindingCompetition.cpp
        Count.cpp
        Generator.cpp)

#Set(BUILD_SHARED_LIBS ON)
#MESSAGE(HALLOOO ${BUILD_SHARED_LIBS})

#create a library for test purposes (a static one for now)
add_library(MIMESim ${MIMESIM_SOURCE_LIST} ${MIMESIM_HEADER_LIST})
# tell the library that the headers are given in this directory
#PUBLICH: use it for specified target and for targets which links with this project
target_include_directories(MIMESim PUBLIC ../include/MIMESim)

# shared libraries need PIC
#set_property(TARGET DCABenchmark PROPERTY POSITION_INDEPENDENT_CODE 1)
target_link_libraries(MIMESim Eigen3::Eigen)
#is linked automatically with gcc-9
target_link_libraries(MIMESim stdc++fs)
target_link_libraries(MIMESim
    gtest
)


add_executable(MIMESim_prog
        main.cpp
        #${DCABENCHMARK_SOURCE_LIST}
        #${DCABECHMARK_HEADER_LIST}
        )
# We need this directory, and users of our library will need it too
#target_include_directories(DCABenchmark_prog PUBLIC ../include)

target_link_libraries(MIMESim_prog MIMESim)


