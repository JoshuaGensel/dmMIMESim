include(GoogleTest)

function(AddTest test)
  add_executable(
    ${test}
    ${test}.cpp
    )
  target_link_libraries(
    ${test}
    gtest_main
    MIMESim
    )
    
    gtest_discover_tests(${test} WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/test)
endfunction()
  
AddTest("Constants_test")
AddTest("Species_test")
