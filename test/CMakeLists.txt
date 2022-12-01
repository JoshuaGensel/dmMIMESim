# include(GoogleTest)

set(TEST_SOURCE_FILES
  Constants_test
  FunctionalSequence_test
  Species_test
  main_test
)

add_executable(all_tests ${TEST_SOURCE_FILES})
target_link_libraries(all_tests
    ${PROJECT_NAME}
    gtest
    gtest_main
)
