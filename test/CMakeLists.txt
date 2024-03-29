include(ExternalProject)

set(TARGET hello-lib)
set(GTEST_TARGET ${TARGET}_test)
file(GLOB_RECURSE TESTS_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)

set_directory_properties(PROPERTIES EP_PREFIX ${CMAKE_BINARY_DIR}/ThirdParty)

ExternalProject_Add(
    googletest
    SVN_REPOSITORY http://googletest.googlecode.com/svn/trunk/
    SVN_REVISION -r 660
    TIMEOUT 10
    CMAKE_ARGS -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
               -Dgtest_force_shared_crt=ON
    INSTALL_COMMAND ""
    UPDATE_COMMAND ""
    LOG_DOWNLOAD ON
    LOG_CONFIGURE ON
    LOG_BUILD ON)

ExternalProject_Get_Property(googletest source_dir)
include_directories(${source_dir}/include)

set_target_properties(googletest PROPERTIES EXCLUDE_FROM_ALL TRUE)
ExternalProject_Get_Property(googletest binary_dir)
if(MSVC)
  set(Suffix ".lib")
else()
  set(Suffix ".a")
  set(Pthread "-pthread")
endif()

add_executable(${GTEST_TARGET} EXCLUDE_FROM_ALL ${TESTS_FILES})
target_link_libraries(
    ${GTEST_TARGET}
    ${TARGET}
    ${binary_dir}/${CMAKE_FIND_LIBRARY_PREFIXES}gtest${Suffix}
    ${binary_dir}/${CMAKE_FIND_LIBRARY_PREFIXES}gtest_main${Suffix}
    ${LIBRARIES}
)

add_dependencies(${GTEST_TARGET} googletest)

if(CMAKE_COMPILER_IS_GNUC OR CMAKE_COMPILER_IS_GNUCXX)
    target_link_libraries(${GTEST_TARGET} pthread)
endif()

add_custom_target(check COMMAND ${EXECUTABLE_OUTPUT_PATH}/${GTEST_TARGET} --gtest_output="xml:${CMAKE_BINARY_DIR}/tests.xml" --gtest_color=yes || true)
add_dependencies(check ${GTEST_TARGET})

if(COVERAGE_ENABLED)
    find_package(gcovr REQUIRED)
    target_link_libraries(${GTEST_TARGET} -lgcov)
    add_custom_target(coverage WORKING_DIRECTORY ${CMAKE_BINARY_DIR} COMMAND ${GCOVR_EXECUTABLE} -r ${CMAKE_SOURCE_PATH}/src --xml -o ${CMAKE_BINARY_DIR}/coverage.xml)
    add_dependencies(coverage check)
endif()

# vim:ft=cmake

