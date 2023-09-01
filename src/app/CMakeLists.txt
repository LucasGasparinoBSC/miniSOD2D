project(miniSOD2D)

if(SHARED_LIBS STREQUAL "ON")
    file(GLOB_RECURSE SRC_FILES ${CMAKE_CURRENT_LIST_DIR}/*.f90)
else()
    file(GLOB_RECURSE SRC_FILES ${CMAKE_SOURCE_DIR}/src/utils/sources/*.f90
                                ${CMAKE_SOURCE_DIR}/src/geometry/sources/*.f90
                                ${CMAKE_SOURCE_DIR}/src/fem/sources/*.f90
                                ${CMAKE_SOURCE_DIR}/src/kernels/sources/*.f90
                                ${CMAKE_CURRENT_LIST_DIR}/*.f90)
endif()

add_definitions(-D__BUILD_VERSION__="${TEST_VERSION}")

add_executable(${PROJECT_NAME} ${SRC_FILES})

if(SHARED_LIBS)
    target_link_libraries(${PROJECT_NAME} tmu_utils tmu_geometry tmu_fem tmu_kernels)
endif()

set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME "miniSOD" VERSION ${TEST_VERSION} SOVERSION ${TEST_VERSION_MAJOR})

install(TARGETS ${PROJECT_NAME} DESTINATION bin)