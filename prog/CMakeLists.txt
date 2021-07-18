set(DEFAULT_PROGS
    anderson
    hubbard2d
    quantum_model
)

# Find or fetch gftools
find_package(gftools CONFIG)
if(gftools_FOUND)
    get_target_property(GFTOOLS_INCLUDE_DIR gftools INTERFACE_INCLUDE_DIRECTORIES)
else()
    message(STATUS "Fetching gftools")
    find_package(Git REQUIRED)
    execute_process(COMMAND "${GIT_EXECUTABLE}" "clone" "--branch" "master"
                    "https://github.com/aeantipov/gftools.git" "gftools"
                    WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
    set(GFTOOLS_INCLUDE_DIR "${CMAKE_BINARY_DIR}/gftools")
endif()
message(STATUS "gftools includes : ${GFTOOLS_INCLUDE_DIR}")

message(STATUS "Building executables")
if(NOT progs_list)
    set(progs_list ${DEFAULT_PROGS})
endif()

set(PROG_SOURCES
    main.cpp
    quantum_model.cpp
    anderson_model.hpp
    hubbard2d_model.hpp
)

foreach(model ${progs_list})
    list (FIND DEFAULT_PROGS ${model} _index)
    if (${_index} GREATER -1)
        string(TOUPPER ${model} def)
        set(def "POMEROL_${def}")
        set(prog_name "${model}.${PROJECT_NAME}")
        add_executable(${prog_name} ${PROG_SOURCES})
        target_include_directories(${prog_name}
                                   PRIVATE ${PROJECT_SOURCE_DIR}/include
                                           ${PROJECT_BINARY_DIR}/include
                                           ${GFTOOLS_INCLUDE_DIR}
        )
        target_compile_definitions(${prog_name} PRIVATE ${def})
        target_link_libraries(${prog_name} ${PROJECT_NAME})
        set_target_properties(${prog_name} PROPERTIES
                              RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
        install(TARGETS ${prog_name} DESTINATION bin)
    else()
        message(WARNING "Unknown program \"${model}\" has been skipped.")
    endif()
endforeach()
