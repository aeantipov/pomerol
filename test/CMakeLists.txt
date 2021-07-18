# Build Catch2 main object file
add_library(catch2 OBJECT catch2/catch2-main.cpp)
target_include_directories(catch2 PRIVATE ${MPI_CXX_INCLUDE_PATH})

# Pomerol tests
set (tests
HamiltonianTest
GF1siteTest
GF2siteTest
GF4siteTest
AndersonTest
AndersonComplexTest
Anderson2PGFTest
Vertex4Test
SusceptibilityTest
)

include_directories(${PROJECT_SOURCE_DIR}/include ${PROJECT_BINARY_DIR}/include)

foreach (test ${tests})
    set(test_src ${test}.cpp)
    add_executable(${test} ${test_src})
    set(test_parameters
        ${MPIEXEC_NUMPROC_FLAG} 1
        ${MPIEXEC_PREFLAGS} ${test} ${MPIEXEC_POSTFLAGS})
    add_test(NAME ${test} COMMAND "${MPIEXEC}" ${test_parameters})
    target_link_libraries(${test} ${PROJECT_NAME} ${MPI_CXX_LIBRARIES} catch2)
endforeach(test)

set(mpi_tests BroadcastTest MPIDispatcherTest)
foreach (test ${mpi_tests})
    set(test_src ${test}.cpp)
    add_executable(${test} ${test_src})
    target_link_libraries(${test} PRIVATE
                          ${PROJECT_NAME} ${MPI_CXX_LIBRARIES} catch2
    )

    foreach (np 2 4 8 16)
        set(test_parameters
            ${MPIEXEC_NUMPROC_FLAG} ${np}
            ${MPIEXEC_PREFLAGS} ${test} ${MPIEXEC_POSTFLAGS})
        add_test(NAME ${test}${np}cpu COMMAND "${MPIEXEC}" ${test_parameters})
    endforeach(np)
endforeach(test)
