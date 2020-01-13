add_executable(solver solver.cu)    
target_include_directories(solver PUBLIC ${SOLVER_SOURCE_DIR}/src/)
target_link_libraries(solver cublas cusparse)

add_executable(staggered_solver staggered_solver.cu)    
target_include_directories(staggered_solver PUBLIC ${SOLVER_SOURCE_DIR}/src/)
target_link_libraries(staggered_solver cublas cusparse)

add_executable(collocated_solver collocated_solver.cu)    
target_include_directories(collocated_solver PUBLIC ${SOLVER_SOURCE_DIR}/src/)
target_link_libraries(collocated_solver cublas cusparse)
