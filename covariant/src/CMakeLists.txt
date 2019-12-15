add_executable(solver solver.cu)    
target_include_directories(solver PUBLIC ${SOLVER_SOURCE_DIR}/src/)
target_link_libraries(solver cublas cusparse)

add_executable(stencil_solver stencil_solver.cu)    
target_include_directories(stencil_solver PUBLIC ${SOLVER_SOURCE_DIR}/src/)
target_link_libraries(stencil_solver cublas cusparse)
