include(CTest)

# Helper: each test is a self-checking executable linked against sw::mp_numerics.
function(mpnumerics_add_test test_name)
	add_executable(${test_name} ${test_name}.cpp)
	target_link_libraries(${test_name} PRIVATE sw::mp_numerics)
	add_test(NAME ${test_name} COMMAND ${test_name})
endfunction()

# Tests are organized by numerical-method category:
#   integration   -- quadrature (Simpson, trapezoidal, ...)
#   roots         -- polynomial / nonlinear root finding
#   interpolation -- polynomial / spline interpolation
#   optimization  -- line search / descent
#   approximation -- Chebyshev / minimax approximation
add_subdirectory(integration)
add_subdirectory(roots)
add_subdirectory(interpolation)
add_subdirectory(optimization)
add_subdirectory(approximation)
