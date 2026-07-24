# Applications (demonstration programs)
# Organized by numerical-method category (roots, integration, interpolation,
# optimization, approximation); each leaf subdirectory has its own CMakeLists.

# Composite Simpson quadrature accuracy across number precisions.
add_subdirectory(integration/simpson_precision)

# Scalar root finding (Newton / secant) accuracy across number precisions.
add_subdirectory(roots/root_precision)

# Newton interpolation accuracy across number precisions.
add_subdirectory(interpolation/interp_precision)
# 1-D golden-section minimization accuracy across number precisions.
add_subdirectory(optimization/optim_precision)

# Chebyshev-series approximation accuracy across number precisions.
add_subdirectory(approximation/cheb_precision)
