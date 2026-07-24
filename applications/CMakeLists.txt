# Applications (demonstration programs)
# Organized by numerical-method category (roots, integration, interpolation,
# optimization, approximation); each leaf subdirectory has its own CMakeLists.

# Composite Simpson quadrature accuracy across number precisions.
add_subdirectory(integration/simpson_precision)

# Scalar root finding (Newton / secant) accuracy across number precisions.
add_subdirectory(roots/root_precision)
