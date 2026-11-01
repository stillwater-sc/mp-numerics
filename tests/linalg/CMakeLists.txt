# Linear-algebra kernels across number systems (stillwater-sc/mtl5#666).
mpnumerics_add_test(test_exact)
mpnumerics_add_test(test_rotations)
mpnumerics_add_test(test_givens_qr)
mpnumerics_add_test(test_householder_qr)
# test_rotations compares float and double with cfloat bit for bit, and
# test_householder_qr rebuilds mtl::householder's v and tau bit for bit from
# householder_beta. A compiler that fuses a*b + c into one FMA (Clang on
# AArch64) would break both. MTL5 is header-only, so the setting covers its
# side of each comparison too.
if(NOT MSVC)
	target_compile_options(test_rotations PRIVATE -ffp-contract=off)
	target_compile_options(test_householder_qr PRIVATE -ffp-contract=off)
endif()
