# Linear-algebra kernels across number systems (stillwater-sc/mtl5#666).
mpnumerics_add_test(test_exact)
mpnumerics_add_test(test_rotations)
# test_rotations compares float and double with cfloat bit for bit, which a
# compiler that fuses a*b + c into one FMA (Clang on AArch64) would break.
if(NOT MSVC)
	target_compile_options(test_rotations PRIVATE -ffp-contract=off)
endif()
