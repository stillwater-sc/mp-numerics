# mp-numerics#13: Givens rotations and Householder reflectors, one at a time.
add_executable(givens_householder_study
	main.cpp study_ieee.cpp study_posit.cpp study_fixpnt.cpp
	study_cfloat_narrow.cpp study_cfloat_mid.cpp study_cfloat_wide.cpp
	study_lns.cpp study_takum.cpp study_dd_qd.cpp)
target_link_libraries(givens_householder_study PRIVATE sw::mp_numerics)
target_include_directories(givens_householder_study PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/..)
mpn_study_provenance(givens_householder_study)
