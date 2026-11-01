# mp-numerics#14: QR and least squares across number systems (part 1).
add_executable(decompositions_study
	main.cpp study_ieee.cpp study_posit.cpp study_fixpnt.cpp)
target_link_libraries(decompositions_study PRIVATE sw::mp_numerics)
target_include_directories(decompositions_study PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/..)
mpn_study_provenance(decompositions_study)
