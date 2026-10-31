# mp-numerics#13: Givens rotations and Householder reflectors, one at a time.
add_executable(givens_householder_study
	main.cpp study_ieee.cpp study_posit.cpp study_fixpnt.cpp
	study_cfloat_narrow.cpp study_cfloat_mid.cpp study_cfloat_wide.cpp
	study_lns.cpp study_takum.cpp study_dd_qd.cpp)
target_link_libraries(givens_householder_study PRIVATE sw::mp_numerics)

# Provenance for the .sysinfo sidecars: the commits of this repo and of the
# MTL5 and Universal sources the study was built against.
# Taken at configure time: reconfigure before a run whose results are kept.
# A "+dirty" suffix marks uncommitted changes in that tree.
function(_gh_commit out dir)
	set(sha "unknown")
	if(dir AND EXISTS "${dir}")
		execute_process(COMMAND git -C "${dir}" rev-parse --short=12 HEAD
			OUTPUT_VARIABLE sha OUTPUT_STRIP_TRAILING_WHITESPACE
			RESULT_VARIABLE rc ERROR_QUIET)
		if(NOT rc EQUAL 0)
			set(sha "unknown")
		else()
			execute_process(COMMAND git -C "${dir}" status --porcelain --untracked-files=no
				OUTPUT_VARIABLE st OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
			if(st)
				set(sha "${sha}+dirty")
			endif()
		endif()
	endif()
	set(${out} "${sha}" PARENT_SCOPE)
endfunction()
_gh_commit(_gh_mpn "${PROJECT_SOURCE_DIR}")
_gh_commit(_gh_mtl5 "${mtl5_SOURCE_DIR}")
_gh_commit(_gh_universal "${universal_SOURCE_DIR}")
string(TOUPPER "${CMAKE_BUILD_TYPE}" _gh_bt)
target_compile_definitions(givens_householder_study PRIVATE
	MPN_GIT_COMMIT="${_gh_mpn}"
	MTL5_GIT_COMMIT="${_gh_mtl5}"
	UNIVERSAL_GIT_COMMIT="${_gh_universal}"
	MPN_BUILD_TYPE="${CMAKE_BUILD_TYPE}"
	MPN_CXX_FLAGS="${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${_gh_bt}}")
