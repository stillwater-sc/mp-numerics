# Linear-algebra kernels across number systems (stillwater-sc/mtl5#666).

# Provenance for the studies' .sysinfo sidecars: the commits of this repo and
# of the MTL5 and Universal sources a study was built against, the build type
# and the flags. Taken at configure time: reconfigure before a run whose
# results are kept. A "+dirty" suffix marks uncommitted changes in that tree.
function(_mpn_commit out dir)
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

function(mpn_study_provenance target)
	_mpn_commit(mpn "${PROJECT_SOURCE_DIR}")
	_mpn_commit(mtl5 "${mtl5_SOURCE_DIR}")
	_mpn_commit(universal "${universal_SOURCE_DIR}")
	string(TOUPPER "${CMAKE_BUILD_TYPE}" bt)
	target_compile_definitions(${target} PRIVATE
		MPN_GIT_COMMIT="${mpn}"
		MTL5_GIT_COMMIT="${mtl5}"
		UNIVERSAL_GIT_COMMIT="${universal}"
		MPN_BUILD_TYPE="${CMAKE_BUILD_TYPE}"
		MPN_CXX_FLAGS="${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${bt}}")
endfunction()

add_subdirectory(givens_householder)
add_subdirectory(decompositions)
