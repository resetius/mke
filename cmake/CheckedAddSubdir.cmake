
macro (add_subdirectory_if_exists _subdir)
	if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${_subdir})
		add_subdirectory (${_subdir})
	endif ()
endmacro(add_subdirectory_if_exists)

