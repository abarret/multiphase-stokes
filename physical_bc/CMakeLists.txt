SET(_dirs apply solve)
FOREACH(_dir ${_dirs})
	ADD_SUBDIRECTORY(${_dir})
ENDFOREACH()
