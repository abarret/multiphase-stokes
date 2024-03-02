SET(_dirs apply solve channel)
FOREACH(_dir ${_dirs})
	ADD_SUBDIRECTORY(${_dir})
ENDFOREACH()
