# flex a .ll file

INCLUDE(FindBison)


MACRO(ADD_BISON_FILES _sources )
    FIND_PACKAGE(Bison QUIET REQUIRED)

    FOREACH (_current_FILE ${ARGN})
      GET_FILENAME_COMPONENT(_in ${_current_FILE} ABSOLUTE)
      GET_FILENAME_COMPONENT(_basename ${_current_FILE} NAME_WE)

      SET(_out ${CMAKE_CURRENT_SOURCE_DIR}/${_basename})

      ADD_CUSTOM_COMMAND(
         OUTPUT ${_out}.cpp ${_out}.hpp
         COMMAND ${BISON_EXECUTABLE}
         ARGS
	 -d
         -o${_out}.cpp
         ${_in}
         DEPENDS ${_in}
      )

      SET(${_sources} ${${_sources}} ${_out} )
   ENDFOREACH (_current_FILE)
ENDMACRO(ADD_BISON_FILES)
