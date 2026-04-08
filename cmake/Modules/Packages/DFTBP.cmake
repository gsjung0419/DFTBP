set(DFTBPLUS_PREFIX "" CACHE PATH "Prefix of the DFTB+ installation used by the DFTBP package")
set(DFTBPLUS_LIBRARY "" CACHE FILEPATH "Path to the DFTB+ library")
set(DFTBPLUS_EXTRA_LIBS "" CACHE STRING "Additional link flags required by DFTB+")

if(NOT DFTBPLUS_LIBRARY)
  unset(DFTBPLUS_LIBRARY CACHE)
  find_library(DFTBPLUS_LIBRARY NAMES dftbplus
    HINTS
      ${DFTBPLUS_PREFIX}/lib
      ${LAMMPS_LIB_SOURCE_DIR}/dftbp/liblink
  )
endif()

if(NOT DFTBPLUS_LIBRARY)
  message(FATAL_ERROR
    "Could not find libdftbplus for the DFTBP package. "
    "Set DFTBPLUS_PREFIX or DFTBPLUS_LIBRARY.")
endif()

set(_dftbplus_link_items)

if(DFTBPLUS_EXTRA_LIBS)
  separate_arguments(_dftbplus_link_items NATIVE_COMMAND "${DFTBPLUS_EXTRA_LIBS}")
elseif(DFTBPLUS_PREFIX AND EXISTS "${DFTBPLUS_PREFIX}/lib/pkgconfig/dftbplus.pc")
  file(STRINGS "${DFTBPLUS_PREFIX}/lib/pkgconfig/dftbplus.pc" _pc_libs REGEX "^Libs:")
  file(STRINGS "${DFTBPLUS_PREFIX}/lib/pkgconfig/dftbplus.pc" _pc_libs_private REGEX "^Libs.private:")
  string(REGEX REPLACE "^Libs: *" "" _pc_libs "${_pc_libs}")
  string(REGEX REPLACE "^Libs.private: *" "" _pc_libs_private "${_pc_libs_private}")
  set(_dftbplus_raw_items "${_pc_libs} ${_pc_libs_private}")
  separate_arguments(_dftbplus_link_items NATIVE_COMMAND "${_dftbplus_raw_items}")
  list(REMOVE_ITEM _dftbplus_link_items "-ldftbplus")
endif()

target_link_libraries(lammps PRIVATE ${DFTBPLUS_LIBRARY} ${_dftbplus_link_items})
