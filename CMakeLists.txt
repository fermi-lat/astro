####################################
# Create target and set properties #
####################################

add_library(
  astro SHARED
  src/EarthCoordinate.cxx
  src/EarthOrbit.cxx
  src/Geomag.cxx
  src/GPS.cxx
  src/HealpixProj.cxx
  src/HTM.cxx
  src/IGRF_data.cxx
  src/IGRField.cxx
  src/JulianDate.cxx
  src/PointingHistory.cxx
  src/PointingInfo.cxx
  src/PointingTransform.cxx
  src/Quaternion.cxx
  src/SkyDir.cxx
  src/SkyProj.cxx
  src/SkyStat.cxx
  src/SolarSystem.cxx
  src/igrf_sub/igrf_sub.cxx
  src/jplephem/bary.cxx
  src/jplephem/dpleph.cxx
)

target_include_directories(
  astro PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:>
)

target_link_libraries(
  astro

  PUBLIC
  facilities
  healpix_cxx::healpix_cxx
  CLHEP::CLHEP

  PRIVATE
  tip
  cfitsio::cfitsio
  WCSlib::WCSlib
)

add_executable(test_astro src/test/test.cxx)
target_link_libraries(test_astro PRIVATE astro tip)


###############################################################
# Installation
###############################################################
install(DIRECTORY astro DESTINATION ${FERMI_INSTALL_INCLUDEDIR})
install(DIRECTORY data/ DESTINATION ${FERMI_INSTALL_REFDATADIR}/astro)

install(
  TARGETS astro test_astro
  # EXPORT fermiTargets
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
