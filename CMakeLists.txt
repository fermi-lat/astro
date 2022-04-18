####################################
# Create target and set properties #
####################################

add_library(astro SHARED
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

target_include_directories(astro PUBLIC
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           $<INSTALL_INTERFACE:>
)

target_link_libraries(astro
    PUBLIC
    facilities
    tip
    CLHEP::CLHEP
    OpenMP::OpenMP_CXX
    PkgConfig::cfitsio
    PkgConfig::wcslib
    PkgConfig::healpix_cxx
)

add_executable(test_astro src/test/test.cxx)
target_link_libraries(test_astro astro)
