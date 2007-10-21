// $Header: /nfs/slac/g/glast/ground/cvs/astro/src/mainpage.h,v 1.5 2006/11/05 20:06:27 burnett Exp $
// Mainpage for doxygen

/*! \mainpage package astro

   \authors  Toby Burnett, Sean Robinson, Theodore Hierath, Gino Tosti


  This package builds a library with useful Astronomical utility definitions

  - astro::JulianDate Represent a Julian Date
  - astro::SkyDir Point in the sky
  - astro::EarthCoordinate Point with respect to the surface of the Earth
  - astro::EarthOrbit Dynamics of an Earth orbit
  - astro::SolarSystem Sky positions of Moon and Sun
  - astro::PointingTransform Transformation between GLAST and Celestial coordinate systems
  - astro::SkyProj Image projections to/from celestial coordinates (wrapper of wcslib)
  - astro::Healpix manage transformation to/from Healpix pixels, as an STL container
  - astro::HealPixel Special class to support multiple pixel sizes 
  - astro::Quaternion Quaternion objects, used to represent rotations
  - astro::IGRField interface to the fortran code for the magnetic field. Used by EarthCoordinate

    <hr>
  \section notes release notes
  release.notes
  \section requirements requirements
  \include requirements

*/

