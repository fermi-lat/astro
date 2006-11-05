/** @file PointingInfo.h
    @brief declare class PointingInfo

    $Header$
*/
#ifndef ASTRO_POINTINGINFO_H
#define ASTRO_POINTINGINFO_H

#include "astro/SkyDir.h"
#include "astro/EarthCoordinate.h"
#include "astro/Quaternion.h"
namespace astro{

class PointingInfo {
public:
    /** ctor sets everyting
        @param position position
        @param orientation
        @param earthPos
    */
    PointingInfo(const CLHEP::Hep3Vector& position, 
        const astro::Quaternion& orientation,
        const astro::EarthCoordinate& earthPos);


    /// default ctor
    PointingInfo(){}

    ~PointingInfo(){}

    // accessors
    astro::SkyDir xAxis()const;
    astro::SkyDir zAxis()const;
    astro::SkyDir zenith()const;
    const CLHEP::Hep3Vector& position()const{return m_position;}
    const astro::EarthCoordinate& earthCoord()const{return m_earth;}

    ///! equivalent rotation
    CLHEP::HepRotation rotation()const{return m_q.rotation();}

    PointingInfo interpolate(const astro::PointingInfo& next, double fraction)const;

private:
    CLHEP::Hep3Vector m_position; ///< position (km)
    astro::Quaternion m_q; ///< orientaion of the SC: xaxis direction
    astro::EarthCoordinate m_earth; ///< Earth position: (lat, lon)
};

}

#endif
