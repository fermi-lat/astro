/** @file PointingInfo.h
    @brief implement class PointingInfo

    $Header: /nfs/slac/g/glast/ground/cvs/astro/src/PointingInfo.cxx,v 1.2 2007/10/23 17:52:21 burnett Exp $
*/

#include "astro/PointingInfo.h"

using namespace astro;

PointingInfo::PointingInfo(const CLHEP::Hep3Vector& position, 
                           const Quaternion& orientation,
                           const EarthCoordinate& earthPos)
: m_position(position)
, m_q(orientation)
, m_earth(earthPos)
{}

astro::SkyDir PointingInfo::xAxis()const
{
    return m_q.rotate(Hep3Vector(1,0,0));
}

astro::SkyDir PointingInfo::zAxis()const
{
    return m_q.rotate(Hep3Vector(0,0,1));
}

astro::SkyDir PointingInfo::zenith()const
{
    return SkyDir(m_position);
}

namespace {
    template<typename T>
    inline T linear_interp(T a, T b, double f){
        return (1-f)*a + f*b;
    }
}
astro::PointingInfo PointingInfo::interpolate(const astro::PointingInfo& next, double f)const
{
    using CLHEP::Hep3Vector;

    // linear interpolation of earth location
    double lat1( earthCoord().latitude() )
         , lat2(next.earthCoord().latitude())
         , lat( linear_interp(lat1, lat2, f) );
    double lon1( earthCoord().longitude() )
         , lon2(next.earthCoord().longitude())
         , lon( linear_interp(lon1,lon2, f) );

    //this piece of code should just handle the "wraparound" cases:
    if(fabs(lon1-lon2) >= 330.){
        //we have gone off one end of the longitude scale.
        double lonlesser=std::max(lon1,lon2);
        double longreater=std::min(lon1,lon2)+360.;
        lon = lonlesser+((longreater-lonlesser)*f);
        while(lon > 360.)lon -= 360.;
    }

    // linear interpolation of position
    Hep3Vector pos1( position()), pos2(next.position());
    double alt1( pos1.mag() )
         , alt2( pos2.mag() )
         , alt( linear_interp(alt1,alt2,f) );
    Hep3Vector position (linear_interp(pos1,pos2,f).unit() * alt );

    // note using the quaternion interpolation (SLERP)
    return PointingInfo(position, m_q.interpolate(next.m_q, f), earthCoord());
 
    return *this; // todo

}

