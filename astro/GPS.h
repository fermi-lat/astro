/** @file GPS.h
@brief declare class GPS

$Header: /nfs/slac/g/glast/ground/cvs/astro/astro/GPS.h,v 1.14 2006/10/06 01:16:36 burnett Exp $
*/
#ifndef ASTRO_GPS_H
#define ASTRO_GPS_H

#include "facilities/Scheduler.h"
#include "facilities/Observer.h"

#include "astro/SkyDir.h"
#include "astro/EarthCoordinate.h"
#include "astro/PointingInfo.h"

#include "CLHEP/Vector/Rotation.h"

#include <string>
namespace astro {

    // forward declarations
class PointingHistory;
class EarthOrbit;

/** 
* \class GPS
* \brief Models the Global Positoning System for a spacecraft. Handles time, position, and orientation.
* 
Represents the Global Positioning System on-board the spacecraft. An Orbit
object is used to compute the spacecraft's position and pointing characteristics.
Time is tracked through this object, and synchronized with the Scheduler for 
discrete event simulation. 

An expansion factor is provided to allow for acceleration. 

*/
class GPS  
{
public:

    enum CoordSystem { 
        GLAST=0,  //! 0: The original direction is in the GLAST frame already
        ZENITH=1, //! 1: rotate from the earth-zenith frame to GLAST
        CELESTIAL=2 //! 2: rotate from the cartesian celestial coordinate system (like a SkyDir)
    };

    enum RockType { 
        NONE,  ///<  No rocking rotation done at all.
        UPDOWN, ///< Satellite will be rocked toward the north pole in the northern hemisphere, opposite in the south.
        SLEWING, ///< (experimental) like UPDOWN, except that rotation at equator happens gradually.
        ONEPERORBIT, ///< LAT rocked northward for one orbit, southward for the next.
        EXPLICIT, ///<  Explicit angles given - this is used only if rotAngles get set.
        POINT, //!  Inertial pointing 
        HISTORY //! This setting is for using a previously generated pointing database to represent the orbit.
    };

    double	time () const; /// <current time

    double lat()const;  /// < latitude (degrees)
    double lon()const; ///  < longitude (degrees)
    double altitude()const; ///< altitude (km)

    astro::SkyDir zAxisDir()const; ///< spacecraft z-axis direction
    astro::SkyDir xAxisDir()const; ///< spacecraft x-axis direction
    astro::SkyDir zenithDir()const;///< local zenith direction
    CLHEP::Hep3Vector position()const;///< return current position;

    /// return a rotation matrix for the requested transformation
    CLHEP::HepRotation transformToGlast(double seconds,CoordSystem index);

    /// expansion of the current orbit
    double      expansion () const; 
    /// sample interval for random orbit distribution
    double     sampleintvl () const; 

    /// position in Earth coordinates
    astro::EarthCoordinate earthpos()const;

    /// pass a specific amount of time
    void    pass ( double );
    /// set the expansion factor for the orbit (-1) = random
    void    expansion ( double );
    /// synchronize w. scheduler
    void    synch ();    
    /// set the sample interval
    void    sampleintvl ( double );

    /// set the direction to point
    void setPointingDirection(const astro::SkyDir& dir);

    /** @brief set the desired pointing history file to use. 
        @param fileName 
        @param offset mission elapsed time for "launch", 
              number to be added to time increments

       The file can be either FT1 or an ascii format
    */
    void setPointingHistoryFile(std::string fileName, double offset=0);

    // notification support, managed by facilities/Observer
    void notifyObservers() { m_notification.notify();}
    Subject&    notification() { return m_notification; }

    // static access/destruction
    static GPS*	instance();
    static void     kill ();

    /// set angle to "rock"
    double rockingDegrees(double rockDegrees);

    /// set pointing strategy
    GPS::RockType setRockType(RockType rockType);

    void    time ( double );// set time

    double endTime()const{return m_endTime==0? 1e30: m_endTime;}

    static int test();

protected:
    // singleton - protect ctor/dtor
    GPS();
    virtual ~GPS();

private:
    static GPS* s_instance;
    astro::EarthOrbit* m_earthOrbit; //orbital position object, from the astro package.

    double  m_expansion;    // orbit expansion factor
    double m_time;	    // global time
    double m_lastQueriedTime; //the last time that was asked for

    // notification
    double  m_sampleintvl;  // interval to sample for each pt. in the orbit - to normalize spectra
    Subject  m_notification; 

    double m_rockDegrees; //number of degrees to "rock" the spacecraft, along the local x axis. 
    RockType m_rockType;//current rocking scheme
    
    astro::PointingHistory* m_history;
    astro::PointingInfo m_currentPoint;

    double m_endTime;

    astro::SkyDir m_point; ///< set for pointing 

    ///! update position, orientaion
    void update(double inputTime);

};
} // namespace
#endif 
