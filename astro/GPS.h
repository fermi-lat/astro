// $Header: /nfs/slac/g/glast/ground/cvs/astro/astro/GPS.h,v 1.4 2004/11/12 01:52:23 burnett Exp $

#if !defined(_H_GPS_CLASS)
#define _H_GPS_CLASS

#include "facilities/Scheduler.h"
#include "facilities/Observer.h"

#include "astro/SkyDir.h"
#include "astro/EarthOrbit.h"
#include "CLHEP/Vector/Rotation.h"

#include <iostream>
#include <string>


/** 
* \class GPS
* \brief Models the Global Positoning System for a spacecraft. Handles time, position, and orientation for the instrument as a whole.
* 
* $Header: /nfs/slac/g/glast/ground/cvs/astro/astro/GPS.h,v 1.4 2004/11/12 01:52:23 burnett Exp $
Represents the Global Positioning System on-board the spacecraft. An Orbit
object is used to compute the spacecraft's position and pointing characteristics.
Time is tracked through this object, and synchronized with the Scheduler for 
discrete event simulation. An expansion factor is provided to allow for acceleration
or randomization of the orbit. If the expansion factor is negative, then the position
of the spacecraft is chosen as a random distribution over an orbit. Otherwise, the expansion
factor represents an acceleration of the spacecraft's orbit. Ie. an expansion factor
of 2 would reduce the orbit period of the spacecraft by 1/2.

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
        NONE,  //!  No rocking rotation done at all.
        UPDOWN, //! Satellite will be rocked toward the north pole in the northern hemisphere, opposite in the south.
        SLEWING, //! (experimental) like UPDOWN, except that rotation at equator happens gradually.
        ONEPERORBIT, //! LAT rocked northward for one orbit, southward for the next.
        EXPLICIT, //!  Explicit angles given - this is used only if rotAngles get set.
        POINT, //!  The Lat points in a given direction.  Here, m_rotangles is in (l,b) format for pointing.
        HISTORY //! This setting is for using a previously generated pointing database to represent the orbit.
    };

    typedef struct{
        astro::SkyDir dirZ;
        astro::SkyDir dirX;
        double lat,lon;
        Hep3Vector position;
        double altitude;
    }POINTINFO;

    typedef std::map<double,GPS::POINTINFO> HistoryMap;
    typedef HistoryMap::const_iterator history_iterator;
    const HistoryMap& getHistory()const{return m_pointingHistory;}

    class Coords {
    public:
        Coords( double alat, double alon, double apitch
            , double ayaw, double aroll, GPStime atime, double phase ); 
        Coords();

        GPStime time () const { return m_time; }
        double  lat () const { return m_lat; }
        double  lon () const { return m_lon; }
        double  phase () const { return m_phase; }

    private:
        GPStime m_time;
        double m_lat, m_lon, m_pitch, m_yaw, m_roll, m_phase;
        double m_altitude;
    };

    // const access

    /// GPS synchronized time for the satellite
    GPStime	time () const; 
    /// present latitude
    double lat()const;//{getPointingCharacteristics(time);return m_lat;} 
    /// present longitude
    double lon()const;//{getPointingCharacteristics(time);return m_lon;}  
    double altitude()const; // rad_geo
    /// pointing characteristics
    double RAX()const;//{getPointingCharacteristics(time);return m_RAX;}
    double RAZ()const;//{getPointingCharacteristics(time);return m_RAZ;}
    double DECX()const;//{getPointingCharacteristics(time);return m_DECX;}
    double DECZ()const;//{getPointingCharacteristics(time);return m_DECZ;}
    double RAZenith()const;//{getPointingCharacteristics(time);return m_RAZenith;}
    double DECZenith()const;//{getPointingCharacteristics(time);return m_DECZenith;}
    /// expansion of the current orbit
    double      expansion () const; 
    /// sample interval for random orbit distribution
    GPStime     sampleintvl () const; 
    /// return the orbit's ascending longitude 
    double     ascendingLon()const; 
    /// access m_rotangles
    std::pair<double,double> rotateAngles(); 


    // set data

    /// get the pointing characteristics of the satellite, given a location and rocking angle.
    void getPointingCharacteristics(double inputTime);

    /// pass a specific amount of time
    void    pass ( double );
    /// set the expansion factor for the orbit (-1) = random
    void    expansion ( double );
    /// synchronize w. scheduler
    void    synch ();    
    /// set the sample interval
    void    sampleintvl ( GPStime );
    /// special to set the ascending longitude	
    void    ascendingLon(double);   
    /// set m_rotangles
    void    rotateAngles(std::pair<double,double> coords); 

    /// set the desired pointing history file to use:
    void setPointingHistoryFile(std::string fileName);

    /// write the explicit history data for re-creation of orbit.
    void setUpHistory();

    /// print time & position
    void    printOn(std::ostream& out) const; 

    // notification
    void notifyObservers() { m_notification.notify();}
    Subject&    notification() { return m_notification; }

    // static access/destruction
    static GPS*	instance();
    static void     kill ();

    //this is the only external rotation function that should survive.
    //all the others are being phased out outside of GPS, whiule this one should
    //take care of the various necessary rotations.
    HepRotation transformToGlast(double seconds,CoordSystem index);

    /// return the rotation for compensation for the rocking angles.
    HepRotation rockingAngleTransform(double seconds);

    ///this transforms glast-local (cartesian) vectors into galactic (cartesian) vectors
    HepRotation transformGlastToGalactic(double seconds);

    HepRotation CELTransform(double seconds);

    HepRotation transformCelToGlast(double seconds);

    double rockingDegrees(double rockDegrees){double ret=m_rockDegrees;
    m_rockDegrees = rockDegrees;
    return ret;}

    int setRockType(RockType rockType);
    int setRockType(int rockType);

    void    time ( GPStime );// set time

    Hep3Vector position(double seconds)const{
        if(m_rockType == HISTORY){
            instance()->setInterpPoint(seconds);
            return m_currentInterpPoint.position;
        }
        double time = m_earthOrbit->dateFromSeconds(seconds);
        return m_earthOrbit->position(time);
        /*return m_position;*/} //interface to EarthOrbit::position()

protected:
    // singleton - protect ctor/dtor
    GPS();
    virtual ~GPS();

    std::pair<double,double> m_rotangles;  //angles for coordinate rotation (rocking angle)

    // friends
    friend class FluxGenerator;

private:
    static GPS* s_instance;
    astro::EarthOrbit* m_earthOrbit; //orbital position object, from the astro package.

    void setInterpPoint(double time);

    double  m_expansion;    // orbit expansion factor
    GPStime m_time;	    // global time
    double m_lastQueriedTime; //the last time that was asked for
    double  m_sampleintvl;  // interval to sample for each pt. in the orbit - to normalize spectra
    double m_lat,m_lon; //position characteristics
    double m_altitude; 
    double m_RAX,m_RAZ,m_DECX,m_DECZ; //pointing characteristics.
    double m_RAZenith,m_DECZenith,m_RAXZenith,m_DECXZenith; //pointing characteristic of the zenith direction.
    Hep3Vector m_position; //current vector position of the LAT.
    // notification
    Subject    m_notification; 
    double m_rockDegrees; //number of degrees to "rock" the spacecraft, along the local x axis. 
    RockType m_rockType;//current rocking scheme
    double m_rockNorth; //internal value for the current number of degrees the craft is rotated at the time.
    std::string m_pointingHistoryFile;//pointing/livetime database history file to use.
    HistoryMap m_pointingHistory;//pointing/livetime database history
    POINTINFO m_currentInterpPoint; //holder object for currently interpotated pointing information

   bool haveFitsFile() const;
   void readFitsData();
   void fitsReportError(FILE *, int) const;
};

inline std::istream&    operator>>(std::istream& i, GPS::Coords& c) {
    double  p, y, r, lt, ln, tm, ph;
    i >> lt; i >> ln; i >> p; i >> y; i >> r; i >> tm; i >> ph;
    c = GPS::Coords( lt, ln, p, y, r, tm, ph );
    return i;
}

inline std::ostream&    operator<<(std::ostream& o, const GPS::Coords& c) {
    o << ' ' << c.lat() << ' ' << c.lon() << ' '  
        << c.time() <<' ' << c.phase();
    return o;
}
#endif // !defined(AFX_GPS_H__F9844433_4E64_11D2_B4DD_00A0C9960210__INCLUDED_)
