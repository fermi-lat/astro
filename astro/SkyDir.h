// $Header: /nfs/slac/g/glast/ground/cvs/Event/Event/Utilities/SkyDir.h,v 1.3 2002/08/12 20:20:39 srobinsn Exp $
#ifndef OrbitModel_SkyDir_H
#define OrbitModel_SkyDir_H


// Include files
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

namespace astro {


/** @class SkyDir
  * @brief Describe an absolute direction
  * @author S. Robinson 
  * <br>$Id:$
  *
  */
class SkyDir
{
public:
    enum CoordSystem { 
        GALACTIC,  //!  fixed direction with respect to the galactic coordinate system (l,b)
        CELESTIAL //! fixed direction with respect to the celestial coordinate system (ra,dec) in the J2000 epoch.
    };
    ///Constructors
    ///(l,b) or (Ra, Dec) instantiation
    SkyDir(double param1=0, double param2=0, CoordSystem inputType = CELESTIAL);
    SkyDir(Hep3Vector);
    
    ///return methods
    Hep3Vector& operator () () {return m_dir;}
    double l () const;
    double b () const;
    double ra () const;
    double dec () const;
    const Hep3Vector& dir () const {return m_dir;}
    
    //!to return the opening angle (in radians) between two objects:
    double SkyDir::difference(const SkyDir& other)const;

    static const HepRotation& celToGal();
private:
    static HepRotation s_celestialToGalactic;

    Hep3Vector m_dir;
    std::pair<double,double> setGalCoordsFromDir() const;

};

inline double SkyDir::difference(const SkyDir& other)const{
    return 2.*asin(0.5*(m_dir-other.dir()).mag());
}

} // namespace astro
#endif    // LHCBEVENT_SKYDIR_H
