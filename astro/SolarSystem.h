// 
#ifndef astro_SolarSystem_H
#define astro_SolarSystem_H


#include "astro/JulianDate.h"
#include "astro/SkyDir.h"

// forward declaration
class SolSystem;

namespace astro {

/** @class SolarSystem
  * @brief position of various solarsystem entities
  * @author G. Tosti original code 
  * @author T. Burnett convert to class
  * $Id: SolarSystem.h,v 1.1.1.1 2002/08/13 00:20:45 burnett Exp $
  *
  * This is a thin wrapper to code in the class SolSystem
  */

class SolarSystem  {

public:
    /**
     * index of body that corresponds to the #define's in SolSystem  
     */
    enum Body{ 
        //VENUS=1, MARS=2, JUPITER=3,SATURN=4, 
        //URANUS=5,NEPTUNE=6, PLUTO=7, 
            Sun=8, Moon=9 };
    SolarSystem();

    /** @brief simple constructor that initializes
     * @param body Body index
     * @param date Julian date for given position 
     */
    SolarSystem(  Body body, JulianDate date );
    ~SolarSystem();

    /**
    * @brief set the body and date; retun the SkyDir 
    */
    SkyDir direction(Body body, JulianDate date) ;

    /**
    * @brief conversion operator that returns the SkyDir 
    */
    operator SkyDir()const { return m_dir; }
private:

    SkyDir m_dir;
    SolSystem* m_ss;
};

} // namespace astro

#endif
