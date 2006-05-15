/** @file Photon.h
@brief Define the Photon class

$Header: /nfs/slac/g/glast/ground/cvs/astro/astro/Photon.h,v 1.1 2006/03/16 04:39:02 burnett Exp $
*/

#ifndef astro_Photon_h
#define astro_Photon_h

#include "astro/SkyDir.h"


namespace astro {

    /** @class Photon
    @brief add energy and event class attributes to a SkyDir.


    */
    class Photon : public astro::SkyDir {
    public:
        Photon(const astro::SkyDir& dir, double energy, double time, int event_class=-1)
            : astro::SkyDir(dir)
            , m_energy(energy)
            , m_time(time)
            , m_event_class(event_class)
        {}
        double energy()const{return m_energy;}
        int eventClass()const{return m_event_class;}
        double time()const{return m_time;}

    private:
        double m_energy;
        double m_time;
        int m_event_class;

    };

}// namespace astro

#endif
