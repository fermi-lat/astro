// $Header: /nfs/slac/g/glast/ground/cvs/astro/src/PointingTransform.cxx,v 1.2 2003/09/30 21:03:51 srobinsn Exp $

// Include files
#include "astro/PointingTransform.h"

namespace astro {

    /** @brief initialize from Z and X directions
    @param inputType SkyDir
    */
	PointingTransform::PointingTransform(SkyDir zDir, SkyDir xDir)
		:m_xDir(xDir),m_zDir(zDir)
	{}
        
	/** @brief transform a local to a celestial direction
    */
	HepRotation PointingTransform::localToCelestial () const{
		const Hep3Vector& xd(m_xDir.dir());
		const Hep3Vector& zd(m_zDir.dir());
		const Hep3Vector& yd( (m_zDir.dir()).cross( m_xDir.dir() ) );
		HepRotation ret(xd,yd,zd);
		return ret;
	}
		
	/** @brief get the celestial direction from the local one
    */
	SkyDir PointingTransform::gDir(Hep3Vector localDir) const{
		Hep3Vector dir(localDir);
		SkyDir ret(localToCelestial()*dir,SkyDir::CELESTIAL);
		return ret;
	}

} // namespace astro