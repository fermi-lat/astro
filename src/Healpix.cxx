/** @file Healpix.cxx
    @brief Healpix class implementation with code from WMAP

    @author B. Lesnick 
    $Header: /nfs/slac/g/glast/ground/cvs/astro/src/Healpix.cxx,v 1.3 2005/03/29 19:13:57 burnett Exp $
*/
/* Local Includes */
#include "astro/Healpix.h"

/* Standard Includes */
#include <numeric> // for accumulate

#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cassert>


using namespace astro;


namespace {
#define NS_MAX 8172  // Max allowed value for nside.

void mk_xy2pix(int *x2pix, int *y2pix) {
  /* =======================================================================
   * subroutine mk_xy2pix
   * =======================================================================
   * sets the array giving the number of the pixel lying in (x,y)
   * x and y are in {1,128}
   * the pixel number is in {0,128**2-1}
   *
   * if  i-1 = sum_p=0  b_p * 2^p
   * then ix = sum_p=0  b_p * 4^p
   * iy = 2*ix
   * ix + iy in {0, 128**2 -1}
   * =======================================================================
   */
  int i, K,IP,I,J,ID;
  
  for(i = 0; i < 127; i++) x2pix[i] = 0;
  for( I=1;I<=128;I++ ) {
    J  = I-1;//            !pixel numbers
    K  = 0;//
    IP = 1;//
    truc : if( J==0 ) {
      x2pix[I-1] = K;
      y2pix[I-1] = 2*K;
    }
    else {
      ID = J % 2;
      J  = J/2;
      K  = IP*ID+K;
      IP = IP*4;
      goto truc;
    }
  }     
  
}


void ang2pix_nest( const long nside, double theta, double phi, long *ipix) {

  /* =======================================================================
   * subroutine ang2pix_nest(nside, theta, phi, ipix)
   * =======================================================================
   * gives the pixel number ipix (NESTED) corresponding to angles theta and phi
   *
   * the computation is made to the highest resolution available (nside=8192)
   * and then degraded to that required (by integer division)
   * this doesn't cost more, and it makes sure that the treatement of round-off 
   * will be consistent for every resolution
   * =======================================================================
   */
  
  double z, za, z0, tt, tp, tmp;
  int    face_num,jp,jm;
  long   ifp, ifm;
  int    ix, iy, ix_low, ix_hi, iy_low, iy_hi, ipf, ntt;
  double piover2 = 0.5*M_PI, pi = M_PI, twopi = 2.0*M_PI;
  int    ns_max = 8192;
  static int x2pix[128], y2pix[128];
  static char setup_done = 0;
  
  if( nside<1 || nside>ns_max ) {
    throw std::runtime_error("nside out of range");
  }
  if( theta<0 || theta>pi ) {
    throw std::runtime_error("theta out of range");
  }
  if( !setup_done ) {
    mk_xy2pix(x2pix,y2pix);
    setup_done = 1;
  }
  
  z  = cos(theta);
  za = fabs(z);
  z0 = 2./3.;
  if( phi>=twopi ) phi = phi - twopi;
  if( phi<0. )    phi = phi + twopi;
  tt = phi / piover2; /* in [0,4[ */
  
  if( za<=z0 ) { /* equatorial region */
    
    /* (the index of edge lines increase when the longitude=phi goes up) */
    jp = (int)floor(ns_max*(0.5 + tt - z*0.75)); /* ascending edge line index */
    jm = (int)floor(ns_max*(0.5 + tt + z*0.75)); /* descending edge line index */
    
    /* finds the face */
    ifp = jp / ns_max; /* in {0,4} */
    ifm = jm / ns_max;
    
    /* Note by Bruce L:  All statements in all functions in this module
                         of the following form 
    
                           if( ifp==ifm ) face_num = (int)fmod(ifp, 4) + 4; 
     
                         which have integer arguments were were revised to use
                         the % operator so that they would compile.
     */ 
    
    if( ifp==ifm ) face_num = (ifp % 4) + 4; /* faces 4 to 7 */
    else if( ifp<ifm ) face_num = ifp % 4;   /* (half-)faces 0 to 3 */
    else face_num = (ifm % 4) + 8;           /* (half-)faces 8 to 11 */
    
    ix = jm % ns_max;
    iy = ns_max - (jp % ns_max) - 1;
  }
  else { /* polar region, za > 2/3 */
    
    ntt = (int)floor(tt);
    if( ntt>=4 ) ntt = 3;
    tp = tt - ntt;
    tmp = sqrt( 3.*(1. - za) ); /* in ]0,1] */
    
    /* (the index of edge lines increase when distance from the closest pole
     * goes up)
     */
    /* line going toward the pole as phi increases */
    jp = (int)floor( ns_max * tp          * tmp ); 

    /* that one goes away of the closest pole */
    jm = (int)floor( ns_max * (1. - tp) * tmp );
    jp = (int)(jp < ns_max-1 ? jp : ns_max-1);
    jm = (int)(jm < ns_max-1 ? jm : ns_max-1);
    
    /* finds the face and pixel's (x,y) */
    if( z>=0 ) {
      face_num = ntt; /* in {0,3} */
      ix = ns_max - jm - 1;
      iy = ns_max - jp - 1;
    }
    else {
      face_num = ntt + 8; /* in {8,11} */
      ix =  jp;
      iy =  jm;
    }
  }
  
  ix_low = ix % 128;
  ix_hi  =     ix/128;
  iy_low = iy % 128;
  iy_hi  =     iy/128;

  ipf = (x2pix[ix_hi]+y2pix[iy_hi]) * (128 * 128)+ (x2pix[ix_low]+y2pix[iy_low]);
  ipf = (long)(ipf / pow(ns_max/nside,2));     /* in {0, nside**2 - 1} */
  *ipix =(long)( ipf + face_num*pow(nside,2)); /* in {0, 12*nside**2 - 1} */
}



void ang2pix_ring( const long nside, double theta, double phi, long *ipix) {
  /*
    c=======================================================================
    c     gives the pixel number ipix (RING) 
    c     corresponding to angles theta and phi
    c=======================================================================
  */
  
  int nl2, nl4, ncap, npix, jp, jm, ipix1;
  double  z, za, tt, tp, tmp;
  int ir, ip, kshift;
  
  double piover2 = 0.5*M_PI;
  double PI=M_PI;
  double twopi=2.0*M_PI;
  double z0=2.0/3.0;
  long ns_max=8192;
  
  if( nside<1 || nside>ns_max ) {
    throw std::runtime_error("nside out of range");
  }
  
  if( theta<0. || theta>PI) {
    throw std::runtime_error("theta out of range");
  }
  
  z = cos(theta);
  za = fabs(z);
  if( phi >= twopi)  phi = phi - twopi;
  if (phi < 0.)     phi = phi + twopi;
  tt = phi / piover2;//  ! in [0,4)
  
  nl2 = 2*nside;
  nl4 = 4*nside;
  ncap  = nl2*(nside-1);// ! number of pixels in the north polar cap
  npix  = 12*nside*nside;
  
  if( za <= z0 ) {
    
    jp = (int)floor(nside*(0.5 + tt - z*0.75)); /*index of ascending edge line*/
    jm = (int)floor(nside*(0.5 + tt + z*0.75)); /*index of descending edge line*/
    
    ir = nside + 1 + jp - jm;// ! in {1,2n+1} (ring number counted from z=2/3)
    kshift = 0;
    if ((ir % 2)==0.) kshift = 1;// ! kshift=1 if ir even, 0 otherwise
    
    ip = (int)floor( (double)(( jp+jm - nside + kshift + 1 ) / 2 )) + 1;// ! in {1,4n}
    if( ip>nl4 ) ip = ip - nl4;
    
    ipix1 = ncap + nl4*(ir-1) + ip ;
  }
  else {
    
    tp = tt - floor(tt);//      !MOD(tt,1.d0)
    tmp = sqrt( 3.*(1. - za) );
    
    jp = (int)floor( nside * tp * tmp );// ! increasing edge line index
    jm = (int)floor( nside * (1. - tp) * tmp );// ! decreasing edge line index
    
    ir = jp + jm + 1;//        ! ring number counted from the closest pole
    ip = (int)floor( tt * ir ) + 1;// ! in {1,4*ir}
    if( ip>4*ir ) ip = ip - 4*ir;
    
    ipix1 = 2*ir*(ir-1) + ip;
    if( z<=0. ) {
      ipix1 = npix - 2*ir*(ir+1) + ip;
    }
  }
  *ipix = ipix1 - 1;// ! in {0, npix-1}
  
}







void mk_pix2xy(int *pix2x, int *pix2y) {

  /* =======================================================================
   * subroutine mk_pix2xy
   * =======================================================================
   * constructs the array giving x and y in the face from pixel number
   * for the nested (quad-cube like) ordering of pixels
   *
   * the bits corresponding to x and y are interleaved in the pixel number
   * one breaks up the pixel number by even and odd bits
   * =======================================================================
   */

  int i, kpix, jpix, IX, IY, IP, ID;
  for (i = 0; i < 1023; i++) pix2x[i]=0;
  
  for( kpix=0;kpix<1024;kpix++ ) {
    jpix = kpix;
    IX = 0;
    IY = 0;
    IP = 1 ;//              ! bit position (in x and y)
    while( jpix!=0 ){// ! go through all the bits
      ID = jpix % 2;//  ! bit value (in kpix), goes in ix
      jpix = jpix/2;
      IX = ID*IP+IX;
      
      ID = jpix % 2;//  ! bit value (in kpix), goes in iy
      jpix = jpix/2;
      IY = ID*IP+IY;
      
      IP = 2*IP;//         ! next bit (in x and y)
    }
    
    pix2x[kpix] = IX;//     ! in 0,31
    pix2y[kpix] = IY;//     ! in 0,31
  }
  
  /* Later */
  return;
}





void nest2ring( long nside, long int ipnest, long *ipring) {
  /*
    c=======================================================================
    subroutine nest2ring(nside, ipnest, ipring)
    c=======================================================================
    c     conversion from NESTED to RING pixel number
    c=======================================================================
  */
      int npix, npface, face_num, ncap, n_before;
      int ipf, ip_low, ip_trunc, ip_med, ip_hi;
      int ix, iy, jrt, jr, nr, jpt, jp, kshift, nl4;
      int ns_max=8192;

      static int pix2x[1024], pix2y[1024];
      static char nest2string_setup_done = 0;

      int jrll[12], jpll[12];// ! coordinate of the lowest corner of each face
      //      data jrll/2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4/ ! in unit of nside
      //      data jpll/1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7/ ! in unit of nside/2
      jrll[0]=2;
      jrll[1]=2;
      jrll[2]=2;
      jrll[3]=2;
      jrll[4]=3;
      jrll[5]=3;
      jrll[6]=3;
      jrll[7]=3;
      jrll[8]=4;
      jrll[9]=4;
      jrll[10]=4;
      jrll[11]=4;
      jpll[0]=1;
      jpll[1]=3;
      jpll[2]=5;
      jpll[3]=7;
      jpll[4]=0;
      jpll[5]=2;
      jpll[6]=4;
      jpll[7]=6;
      jpll[8]=1;
      jpll[9]=3;
      jpll[10]=5;
      jpll[11]=7;

      if( nside<1 || nside>ns_max ) {
	throw std::runtime_error("nside out of range");
      }
      npix = 12 * nside*nside;
      if( ipnest<0 || ipnest>npix-1 ) {
	throw std::runtime_error("ipnest out of range");
      }

      //c     initiates the array for the pixel number -> (x,y) mapping
      if ( !nest2string_setup_done ) {
	 mk_pix2xy(pix2x,pix2y);
	 nest2string_setup_done = 1;
      }

      ncap  = 2*nside*(nside-1);// ! number of points in the North Polar cap
      nl4   = 4*nside;

      //c     finds the face, and the number in the face
      npface = nside*nside;
      //cccccc      ip = ipnest - 1         ! in {0,npix-1}

      face_num = ipnest/npface;//  ! face number in {0,11}
      ipf = ipnest % npface;//  ! pixel number in the face {0,npface-1}

      //c     finds the x,y on the face (starting from the lowest corner)
      //c     from the pixel number
      ip_low = ipf % 1024;//       ! content of the last 10 bits
      ip_trunc =   ipf/1024;//        ! truncation of the last 10 bits
      ip_med = ip_trunc % 1024;//  ! content of the next 10 bits
      ip_hi  =     ip_trunc/1024;//   ! content of the high weight 10 bits

      //      ix = 1024*pix2x[ip_hi] + 32*pix2x[ip_med] + pix2x[ip_low];
      //      iy = 1024*pix2y[ip_hi] + 32*pix2y[ip_med] + pix2y[ip_low];
      ix = 1024*pix2x[ip_hi] + 32*pix2x[ip_med] + pix2x[ip_low];
      iy = 1024*pix2y[ip_hi] + 32*pix2y[ip_med] + pix2y[ip_low];
      //      cout << "ix = " << ix << " iy = " << iy << endl;

      //c     transforms this in (horizontal, vertical) coordinates
      jrt = ix + iy;//  ! 'vertical' in {0,2*(nside-1)}
      jpt = ix - iy;//  ! 'horizontal' in {-nside+1,nside-1}

      //c     computes the z coordinate on the sphere
      //      jr =  jrll[face_num+1]*nside - jrt - 1;//   ! ring number in {1,4*nside-1}
      jr =  jrll[face_num]*nside - jrt - 1;

      nr = nside;//                  ! equatorial region (the most frequent)
      n_before = ncap + nl4 * (jr - nside);
      kshift = (jr - nside) % 2;
      if( jr<nside ) {//then     ! north pole region
         nr = jr;
         n_before = 2 * nr * (nr - 1);
         kshift = 0;
      }
      else if( jr>3*nside ) {//then ! south pole region
         nr = nl4 - jr;
         n_before = npix - 2 * (nr + 1) * nr;
         kshift = 0;
      }

      //c     computes the phi coordinate on the sphere, in [0,2Pi]
      //      jp = (jpll[face_num+1]*nr + jpt + 1 + kshift)/2;//  ! 'phi' number in the ring in {1,4*nr}
      jp = (jpll[face_num]*nr + jpt + 1 + kshift)/2;

      if( jp>nl4 ) jp = jp - nl4;
      if( jp<1 )   jp = jp + nl4;

      *ipring = n_before + jp - 1;// ! in {0, npix-1}

}


long nside2npix(const long nside) {
  return 12*nside*nside;
}


void pix2ang_nest( long nside, long ipix, double *theta, double *phi) {

  /*
    c=======================================================================
    subroutine pix2ang_nest(nside, ipix, theta, phi)
    c=======================================================================
    c     gives theta and phi corresponding to pixel ipix (NESTED) 
    c     for a parameter nside
    c=======================================================================
  */
    
      int npix, npface, face_num;
      int  ipf, ip_low, ip_trunc, ip_med, ip_hi;
      int     ix, iy, jrt, jr, nr, jpt, jp, kshift, nl4;
      double z, fn, fact1, fact2;
      double piover2=0.5*M_PI;
      int ns_max=8192;
      
      static int pix2x[1024], pix2y[1024];
      //      common /pix2xy/ pix2x, pix2y
      
      int jrll[12], jpll[12];// ! coordinate of the lowest corner of each face
      //      data jrll/2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4/ ! in unit of nside
      //      data jpll/1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7/ ! in unit of nside/2
      jrll[0]=2;
      jrll[1]=2;
      jrll[2]=2;
      jrll[3]=2;
      jrll[4]=3;
      jrll[5]=3;
      jrll[6]=3;
      jrll[7]=3;
      jrll[8]=4;
      jrll[9]=4;
      jrll[10]=4;
      jrll[11]=4;
      jpll[0]=1;
      jpll[1]=3;
      jpll[2]=5;
      jpll[3]=7;
      jpll[4]=0;
      jpll[5]=2;
      jpll[6]=4;
      jpll[7]=6;
      jpll[8]=1;
      jpll[9]=3;
      jpll[10]=5;
      jpll[11]=7;
      
      
      if( nside<1 || nside>ns_max ) {
	throw std::runtime_error("nside out of range");
      }
      npix = 12 * nside*nside;
      if( ipix<0 || ipix>npix-1 ) {
	throw std::runtime_error("ipx out of range");
      }

      /* initiates the array for the pixel number -> (x,y) mapping */
      if( pix2x[1023]<=0 ) mk_pix2xy(pix2x,pix2y);

      fn = 1.*nside;
      fact1 = 1./(3.*fn*fn);
      fact2 = 2./(3.*fn);
      nl4   = 4*nside;

      //c     finds the face, and the number in the face
      npface = nside*nside;

      face_num = ipix/npface;//  ! face number in {0,11}
      ipf = ipix % npface;//  ! pixel number in the face {0,npface-1}

      //c     finds the x,y on the face (starting from the lowest corner)
      //c     from the pixel number
      ip_low = ipf % 1024;//       ! content of the last 10 bits
      ip_trunc =   ipf/1024 ;//       ! truncation of the last 10 bits
      ip_med = ip_trunc % 1024;//  ! content of the next 10 bits
      ip_hi  =     ip_trunc/1024   ;//! content of the high weight 10 bits

      ix = 1024*pix2x[ip_hi] + 32*pix2x[ip_med] + pix2x[ip_low];
      iy = 1024*pix2y[ip_hi] + 32*pix2y[ip_med] + pix2y[ip_low];

      //c     transforms this in (horizontal, vertical) coordinates
      jrt = ix + iy;//  ! 'vertical' in {0,2*(nside-1)}
      jpt = ix - iy;//  ! 'horizontal' in {-nside+1,nside-1}

      //c     computes the z coordinate on the sphere
      //      jr =  jrll[face_num+1]*nside - jrt - 1;//   ! ring number in {1,4*nside-1}
      jr =  jrll[face_num]*nside - jrt - 1;
      //      cout << "face_num=" << face_num << endl;
      //      cout << "jr = " << jr << endl;
      //      cout << "jrll(face_num)=" << jrll[face_num] << endl;
      //      cout << "----------------------------------------------------" << endl;
      nr = nside;//                  ! equatorial region (the most frequent)
      z  = (2*nside-jr)*fact2;
      kshift = (jr - nside) % 2;
      if( jr<nside ) { //then     ! north pole region
         nr = jr;
         z = 1. - nr*nr*fact1;
         kshift = 0;
      }
      else {
	if( jr>3*nside ) {// then ! south pole region
         nr = nl4 - jr;
         z = - 1. + nr*nr*fact1;
         kshift = 0;
	}
      }
      *theta = acos(z);
      
      //c     computes the phi coordinate on the sphere, in [0,2Pi]
      //      jp = (jpll[face_num+1]*nr + jpt + 1 + kshift)/2;//  ! 'phi' number in the ring in {1,4*nr}
      jp = (jpll[face_num]*nr + jpt + 1 + kshift)/2;
      if( jp>nl4 ) jp = jp - nl4;
      if( jp<1 )   jp = jp + nl4;

      *phi = (jp - (kshift+1)*0.5) * (piover2 / nr);

}



void pix2ang_ring( long nside, long ipix, double *theta, double *phi) {
  /*
    c=======================================================================
    c     gives theta and phi corresponding to pixel ipix (RING) 
    c     for a parameter nside
    c=======================================================================
  */
  
  int nl2, nl4, npix, ncap, iring, iphi, ip, ipix1;
  double  fact1, fact2, fodd, hip, fihip;
  double PI=M_PI;
  //      PARAMETER (pi     = 3.1415926535897932384626434d0)
  //      parameter (ns_max = 8192) ! 2^13 : largest nside available
  
  int ns_max=8192;
  
  if( nside<1 || nside>ns_max ) {
    throw std::runtime_error("nside out of range");
  }
  npix = 12*nside*nside;      // ! total number of points
  if( ipix<0 || ipix>npix-1 ) {
    throw std::runtime_error("ipx out of range");
  }
  
  ipix1 = ipix + 1; // in {1, npix}
  nl2 = 2*nside;
  nl4 = 4*nside;
  ncap = 2*nside*(nside-1);// ! points in each polar cap, =0 for nside =1
  fact1 = 1.5*nside;
  fact2 = 3.0*nside*nside;
  
  if( ipix1 <= ncap ) {  //! North Polar cap -------------
    
    hip   = ipix1/2.;
    fihip = floor(hip);
    iring = (int)floor( sqrt( hip - sqrt(fihip) ) ) + 1;// ! counted from North pole
    iphi  = ipix1 - 2*iring*(iring - 1);
    
    *theta = acos( 1. - iring*iring / fact2 );
    *phi   = (1.*iphi - 0.5) * PI/(2.*iring);
  }
  else if( ipix1 <= nl2*(5*nside+1) ) {//then ! Equatorial region ------
    
    ip    = ipix1 - ncap - 1;
    iring = (int)floor( (double) (ip / nl4) ) + nside;// ! counted from North pole
    iphi  = (ip % nl4) + 1;
    
    fodd  = 0.5 * (1 + fmod((double)(iring+nside),2));//  ! 1 if iring+nside is odd, 1/2 otherwise
    *theta = acos( (nl2 - iring) / fact1 );
    *phi   = (1.*iphi - fodd) * PI /(2.*nside);
  }
  else {//! South Polar cap -----------------------------------
    
    ip    = npix - ipix1 + 1;
    hip   = ip/2.;
/* bug corrige floor instead of 1.* */
    fihip = floor(hip);
    iring = (int)floor( sqrt( hip - sqrt(fihip) ) ) + 1;//     ! counted from South pole
    iphi  = (int)(4.*iring + 1 - (ip - 2.*iring*(iring-1)));
    
    *theta = acos( -1. + iring*iring / fact2 );
    *phi   = (1.*iphi - 0.5) * PI/(2.*iring);
  }
}
void ring2nest( long nside, long ipring, long *ipnest) {
  /*
    c=======================================================================
    subroutine ring2nest(nside, ipring, ipnest)
    c=======================================================================
    c     conversion from RING to NESTED pixel number
    c=======================================================================
  */
  
  double fihip, hip;
  int npix, nl2, nl4, ncap, ip, iphi, ipt, ipring1;
  int     kshift, face_num, nr;
  int irn, ire, irm, irs, irt, ifm , ifp;
  int ix, iy, ix_low, ix_hi, iy_low, iy_hi, ipf;
  int ns_max=8192;
  
  static int x2pix[128], y2pix[128];
  //      common    /xy2pix/ x2pix,y2pix

  int jrll[12], jpll[12];// ! coordinate of the lowest corner of each face
  //      data jrll/2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4/ ! in unit of nside
  //      data jpll/1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7/ ! in unit of nside/2
  jrll[0]=2;
  jrll[1]=2;
  jrll[2]=2;
  jrll[3]=2;
  jrll[4]=3;
  jrll[5]=3;
  jrll[6]=3;
  jrll[7]=3;
  jrll[8]=4;
  jrll[9]=4;
  jrll[10]=4;
  jrll[11]=4;
  jpll[0]=1;
  jpll[1]=3;
  jpll[2]=5;
  jpll[3]=7;
  jpll[4]=0;
  jpll[5]=2;
  jpll[6]=4;
  jpll[7]=6;
  jpll[8]=1;
  jpll[9]=3;
  jpll[10]=5;
  jpll[11]=7;
  
  if( nside<1 || nside>ns_max ) {
    throw std::runtime_error("nside out of range");
  }
  npix = 12 * nside*nside;
  if( ipring<0 || ipring>npix-1 ) {
    throw std::runtime_error("ipring out of range");
  }
  if( x2pix[127]<=0 ) mk_xy2pix(x2pix,y2pix);
  
  nl2 = 2*nside;
  nl4 = 4*nside;
  npix = 12*nside*nside;//      ! total number of points
  ncap = 2*nside*(nside-1);// ! points in each polar cap, =0 for nside =1
  ipring1 = ipring + 1;
  
  //c     finds the ring number, the position of the ring and the face number
  if( ipring1<=ncap ) {//then
    
    hip   = ipring1/2.;
    fihip = (int)floor ( hip );
    irn   = (int)floor( sqrt( hip - sqrt(fihip) ) ) + 1;// ! counted from North pole
    iphi  = ipring1 - 2*irn*(irn - 1);
    
    kshift = 0;
    nr = irn   ;//               ! 1/4 of the number of points on the current ring
    face_num = (iphi-1) / irn;// ! in {0,3}
  }
  else if( ipring1<=nl2*(5*nside+1) ) {//then
    
    ip    = ipring1 - ncap - 1;
    irn   = (int)floor( (double)(ip / nl4) ) + nside;//               ! counted from North pole
    iphi  = (ip % nl4) + 1;
    
    kshift  = ((irn+nside) % 2);//  ! 1 if irn+nside is odd, 0 otherwise
    nr = nside;
    ire =  irn - nside + 1;// ! in {1, 2*nside +1}
    irm =  nl2 + 2 - ire;
    ifm = (iphi - ire/2 + nside -1) / nside;// ! face boundary
    ifp = (iphi - irm/2 + nside -1) / nside;
    if( ifp==ifm ) {//then          ! faces 4 to 7
      face_num = (ifp % 4) + 4;
    }
    else if( ifp + 1==ifm ) {//then ! (half-)faces 0 to 3
      face_num = ifp;
    }
    else if( ifp - 1==ifm ) {//then ! (half-)faces 8 to 11
      face_num = ifp + 7;
    }
  }
  else {
    
    ip    = npix - ipring1 + 1;
    hip   = ip/2.;
    fihip = floor ( hip );
    irs   = (int)floor( sqrt( hip - sqrt(fihip) ) ) + 1;//  ! counted from South pole
    iphi  = 4*irs + 1 - (ip - 2*irs*(irs-1));
    
    kshift = 0;
    nr = irs;
    irn   = nl4 - irs;
    face_num = (iphi-1) / irs + 8;// ! in {8,11}
  }
  
  //c     finds the (x,y) on the face
  //  irt =   irn  - jrll[face_num+1]*nside + 1;//       ! in {-nside+1,0}
  //  ipt = 2*iphi - jpll[face_num+1]*nr - kshift - 1;// ! in {-nside+1,nside-1}
  irt =   irn  - jrll[face_num]*nside + 1;//       ! in {-nside+1,0}
  ipt = 2*iphi - jpll[face_num]*nr - kshift - 1;


  if( ipt>=nl2 ) ipt = ipt - 8*nside;// ! for the face #4
  
  ix =  (ipt - irt ) / 2;
  iy = -(ipt + irt ) / 2;
  
  ix_low = (ix % 128);
  ix_hi  = ix/128;
  iy_low = (iy % 128);
  iy_hi  = iy/128;
  //  cout << "ix_low = " << ix_low << " ix_hi = " << ix_hi << endl;
  //  cout << "iy_low = " << iy_low << " iy_hi = " << iy_hi << endl;
  //  ipf =  (x2pix[ix_hi +1]+y2pix[iy_hi +1]) * (128 * 128)
  //    + (x2pix[ix_low+1]+y2pix[iy_low+1]);//        ! in {0, nside**2 - 1}
  ipf =  (x2pix[ix_hi]+y2pix[iy_hi]) * (128 * 128)
    + (x2pix[ix_low]+y2pix[iy_low]);


  //  cout << "ipf = " << ipf << endl;
  //  for( int i(0);i<128;i++ ) cout << x2pix[i] << " || " << y2pix[i] << endl;
  *ipnest = ipf + face_num* nside *nside;//   ! in {0, 12*nside**2 - 1}
  
}

    /*=======================================================================
    !     gives the x, y coords in a face from pixel number within the face (NESTED) 
    !
    !     Benjamin D. Wandelt 13/10/97
    !     
    !     using code from HEALPIX toolkit by K.Gorski and E. Hivon
    !======================================================================= */
void pix2xy_nest(long nside, long ipf, long * ix, long * iy)
{
    long  ip_low, ip_trunc, ip_med, ip_hi;

    if (nside < 1 || nside > NS_MAX)
        throw std::runtime_error("nside out of range");
 
    if (ipf < 0 || ipf > (nside*nside) -1)
        throw std::runtime_error("input pixel out of range"); 

    static int pix2x[1024], pix2y[1024];
    static char setup_done = 0;
    
    if( !setup_done )
    {
        mk_pix2xy(pix2x, pix2y);
        setup_done = 1;
    }

    ip_low = ipf % 1024;      // content of the last 10 bits
    ip_trunc = ipf/1024;       // truncation of the last 10 bits
    ip_med = ip_trunc % 1024; // content of the next 10 bits
    ip_hi = ip_trunc/1024;    // content of the high weight 10 bits

    *ix = 1024*pix2x[ip_hi] + 32*pix2x[ip_med] + pix2x[ip_low];
    *iy = 1024*pix2y[ip_hi] + 32*pix2y[ip_med] + pix2y[ip_low];
}
void
xy2pix_nest(long nside, long ix, long iy, long face_num, long * ipix)
{
    /*=======================================================================
    !     gives the pixel number ipix (NESTED) 
    !     corresponding to ix, iy and face_num
    !
    !     Benjamin D. Wandelt 13/10/97
    !     using code from HEALPIX toolkit by K.Gorski and E. Hivon
    !======================================================================= */
    int  ix_low, ix_hi, iy_low, iy_hi, ipf;

    if (nside < 1 || nside > NS_MAX)
        throw std::runtime_error("nside out of range");
    if (ix < 0 || ix > (nside-1))
        throw std::runtime_error("ix out of range");
    if (iy < 0 || iy > (nside-1))
        throw std::runtime_error("iy out of range");
    
    static int x2pix[128], y2pix[128];
    static bool setup_done = false;
    
    if( !setup_done )
    {
        mk_xy2pix(x2pix, y2pix);
        setup_done = true;
    }

    ix_low = ix % 128;
    ix_hi  = ix/128;
    iy_low = iy % 128;
    iy_hi  = iy/128;

    ipf =  (x2pix[ix_hi]+y2pix[iy_hi]) * (128 * 128) 
              + (x2pix[ix_low]+y2pix[iy_low]);

    * ipix = (long) (ipf + (face_num * nside * nside));    // in {0, 12*nside**2 - 1}
}

int swapLSBMSB(int i)
{
  /* ====================================================================
  !     Returns i with even and odd bit positions interchanged
  ! Benjamin D. Wandelt October 1997
  !==================================================================== */
    int msb,lsb, magic1=89478485, magic2=178956970;
    
    lsb = i & magic1;
    msb = i & magic2;

    return((msb/2)+(lsb*2));
}
int invLSB( int i)
{
  /* ====================================================================
  !     Returns i with even (0,2,4,...) bits inverted
  ! Benjamin D. Wandelt October 1997
  !==================================================================== */
    int magic1=89478485;
    
    return(i ^ magic1);
}
int invMSB( int i)
{
  /* ====================================================================
  !     Returns i with odd (1,3,5,...) bits inverted
  ! Benjamin D. Wandelt October 1997
  !==================================================================== */
    int magic2=178956970;
    
    return(i ^ magic2);
}  
 //=======================================================================
// The following is a routine which finds the 7 or 8 neighbours of 
// any pixel in the nested scheme of the HEALPIX pixelisation.
//====================================================================
//====================================================================
//   Returns list n(8) of neighbours of pixel ipix (in NESTED scheme)
//   the neighbours are ordered in the following way:
//   First pixel is the one to the south (the one west of the south
// direction is taken
// for the pixels which don't have a southern neighbour). From
// then on the neighbours are ordered in the clockwise direction
// about the pixel with number ipix.
//     
//   nneigh is the number of neighbours (mostly 8, 8 pixels have 7 neighbours)
//
//   Benjamin D. Wandelt October 1997
//   Added to pix_tools in March 1999
//====================================================================
void neighbours_nest(long nside, long ipix, long * n, int * nneigh)
{
    long npix, ipf, ipo, ix, ixm, ixp, iy, iym, iyp, ixo, iyo;
    long face_num, other_face;
    long ia, ib, ibp, ibm, ib2, icase, nsidesq;
    long local_magic1, local_magic2;

    if (nside < 1 || nside > NS_MAX)
        throw std::runtime_error("nside out of range");
    
    nsidesq = nside * nside;
    npix = 12 * nsidesq;       // total number of points
    if (ipix < 0 || ipix > npix-1)
        throw std::runtime_error("input pixel out of range");

    //     initiates array for (x,y)-> pixel number -> (x,y) mapping
    static int x2pix[128], y2pix[128];
    static bool setup_done = false;
    
    if( !setup_done )
    {
        mk_xy2pix(x2pix, y2pix);
        setup_done = true;
    }

    local_magic1 = (nsidesq-1)/3;
    local_magic2 = 2* local_magic1;
    face_num = ipix / nsidesq;

    // ipf = modulo (ipix,nsidesq)   // Pixel number in face
    ipf = ipix % nsidesq;

    pix2xy_nest(nside, ipf, &ix, &iy);
    ixm=ix-1;
    ixp=ix+1;
    iym=iy-1;
    iyp=iy+1;

    *nneigh=8;                 // Except in special cases below
    bool special = false;

    //     Exclude corners
    if(ipf==local_magic2)      // WestCorner
    {
       icase=5;             
       special = true;
    }
    else if(ipf==(nsidesq-1))  // NorthCorner
    {
       icase=6;             
       special = true;
    }
    else if(ipf==0)            // SouthCorner
    {
       icase=7;                
       special = true;
    }
    else if(ipf==local_magic1)      // EastCorner
    {
       icase=8;      
       special = true;
    }

    //     Detect edges
    else if((ipf & local_magic1)==local_magic1)  // NorthEast
    {
       icase=1; 
       special = true;
    }
    else if((ipf & local_magic1)==0)       // SouthWest
    {
       icase=2;    
       special = true;
    }
    else if((ipf & local_magic2)==local_magic2)  // NorthWest
    {
       icase=3;                 
       special = true;
    }
    else if((ipf & local_magic2)==0)       // SouthEast
    {
       icase=4;                   
       special = true;
    }

    //     Inside a face
    if (! special)
    {
        xy2pix_nest(nside, ixm, iym, face_num, n);
        xy2pix_nest(nside, ixm, iy , face_num, n + 1);  
        xy2pix_nest(nside, ixm, iyp, face_num, n + 2);
        xy2pix_nest(nside, ix , iyp, face_num, n + 3);
        xy2pix_nest(nside, ixp, iyp, face_num, n + 4);
        xy2pix_nest(nside, ixp, iy , face_num, n + 5);
        xy2pix_nest(nside, ixp, iym, face_num, n + 6);
        xy2pix_nest(nside, ix , iym, face_num, n + 7);
    }

    else  // Handle special cases.
    {
        ia= face_num/4;            // in {0,2}
        ib= face_num % 4;       // in {0,3}
        ibp=(ib+1) % 4;
        ibm=(ib+4-1) % 4;
        ib2=(ib+2) % 4;

        if(ia==0)           // North Pole region
        {
            switch (icase)
            {
                case 1:              // NorthEast edge
                    other_face=0+ibp;
                    xy2pix_nest(nside, ix , iym, face_num, n + 7);
                    xy2pix_nest(nside, ixm, iym, face_num, n);
                    xy2pix_nest(nside, ixm, iy , face_num, n + 1);
                    xy2pix_nest(nside, ixm, iyp, face_num, n + 2);
                    xy2pix_nest(nside, ix , iyp, face_num, n + 3);        
                    ipo = swapLSBMSB(ipf) % nsidesq;    // East-West flip
                    pix2xy_nest(nside,ipo, & ixo, & iyo);
                    xy2pix_nest(nside, ixo+1 , iyo, other_face, n + 4);
                    *(n + 5)=other_face*nsidesq+ipo;
                    xy2pix_nest(nside, ixo-1, iyo, other_face, n + 6);
                break;
                case 2:              // SouthWest edge
                    other_face=4+ib;
                    ipo = invLSB(ipf) % nsidesq;        // SW-NE flip
                    pix2xy_nest(nside,ipo,&ixo,&iyo);
                    xy2pix_nest(nside, ixo, iyo-1, other_face, n);
                    *(n + 1)=other_face*nsidesq+ipo;
                    xy2pix_nest(nside, ixo, iyo+1, other_face, n + 2);
                    xy2pix_nest(nside, ix , iym, face_num, n + 7);
                    xy2pix_nest(nside, ix , iyp, face_num, n + 3);
                    xy2pix_nest(nside, ixp, iym, face_num, n + 6);
                    xy2pix_nest(nside, ixp, iy , face_num, n + 5);
                    xy2pix_nest(nside, ixp, iyp, face_num, n + 4);
                break;
                case 3:              // NorthWest edge
                    other_face=0+ibm;
                    ipo = swapLSBMSB(ipf) % nsidesq;    // East-West flip
                    pix2xy_nest(nside,ipo,&ixo,&iyo);
                    xy2pix_nest(nside, ixo, iyo-1, other_face, n + 2);
                    *(n + 3)=other_face*nsidesq+ipo;
                    xy2pix_nest(nside, ixo, iyo+1, other_face, n + 4);
                    xy2pix_nest(nside, ixm, iym, face_num, n);
                    xy2pix_nest(nside, ixm, iy , face_num, n + 1);
                    xy2pix_nest(nside, ix , iym, face_num, n + 7);
                    xy2pix_nest(nside, ixp, iym, face_num, n + 6);
                    xy2pix_nest(nside, ixp, iy , face_num, n + 5);
                break;
                case 4:              //SouthEast edge
                    other_face=4+ibp;
                    xy2pix_nest(nside, ixm, iy , face_num, n+1);   
                    xy2pix_nest(nside, ixm, iyp, face_num, n+2);
                    xy2pix_nest(nside, ix , iyp, face_num, n+3);
                    xy2pix_nest(nside, ixp, iyp, face_num, n+4);
                    xy2pix_nest(nside, ixp, iy , face_num, n+5);
                    ipo = invMSB(ipf) % nsidesq; // SE-NW flip
                    pix2xy_nest(nside,ipo,&ixo,&iyo);
                    xy2pix_nest(nside, ixo+1, iyo, other_face, n+6);
                    *(n + 7)=other_face*nsidesq+ipo;
                    xy2pix_nest(nside, ixo-1, iyo, other_face, n);
                break;
                case 5:              //West corner
                    *nneigh=7;
                    other_face=4+ib;
                    *(n + 1)=other_face*nsidesq+nsidesq-1;
                    *n =*(n + 1)-2;
                    other_face=0+ibm;
                    *(n + 2)=other_face*nsidesq+local_magic1;
                    *(n + 3)=*(n + 2)+2;
                    *(n + 4)=ipix+1;
                    *(n + 5)=ipix-1;
                    *(n + 6)=ipix-2;
                break;
                case 6:              //North corner
                    *n=ipix-3;
                    *(n + 1)=ipix-1;
                    *(n + 7)=ipix-2;
                    other_face=0+ibm;
                    *(n + 3)=other_face*nsidesq+nsidesq-1;
                    *(n + 2)=*(n + 3)-2;
                    other_face=0+ib2;
                    *(n + 4)=other_face*nsidesq+nsidesq-1;
                    other_face=0+ibp;
                    *(n + 5)=other_face*nsidesq+nsidesq-1;
                    *(n + 6)=*(n + 5)-1;
                break;
                case 7:              //South corner
                    other_face=8+ib;
                    *n = other_face*nsidesq+nsidesq-1;
                    other_face=4+ib;
                    *(n + 1)=other_face*nsidesq+local_magic1;
                    *(n + 2)=*(n + 1)+2;
                    *(n + 3)=ipix+2;
                    *(n + 4)=ipix+3;
                    *(n + 5)=ipix+1;
                    other_face=4+ibp;
                    *(n + 7)=other_face*nsidesq+local_magic2;
                    *(n + 6)=*(n + 7)+1;
                break;
                case 8:              // East corner
                    *nneigh=7;
                    *(n + 1)=ipix-1;
                    *(n + 2)=ipix+1;
                    *(n + 3)=ipix+2;
                    other_face=0+ibp;
                    *(n + 5)=other_face*nsidesq+local_magic2;
                    *(n + 4)=*(n + 5)+1;
                    other_face=4+ibp;
                    *(n + 6)=other_face*nsidesq+nsidesq-1;
                    *n=*(n + 6)-1;
                break;
            } // End switch
        } //  North pole region

        else if(ia==1)       // Equatorial region
        {
            switch (icase)
            {
                case 1:              // NorthEast edge
                    other_face=0+ib;
                    xy2pix_nest(nside, ix , iym, face_num, n+7);
                    xy2pix_nest(nside, ixm, iym, face_num, n);
                    xy2pix_nest(nside, ixm, iy , face_num, n+1);
                    xy2pix_nest(nside, ixm, iyp, face_num, n+2);
                    xy2pix_nest(nside, ix , iyp, face_num, n+3);
                    ipo = invLSB(ipf) % nsidesq;    // NE-SW flip
                    pix2xy_nest(nside,ipo,&ixo,&iyo);
                    xy2pix_nest(nside, ixo , iyo+1, other_face, n+4);
                    *(n + 5)=other_face*nsidesq+ipo;
                    xy2pix_nest(nside, ixo, iyo-1, other_face, n+6);
                break;
                case 2:              // SouthWest edge
                    other_face=8+ibm;
                    ipo = invLSB(ipf) % nsidesq;        // SW-NE flip
                    pix2xy_nest(nside,ipo,&ixo,&iyo);
                    xy2pix_nest(nside, ixo, iyo-1, other_face, n);
                    *(n + 1)=other_face*nsidesq+ipo;
                    xy2pix_nest(nside, ixo, iyo+1, other_face, n+2);
                    xy2pix_nest(nside, ix , iym, face_num, n+7);
                    xy2pix_nest(nside, ix , iyp, face_num, n+3);
                    xy2pix_nest(nside, ixp, iym, face_num, n+6);
                    xy2pix_nest(nside, ixp, iy , face_num, n+5);
                    xy2pix_nest(nside, ixp, iyp, face_num, n+4);
                break;
                case 3:              // NorthWest edge
                    other_face=0+ibm;
                    ipo = invMSB(ipf) % nsidesq;    // NW-SE flip
                    pix2xy_nest(nside,ipo,&ixo,&iyo);
                    xy2pix_nest(nside, ixo-1, iyo, other_face, n+2);
                    *(n + 3)=other_face*nsidesq+ipo;
                    xy2pix_nest(nside, ixo+1, iyo, other_face, n+4);
                    xy2pix_nest(nside, ixm, iym, face_num, n);
                    xy2pix_nest(nside, ixm, iy , face_num, n+1);
                    xy2pix_nest(nside, ix , iym, face_num, n+7);
                    xy2pix_nest(nside, ixp, iym, face_num, n+6);
                    xy2pix_nest(nside, ixp, iy , face_num, n+5);
                break;
                case 4:              // SouthEast edge
                    other_face=8+ib;
                    xy2pix_nest(nside, ixm, iy , face_num, n+1)   ;
                    xy2pix_nest(nside, ixm, iyp, face_num, n+2);
                    xy2pix_nest(nside, ix , iyp, face_num, n+3);
                    xy2pix_nest(nside, ixp, iyp, face_num, n+4);
                    xy2pix_nest(nside, ixp, iy , face_num, n+5);
                    ipo = invMSB(ipf) % nsidesq; // SE-NW flip
                    pix2xy_nest(nside,ipo,&ixo,&iyo);
                    xy2pix_nest(nside, ixo+1, iyo, other_face, n+6);
                    *(n + 7)=other_face*nsidesq+ipo;
                    xy2pix_nest(nside, ixo-1, iyo, other_face, n);
                break;
                case 5:              // West corner
                    other_face=8+ibm;
                    *(n + 1)=other_face*nsidesq+nsidesq-1;
                    *n=*(n + 1)-2;
                    other_face=4+ibm;
                    *(n + 2)=other_face*nsidesq+local_magic1;
                    other_face=0+ibm;
                    *(n + 3)=other_face*nsidesq;
                    *(n + 4)=*(n + 3)+1;
                    *(n + 5)=ipix+1;
                    *(n + 6)=ipix-1;
                    *(n + 7)=ipix-2;
                break;
                case 6:             // North corner
                    *nneigh=7;
                    *n=ipix-3;
                    *(n + 1)=ipix-1;
                    other_face=0+ibm;
                    *(n + 3)=other_face*nsidesq+local_magic1;
                    *(n + 2)=*(n + 3)-1;
                    other_face=0+ib;
                    *(n + 4)=other_face*nsidesq+local_magic2;
                    *(n + 5)=*(n + 4)-2;
                    *(n + 6)=ipix-2;
                break;
                case 7:              // South corner
                    *nneigh=7;
                    other_face=8+ibm;
                    *n=other_face*nsidesq+local_magic1;
                    *(n + 1)=(*n) +2;
                    *(n + 2)=ipix+2;
                    *(n + 3)=ipix+3;
                    *(n + 4)=ipix+1;
                    other_face=8+ib;
                    *(n + 6)=other_face*nsidesq+local_magic2;
                    *(n + 5)=*(n + 6)+1;
                break;
                case 8:              // East corner
                    other_face=8+ib;
                    *(n + 7)=other_face*nsidesq+nsidesq-1;
                    *n=*(n + 7)-1;
                    *(n + 1)=ipix-1;
                    *(n + 2)=ipix+1;
                    *(n + 3)=ipix+2;
                    other_face=0+ib;
                    *(n + 5)=other_face*nsidesq;
                    *(n + 4)=*(n + 5)+2;
                    other_face=4+ibp;
                    *(n + 6)=other_face*nsidesq+local_magic2;
                break;
            } // End switch
        }
        
        else                    // South Pole region
        {
            switch (icase)
            {
                case 1:              // NorthEast edge
                    other_face=4+ibp;
                    xy2pix_nest(nside, ix , iym, face_num, n+7);
                    xy2pix_nest(nside, ixm, iym, face_num, n);
                    xy2pix_nest(nside, ixm, iy , face_num, n+1);
                    xy2pix_nest(nside, ixm, iyp, face_num, n+2);
                    xy2pix_nest(nside, ix , iyp, face_num, n+3);
                    ipo = invLSB(ipf) % nsidesq;    // NE-SW flip
                    pix2xy_nest(nside,ipo,&ixo,&iyo);
                    xy2pix_nest(nside, ixo , iyo+1, other_face, n+4);
                    *(n + 5)=other_face*nsidesq+ipo;
                    xy2pix_nest(nside, ixo, iyo-1, other_face, n+6);
                break;
                case 2:              // SouthWest edge
                    other_face=8+ibm;
                    ipo = swapLSBMSB(ipf) % nsidesq;        // W-E flip
                    pix2xy_nest(nside,ipo,&ixo,&iyo);
                    xy2pix_nest(nside, ixo-1, iyo, other_face, n);
                    *(n + 1)=other_face*nsidesq+ipo;
                    xy2pix_nest(nside, ixo+1, iyo, other_face, n+2);
                    xy2pix_nest(nside, ix , iym, face_num, n+7);
                    xy2pix_nest(nside, ix , iyp, face_num, n+3);
                    xy2pix_nest(nside, ixp, iym, face_num, n+6);
                    xy2pix_nest(nside, ixp, iy , face_num, n+5);
                    xy2pix_nest(nside, ixp, iyp, face_num, n+4);
                break;
                case 3:              // NorthWest edge
                    other_face=4+ib;
                    ipo = invMSB(ipf) % nsidesq;    // NW-SE flip
                    pix2xy_nest(nside,ipo,&ixo,&iyo);
                    xy2pix_nest(nside, ixo-1, iyo, other_face, n+2);
                    *(n + 3)=other_face*nsidesq+ipo;
                    xy2pix_nest(nside, ixo+1, iyo, other_face, n+4);
                    xy2pix_nest(nside, ixm, iym, face_num, n);
                    xy2pix_nest(nside, ixm, iy , face_num, n+1);
                    xy2pix_nest(nside, ix , iym, face_num, n+7);
                    xy2pix_nest(nside, ixp, iym, face_num, n+6);
                    xy2pix_nest(nside, ixp, iy , face_num, n+5);
                break;
                case 4:              // SouthEast edge
                    other_face=8+ibp;
                    xy2pix_nest(nside, ixm, iy , face_num, n+1);
                    xy2pix_nest(nside, ixm, iyp, face_num, n+2);
                    xy2pix_nest(nside, ix , iyp, face_num, n+3);
                    xy2pix_nest(nside, ixp, iyp, face_num, n+4);
                    xy2pix_nest(nside, ixp, iy , face_num, n+5);
                    ipo = swapLSBMSB(ipf) % nsidesq; // E-W flip
                    pix2xy_nest(nside,ipo,&ixo,&iyo);
                    xy2pix_nest(nside, ixo, iyo+1, other_face, n+6);
                    *(n + 7)=other_face*nsidesq+ipo;
                    xy2pix_nest(nside, ixo, iyo-1, other_face, n);
                break;
                case 5:              // West corner
                    *nneigh=7;
                    other_face=8+ibm;
                    *(n + 1)=other_face*nsidesq+local_magic1;
                    *n=*(n + 1)-1;
                    other_face=4+ib;
                    *(n + 2)=other_face*nsidesq;
                    *(n + 3)=*(n + 2)+1;
                    *(n + 4)=ipix+1;
                    *(n + 5)=ipix-1;
                    *(n + 6)=ipix-2;
                break;
                case 6:              // North corner
                    *n=ipix-3;
                    *(n + 1)=ipix-1;
                    other_face=4+ib;
                    *(n + 3)=other_face*nsidesq+local_magic1;
                    *(n + 2)=*(n + 3)-1;
                    other_face=0+ib;
                    *(n + 4)=other_face*nsidesq;
                    other_face=4+ibp;
                    *(n + 5)=other_face*nsidesq+local_magic2;
                    *(n + 6)=*(n + 5)-2;
                    *(n + 7)=ipix-2;
                break;
                case 7:              // South corner
                    other_face=8+ib2;
                    *n=other_face*nsidesq;
                    other_face=8+ibm;
                    *(n + 1)=other_face*nsidesq;
                    *(n + 2)=*(n + 1)+1;
                    *(n + 3)=ipix+2;
                    *(n + 4)=ipix+3;
                    *(n + 5)=ipix+1;
                    other_face=8+ibp;
                    *(n + 7)=other_face*nsidesq;
                    *(n + 6)=*(n + 7)+2;
                break;
                case 8:              // East corner
                    *nneigh=7;
                    other_face=8+ibp;
                    *(n + 6)=other_face*nsidesq+local_magic2;
                    *n=*(n + 6)-2;
                    *(n + 1)=ipix-1;
                    *(n + 2)=ipix+1;
                    *(n + 3)=ipix+2;
                    other_face=4+ibp;
                    *(n + 5)=other_face*nsidesq;
                    *(n + 4)=*(n + 5)+2;
                break;
            } // End switch
        }

    } // Specail cases
}

} // anonymous namespace
//=========================================================================================
//  C++ interface funtions
//=========================================================================================
Healpix::Healpix(long nside, Healpix::Ordering ord, SkyDir::CoordSystem coordsys)
    : m_nside(nside),
      m_ord(ord),
      m_coordsys(coordsys)
{}

void Healpix::pix2ang(long index, double &theta, double &phi)const
{
    if( !nested()) pix2ang_ring( nside(), index, &theta, &phi);
    else pix2ang_nest( nside(), index, &theta, &phi);
}

void Healpix::ang2pix(double theta, double phi, long &index)const
{
    if(!nested() ) ang2pix_ring( nside(), theta, phi, &index);
    else ang2pix_nest( nside(), theta, phi, &index);
}        


Healpix::Pixel::Pixel(const astro::SkyDir &dir, const Healpix& hp)
: m_healpix(&hp)
{
    // get theta, phi (radians) in appropriate coordinate system
    double theta, phi;
    if( hp.coordsys()==astro::SkyDir::EQUATORIAL){
        theta = M_PI/2- dir.dec()*M_PI/180.;
        phi = dir.ra()*M_PI/180;
    }else{  // galactic
        theta = M_PI/2- dir.b()*M_PI/180.;
        phi = dir.l()*M_PI/180;
    }
    // and look up the pixel number
    m_healpix->ang2pix( theta, phi, m_index);
}

Healpix::Pixel::operator astro::SkyDir ()const
{
    double theta, phi;  
    m_healpix->pix2ang( m_index, theta, phi);
    // convert to ra, dec (or l,b)
    return astro::SkyDir( phi*180/M_PI, (M_PI/2-theta)*180/M_PI, m_healpix->coordsys() );
}

void Healpix::Pixel::neighbors(std::vector<Healpix::Pixel> & p) const
{
    long n[8];
    int nbr_neighbors;
    
    p.clear();
    if (!(this->m_healpix->nested()))
        throw std::runtime_error("Nested ordering required to determine neighbors.");
    
    neighbours_nest(this->m_healpix->m_nside, m_index, n, &nbr_neighbors);
    for (int i = 0; i < nbr_neighbors; ++i)
    {
        p.push_back(Healpix::Pixel(n[i], *(this->m_healpix)));
    }
}

double Healpix::integrate(const astro::SkyFunction& f)const
{
    return std::accumulate(begin(), end(), 0., Integrand(f));
}



