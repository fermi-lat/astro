/** @file PointingHistory.cxx
    @brief implement PointingHistory

    $Header: /nfs/slac/g/glast/ground/cvs/astro/src/PointingHistory.cxx,v 1.11 2008/06/29 11:06:52 burnett Exp $

    */

#include "astro/PointingHistory.h"
using namespace astro;

#include "tip/IFileSvc.h"
#include "tip/Table.h"

#include <fstream>
#include <sstream>

namespace {
    static double time_tol(10); // seconds allow beyond the end
}

PointingHistory::PointingHistory(const std::string& filename, double offset)
: m_selected(-1)
, m_startTime(-1)
{
    if( haveFitsFile(filename) ){
        readFitsData(filename);
    }else{
        readTextData(filename, offset);
    }
}
void PointingHistory::readTextData(std::string filename, double offset)
{
    std::ifstream input_file;
    input_file.open(filename.c_str());

    if(!input_file.is_open()) {
        std::cerr << "ERROR:  Unable to open:  " << filename << std::endl;
        throw std::invalid_argument("PointingHistory: could not open pointing history file");
    }else{
        //initialize the key structure:
        while (!input_file.eof()){
            std::string line; std::getline(input_file, line);
            double start; std::stringstream buf(line); buf >>start;
            double x,y,z;
            buf >> x >> y >> z;
            Hep3Vector position(x,y,z);

            double raz, decz, rax, decx;
            buf >> raz >> decz >> rax >> decx;
            Quaternion orientation(SkyDir(raz,decz)(), SkyDir(rax,decx)());

            double razenith, deczenith;
            buf >> razenith >> deczenith; // ignore since redundant with position
            double lat, lon, alt;
            buf >> lon >> lat >> alt;     // these are ignored, too.

            EarthCoordinate earth(position, start);

            m_endTime = start + offset;
            m_data[m_endTime] = PointingInfo(position, orientation, earth);
            if( m_startTime<0) m_startTime=m_endTime;

        }
    }
}

const astro::PointingInfo& PointingHistory::operator()(double time)const throw(TimeRangeError){

    if( time!=m_selected){

        std::map<double,astro::PointingInfo>::const_iterator iter=m_data.upper_bound(time);
        bool tooEarly( time< (*(m_data.begin())).first );
        bool tooLate( iter==m_data.end() && time > m_endTime);

        if(tooEarly && time==0) {
            // special case: just return starting value
            time = (*(m_data.begin())).first;
            tooEarly = false;
        }

        if ( tooEarly || tooLate) {
            std::ostringstream message;
            message << "PointingHistory: Time out of Range!:\n"
                << "Time (" << time 
                << ") out of range of times in the pointing database: ("
                << m_startTime <<", "<< m_endTime <<")"
                << std::endl;
            throw TimeRangeError(message.str());
        }
        const PointingInfo & h2 = iter->second;
        double time2( iter->first);
        --iter;
        double time1(iter->first);
        const PointingInfo & h1 =iter->second;
        double prop( (time-time1)/(time2-time1) );

        m_currentPoint = h1.interpolate(h2, prop);

    }
    return m_currentPoint;
}

bool PointingHistory::haveFitsFile(std::string filename) const {
   std::ifstream file(filename.c_str());
   std::string line;
   std::getline(file, line, '\n');
   file.close();
   if (line.find("SIMPLE") == 0) {
      return true;
   }
   return false;
}

void PointingHistory::readFitsData(std::string filename) {

    std::string sctable("SC_DATA"); //?
    const tip::Table * scData = 
        tip::IFileSvc::instance().readTable(filename, sctable);
    double start_time, stop_time, raz, decz, rax, decx, lat,lon;
    std::vector<double> sc_pos(3);
    tip::Table::ConstIterator it = scData->begin();
    tip::ConstTableRecord & interval = *it;
    double maxlatdiff(0);
    bool have_quaternion(false);
    double qsj_1, qsj_2, qsj_3, qsj_4;
#if 0 // back this out
    try{
        interval["QSJ_1"].get(qsj_1);
        have_quaternion= true;
    }catch(const std::exception&){
    }
#endif
    for ( ; it != scData->end(); ++it) {

        Quaternion orientation;
        if( have_quaternion ) {
            interval["QSJ_1"].get(qsj_1);
            interval["QSJ_2"].get(qsj_2);
            interval["QSJ_3"].get(qsj_3);
            interval["QSJ_4"].get(qsj_4);
            orientation = Quaternion(Hep3Vector(qsj_1,qsj_2,qsj_3), qsj_4);

        }else{
            interval["ra_scz"   ].get(raz);
            interval["dec_scz"  ].get(decz);      
            interval["ra_scx"   ].get(rax);
            interval["dec_scx"  ].get(decx);  
            astro::SkyDir zaxis(raz, decz);
            astro::SkyDir xaxis(rax, decx);
            orientation = Quaternion(zaxis(), xaxis());
        }
        interval["lat_geo"  ].get(lat);
        interval["lon_geo"  ].get(lon);
        interval["start"    ].get(start_time);
        interval["stop"     ].get(stop_time);
        interval["sc_position"].get(sc_pos);
        CLHEP::Hep3Vector position(sc_pos[0]/1e3, sc_pos[1]/1e3, sc_pos[2]/1e3);
        EarthCoordinate earthpos(position, start_time);
    
        // check consistency of latitude, longitude: EarthCoordinate computes from the MET and position
        double check_lat(earthpos.latitude()-lat), check_lon(earthpos.longitude()-lon);
        if( check_lat>maxlatdiff) maxlatdiff = check_lat;

        static double lat_tol(0.5), lon_tol(0.5); // was 0.2, 0.01

        if( fabs(check_lat)>lat_tol || fabs(check_lon)>lon_tol && abs(lon)<179 ){
            std::stringstream error; 
            error << "PointingHistory::readFitsData: apparent inconsistency for Earth position, time=" 
                << start_time
                << ", lat, lon diff: " << check_lat << ", " << check_lon;
            std::cerr << error.str() << std::endl;
#if 0 // disable for now
            throw std::runtime_error(error.str());
#endif
        }

        m_data[start_time] = 
            PointingInfo( position, orientation, earthpos);
        if( m_startTime<0) m_startTime = start_time;
    }
    m_endTime = stop_time+time_tol;
    // an extra entry to allow query for the last interval +slop
    m_data[m_endTime] = m_data[start_time];
}
