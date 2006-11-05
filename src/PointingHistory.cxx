/** @file PointingHistory.cxx
    @brief implement PointingHistory

    $Header$

    */

#include "astro/PointingHistory.h"
using namespace astro;

#include "tip/IFileSvc.h"
#include "tip/Table.h"

#include <fstream>
#include <sstream>

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
            buf >> lat >> lon >> alt;
            EarthCoordinate earth(lat, lon, alt);

            m_endTime = start + offset;
            m_data[m_endTime] = PointingInfo(position, orientation, earth);
            if( m_startTime<0) m_startTime=m_endTime;

        }
    }
}

const astro::PointingInfo& PointingHistory::operator()(double time)const{

    if( time!=m_selected){
    std::map<double,astro::PointingInfo>::const_iterator iter=m_data.upper_bound(time);
    bool tooEarly( time< (*(m_data.begin())).first );
    bool tooLate( iter==m_data.end() );
    if ( tooEarly || tooLate) {
       std::ostringstream message;
       message << "PointingHistory: Time out of Range!:\n"
               << "Time (" << time 
               << ") out of range of times in the pointing database. "
               << std::endl;
       throw std::runtime_error(message.str());
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
    for ( ; it != scData->end(); ++it) {

        interval["ra_scz"   ].get(raz);
        interval["dec_scz"  ].get(decz);      
        interval["ra_scx"   ].get(rax);
        interval["dec_scx"  ].get(decx);      
        interval["lat_geo"  ].get(lat);
        interval["lon_geo"  ].get(lon);
        interval["start"    ].get(start_time);
        interval["stop"     ].get(stop_time);
        interval["sc_position"].get(sc_pos);
        astro::SkyDir zaxis(raz, decz);
        astro::SkyDir xaxis(rax, decx);
        CLHEP::Hep3Vector position(sc_pos[0]/1e3, sc_pos[1]/1e3, sc_pos[2]/1e3);
        Quaternion orientation(zaxis(), xaxis());

        m_data[start_time] = 
            PointingInfo( position, orientation, EarthCoordinate(lat, lon) );
        if( m_startTime<0) m_startTime = start_time;
    }
    m_endTime = stop_time;
}