/** @file SkyStat.cxx
@brief Implement SkyStat methods

@author Bruce Lesnick
*/

#include "astro/SkyDir.h"
#include "astro/SkyFunction.h"
#include "astro/HTM.h"
#include "astro/SkyStat.h"
#include <cmath>

using namespace astro;

SkyStat::SkyStat(const astro::SkyFunction & sf, const int level, const astro::HTM * h)
: m_h(h)
, m_sf(sf), m_level(level)
, m_ave(0), m_rms(0), m_min(0), m_max(0), m_rejected(0)
{    
    if (h != 0)
        m_HTM_Created = false;
    else
    {
        m_h = new astro::HTM (m_level);
        m_HTM_Created = true;
    }
    calculate_values();   
}


auto SkyStat::calculate_values() -> int
{
    int count = 0;
    for( HTM::const_iterator it = m_h->begin(m_level); it != m_h->end(m_level); ++it)
    {
        const HTM::Node & n = *it;
        double work = m_sf(n.dir());
        if( fabs(work)> 1e30 ) {
            m_rejected++;
            break;
        }
        ++count;
        m_ave += work;
        m_rms += work * work;
        if (work < m_min || it == m_h->begin(m_level))
            m_min = work;
        if (work > m_max || it == m_h->begin(m_level))
            m_max = work;
    }
    
    m_ave = m_ave / count;
    m_sigma = sqrt(m_rms/count  - (m_ave * m_ave));
    m_rms = sqrt(m_rms / count);
    
    return 0;
}

SkyStat::~SkyStat()
{
    if (m_HTM_Created)
        delete m_h;
}

auto SkyStat::ave()const -> double
{
    return m_ave;
}

auto SkyStat::rms()const -> double
{
    return m_rms;
}
auto SkyStat::sigma()const -> double
{
    return m_sigma;
}
auto SkyStat::min()const -> double
{
    return m_min;
}
auto SkyStat::max()const -> double
{
    return m_max;
}
