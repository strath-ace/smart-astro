/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------ Copyright (C) 2018 University of Strathclyde and Authors ------
-------------------- e-mail: c.greco@strath.ac.uk --------------------
----------------------- Author: Cristian Greco -----------------------
*/


#include "analyticalPlanetsEphemeris.h"

using namespace smartastro;
using namespace smartastro::ephemerides;

/**
 * Default constructor
 *
 * @param Input struct of parameters:
 *
 * Planet possible choiches:
 * MERCURY, VENUS, EARTH, MARS, JUPITER, SATURN, URANUS, NEPTUNE
 *
 * Mu shall be in [km^3/s^2] !
 *
 */
analyticalPlanetsEphemeris::analyticalPlanetsEphemeris(const ephemerisParams* pParams) :
        keplerEphemeris(pParams), m_pParams(static_cast<const analyticalPlanetsParams*>(pParams))
{
    // Sanity check
    if ( m_pParams == nullptr )
        smartastro_throw("Dynamic cast failed: " +
                         "input parameter struct was not pointing to analyticalPlanetsParams");
    
    // Check planet
    initialise();
}


/**
 * Default destructor
 *
 */
analyticalPlanetsEphemeris::~analyticalPlanetsEphemeris()
{
    // nothing to delete
}



/**
 * getKeplerianState: Function that returns object keplerian State for time t jd
 *
 * @param mjd: Julian Date at which the Cartesian State is desired
 * @return Keplerian state at jd
 *         a [km], e [-], i [rad], RAAN [rad], argument of pericenter [rad], true anomaly [rad]
 *
 */
std::vector<double> analyticalPlanetsEphemeris::getKeplerianState( const double& jd )
{
    // Outpiut vector
    std::vector<double> kep(6);

    // Call external routine which uses mjd2000
    double mjd2000;
    smartastro::astrocore::conversion_time::jd2mjd2000(jd,mjd2000);

    analytical_planets::get_orbel(mjd2000,m_ID,kep);

    // Return output value
    return kep;
}


/**
 * Check input planet requested [should be adapted to case-insensitive]
 */
void analyticalPlanetsEphemeris::initialise()
{

    // Check list of planets
    if ( caseInsensitiveStringEqual(m_pParams->planet,"MERCURY") )
        m_ID = 1;
    else if ( caseInsensitiveStringEqual(m_pParams->planet,"VENUS") )
        m_ID = 2;
    else if ( caseInsensitiveStringEqual(m_pParams->planet,"EARTH") )
        m_ID = 3;
    else if ( caseInsensitiveStringEqual(m_pParams->planet,"MARS") )
        m_ID = 4;
    else if ( caseInsensitiveStringEqual(m_pParams->planet,"JUPITER") )
        m_ID = 5;
    else if ( caseInsensitiveStringEqual(m_pParams->planet,"SATURN") )
        m_ID = 6;
    else if ( caseInsensitiveStringEqual(m_pParams->planet,"URANUS") )
        m_ID = 7;
    else if ( caseInsensitiveStringEqual(m_pParams->planet,"NEPTUNE") )
        m_ID = 8;
    else
        smartastro_throw("Input planet name is not in solar system list");

}