/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------ Copyright (C) 2018 University of Strathclyde and Authors ------
-------------------- e-mail: c.greco@strath.ac.uk --------------------
----------------------- Author: Cristian Greco -----------------------
*/


#include "spiceEphemeris.h"

using namespace smartastro;
using namespace smartastro::ephemerides;



/**
 * Default constructor
 *
 * @param Input parameters as defined in spiceEphemerisParams
 */
spiceEphemeris::spiceEphemeris(const ephemerisParams* pParams) :
    base_ephemeris(pParams), m_pParams(static_cast<const spiceEphemerisParams*>(pParams))
{
    // Sanity check
    if ( m_pParams == nullptr )
        smartastro_throw("Dynamic cast failed: " +
                         "input parameter struct was not pointing to spiceEphemerisParams");
}


/**
 * Default destructor
 *
 */
spiceEphemeris::~spiceEphemeris()
{

}



#ifdef __USE_CSPICE

/**
 * getCartesianState: Function that returns object position-velocity for time t
 *
 * @param mjd2000: time at which the Cartesian State is desired in mjd2000 (time system is defined in derived classes)
 * @return Cartesian state at time mjd2000
 *
 */

std::vector<double> spiceEphemeris::getCartesianState( const double& mjd2000 ) const
{
    // Variables
    double jd;
    SpiceDouble et, lt, state[6];
    std::vector<double> car(6);

    // Load kernels
    furnsh_c(spiceKernels::leap.c_str()); // Time leap
    for (unsigned int i = 0 ; i < m_pParams->kernelToLoad.size(); i++)
        furnsh_c(m_pParams->kernelToLoad[i].c_str());

    // Convert mjd2000 in jd
    smartastro::astrocore::conversion_time::mjd20002jd(mjd2000,jd);

    // Convert mjd2000 in seconds after 2000-01-01 noon
    std::string jdString ("jd " + std::to_string((double)jd));
    str2et_c(jdString.c_str(),&et);

    // Get Cartesian state
    spkezr_c( m_pParams->target.c_str(), et,
              m_pParams->referenceFrame.c_str(),
              m_pParams->abberrationCorrection.c_str(),
              m_pParams->observer.c_str(),state,&lt);

    for (unsigned int i = 0 ; i < 6; i++)
        car[i] = state[i];

    return car;
}


#endif