/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------ Copyright (C) 2018 University of Strathclyde and Authors ------
-------------------- e-mail: c.greco@strath.ac.uk --------------------
----------------------- Author: Cristian Greco -----------------------
*/


#include "keplerEphemeris.h"

using namespace smartastro::ephemerides;


/**
 * Default constructor
 *
 * @param referenceFrame: name of reference frame to which ephemerides are referred
 */
keplerEphemeris::keplerEphemeris(const ephemerisParams* pParams) :
        base_ephemeris(pParams), m_pParams(static_cast<const keplerEphemerisParams*>(pParams))
{
    // Sanity check
    if ( m_pParams == nullptr )
        smartastro_throw("Dynamic cast failed: " +
                         "input parameter struct was not pointing to keplerEphemerisParams");
}


/**
 * Default destructor
 *
 */
keplerEphemeris::~keplerEphemeris()
{
    // nothing to destroy
}


/**
 * getCartesianState: Function that returns object position-velocity for time t
 *
 * Virtual function to be implemented in derived class
 *
 * @param t: time at which the Cartesian State is desired (time system is defined in derived classes)
 * @return Cartesian state at time t
 *
 */
std::vector<double> keplerEphemeris::getCartesianState( const double& t ) const
{
    // Sanity check
    if(!m_pParams->centralBodyGravitationalParameter)
        smartastro_throw("Centra body gravitational parameter not set");

    // Get Keplerian state
    std::vector<double> kep = getKeplerianState(t);

    // Convert to Cartesian state
    std::vector<double> car(6);
    smartastro::astrocore::conversion_coordinates::kep2car(kep,*m_pParams->centralBodyGravitationalParameter,car);

    // Return Cartesian state
    return car;
}