/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------ Copyright (C) 2018 University of Strathclyde and Authors ------
-------------------- e-mail: c.greco@strath.ac.uk --------------------
----------------------- Author: Cristian Greco -----------------------
*/


#include "Ephemerides/integratedEphemeris.h"

using namespace smartastro::ephemerides;

/**
 * Default constructor
 *
 * @param Input parameters as defined in integratedEphemerisParams
 */
integratedEphemeris::integratedEphemeris(const ephemerisParams* pParams) :
        base_ephemeris(pParams), m_pParams(static_cast<const integratedEphemerisParams*>(pParams)), m_currT(-1.0e10)
{
    // Sanity check
    if ( m_pParams == nullptr )
        smartastro_throw("Dynamic cast failed: " +
                         "input parameter struct was not pointing to integratedEphemerisParams");

    if (!m_pParams->integrate)
        smartastro_throw("Integrator not set");
}


/**
 * Default destructor
 *
 */
integratedEphemeris::~integratedEphemeris()
{
    // Nothing to delete
}


/**
 * getCartesianState: Function that returns object position-velocity for time t by numerical integration
 *
 * @param t: time at which the Cartesian State is desired
 * @return Cartesian state at time t
 *
 */
std::vector<double> integratedEphemeris::getCartesianState( const double& t )
{
    // Initial time and state vector
    double ti;
    std::vector<double> xi;

    // Check which is closer
    if ( std::abs(t-m_pParams->ti)<=std::abs(t-m_currT) )
    {
        ti = m_pParams->ti;
        xi = m_pParams->xi;
    }
    else
    {
        ti = m_currT;
        xi = m_currX;
    }

    // Integrate to final state x
    std::vector<double> x(xi);
    m_pParams->integrate(ti,t,m_pParams->nsteps,xi,x);

    // Set current time and state to last integrated one
    m_currT = t;
    m_currX.clear();
    m_currX = x;

    // Return value
    return x;
}