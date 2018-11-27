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
using namespace std;

/**
 * Default constructor
 *
 * @param Input parameters as defined in integratedEphemerisParams
 */
integratedEphemeris::integratedEphemeris(const ephemerisParams* pParams) :
        base_ephemeris(pParams), m_pParams(static_cast<const integratedEphemerisParams*>(pParams)),
        m_histT(vector<double>(1,m_pParams->ti)), m_histX(vector<vector<double>>(1,m_pParams->xi))
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
vector<double> integratedEphemeris::getCartesianState( const double& t )
{
    // Initial time and state vector
    double ti;
    vector<double> xi;

    // Check which is closer
    unsigned int nT = m_histT.size();
    int cc = 0;

    // First, last and between times
    if (t<m_histT[0])
    {
        ti = m_histT[0];
        xi = m_histX[0];
    }
    else if (t>m_histT.back())
    {
        cc = nT;
        ti = m_histT.back();
        xi = m_histX.back();
    }
    else if (abs(t-m_histT.back())<1.0e-8)
        return m_histX.back();
    else
    {
        bool closerFound = false;
        while (cc<nT-1 && !closerFound)
        {
            if (abs(t-m_histT[cc])<1.0e-8)
                return m_histX[cc];

            if (t>m_histT[cc] && t<m_histT[cc+1])
            {
                ti = m_histT[cc];
                xi = m_histX[cc];
                closerFound = true;
            }
            cc++;
        }
    }

    // Integrate to final state x
    vector<double> x(xi);
    int nsteps = ceil(abs((t-ti)/m_pParams->hstep));
    m_pParams->integrate(ti,t,nsteps,xi,x);

    // Set last integrated time and state in history
    m_histT.insert(m_histT.begin()+cc,t);
    m_histX.insert(m_histX.begin()+cc,x);

    // Return value
    return x;
}



/**
 * Get time history
 */
std::vector<double> integratedEphemeris::getTimeHistory() const
{
    return m_histT;
}

/**
 * Get state history
 */
std::vector<std::vector<double>> integratedEphemeris::getStateHistory() const
{
    return m_histX;
}

/**
 * Set time history
 */
void integratedEphemeris::setTimeHistory(const std::vector<double> &THist)
{
    m_histT = THist;
}

/**
 * Set state history
 */
void integratedEphemeris::setStateHistory(const std::vector<std::vector<double>> &XHist)
{
    m_histX = XHist;
}

/**
 * Reset current time and state
 */
void integratedEphemeris::resetCurrentTimeState()
{
    m_histX.clear();
    m_histT.clear();
}