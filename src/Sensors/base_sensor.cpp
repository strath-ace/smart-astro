/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------ Copyright (C) 2018 University of Strathclyde and Authors ------
-------------------- e-mail: c.greco@strath.ac.uk --------------------
----------------------- Author: Cristian Greco -----------------------
*/


#include "base_sensor.h"

using namespace smartastro::sensors;


/**
 * Default constructor
 *
 * @param Parameters structure
 */
base_sensor::base_sensor(const baseSensorParams* pParams) :
        m_pParams(pParams), m_currTime(std::nan("NaN"))
{
    // Check input
    checkInput();
}


/**
 * Default destructor
 *
 */
base_sensor::~base_sensor()
{
    // Nothing to delete
}


/**
 * getMeasurements: Function that returns measurement of target state targetState
 *
 */
std::vector<OutputSensor> base_sensor::getMeasurements( const std::vector<double>& targetState,
                                                                 const double& jd )
{
    // Update sensor state
    if (!std::isnan(jd))
        update(jd);

    // Compute observations
    unsigned int nObs = m_pParams->vGetObservations.size();

    // Allocate output
    std::vector<OutputSensor> outputObs (nObs);

    // Assign output
    for (unsigned int i = 0 ; i < nObs; i++)
    {
        outputObs[i].sensor = m_pParams->name;
        outputObs[i].time = m_currTime;
        outputObs[i].number = i;
        outputObs[i].value = m_pParams->vGetObservations[i](m_currSensorState,targetState);
    }

    // Return output
    return outputObs;
}


/**
 * Update class member at time jd
 *
 * @param jd = time in Julian Dates
 */
void base_sensor::update(const double& jd)
{
    // If nan throw an error
    if (std::isnan(jd))
        smartastro_throw("Input time is NaN");

    // Current time
    m_currTime = jd;

    // Current state
    m_currSensorState.clear();
    m_currSensorState = m_pParams->getSensorEphemeris(m_currTime);
}



/**
 * Private routines
 */


/**
 * Check routine that ensures good input
 */
void base_sensor::checkInput()
{
    if(!m_pParams->getSensorEphemeris)
        smartastro_throw("Sensor ephemeris not set");

    if (m_pParams->vGetObservations.empty())
        smartastro_throw("Empty vector of observations");
}

