/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------ Copyright (C) 2018 University of Strathclyde and Authors ------
-------------------- e-mail: c.greco@strath.ac.uk --------------------
----------------------- Author: Cristian Greco -----------------------
*/


#include "groundStation.h"

using namespace smartastro;
using namespace smartastro::sensors;


/**
 * Default constructor
 *
 * @param Parameters structure
 */
groundStation::groundStation(const groundStationParams* pParams):
    m_pParams(pParams)
{
    // Check input
    checkInput();
}


/**
 * Default destructor
 *
 */
groundStation::~groundStation()
{
    // Nothing to destroy
}


/**
 * getMeasurements: Function that returns measurement of target for time t
 *
 * @param t: time at which the measurements are desired
 * @return Measurements at Julian Date jd
 *
 */
std::vector<observations::OutputMeasurementValue> groundStation::getMeasurements( const double& jd )
{
    // Number of measurements
    unsigned int nObs = m_pParams->observations.size();

    // Output
    std::vector<observations::OutputMeasurementValue> outputObs(nObs);

    // Get measurements
    for ( unsigned int i = 0 ; i < nObs; i++ )
    {
        // Station of measurement
        outputObs[i].station = m_pParams->name;

        // Type of measurement
        outputObs[i].type  = m_pParams->observations[i]->getObservationType();

        // Time of measurement
        outputObs[i].time  = jd;

        // Value of observation
        outputObs[i].value = m_pParams->observations[i]->getObservation(jd);
    }

    return outputObs;
}



/**
 * Check routine that ensures good input
 */
void groundStation::checkInput()
{
//    if (!(m_pParams->getStationEphemeris))
//        smartastro_throw("Station ephemerides not set");
//
//    if (!(m_pParams->getTargetEphemeris))
//        smartastro_throw("Target ephemerides not set");
//
//    if (m_pParams->obsTypes.empty())
//        smartastro_throw("No observation types requested");

    if (m_pParams->observations.empty())
        smartastro_throw("No observation present");

}