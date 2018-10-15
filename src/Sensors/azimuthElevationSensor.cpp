/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------ Copyright (C) 2018 University of Strathclyde and Authors ------
-------------------- e-mail: c.greco@strath.ac.uk --------------------
----------------------- Author: Cristian Greco -----------------------
*/


#include "Sensors/azimuthElevationSensor.h"

using namespace smartastro::sensors;


/**
 * Default constructor
 *
 */
azimuthElevationSensor::azimuthElevationSensor(const base_sensor::sensorParams* pParams) :
        base_sensor(pParams)
{
    // Nothing to initialise
}


/**
 * Default destructor
 *
 */
azimuthElevationSensor::~azimuthElevationSensor()
{
    // Nothing to delete
}


/**
 * getMeasurement: Function that returns measurements at time t
 *
 * @param t: time at which the measurement(time system is defined in derived classes)
 * @return Measurement vector
 *
 */
std::vector<double> azimuthElevationSensor::getPerfectMeasurement(const double &t)
{
    // Return value
    std::vector<double> azimuthElevation (2);

    if (m_absoluteEphemeris) {

        // Compute state of observer sensor and target at time t
        std::vector<double> sensorState = m_pParams->getSensorEphemeris(t);
        std::vector<double> targetState = m_pParams->getTargetEphemeris(t);

        // Compute azimuthElevation [sanity check on dimensions inside]
        azimuthElevation = smartastro::sensors::computeAzimuthElevation(sensorState, targetState);
    }
    else if (m_relativeEphemeris)
    {
        // Compute state of observer sensor and target at time t
        std::vector<double> relativeState = m_pParams->getSensorTargetRelativeEphemeris(t);

        // Compute azimuthElevation [sanity check on dimensions inside]
        azimuthElevation = smartastro::sensors::computeAzimuthElevation(relativeState);
    }
    else
        smartastro_throw("AzimuthElevation requested but function pointers not assigned");

    return azimuthElevation;
}



/**
 * Compute azimuthElevation given two state vectors in the same reference frame
 *
 * @return AzimuthElevation between two state vectors
 */
std::vector<double> smartastro::sensors::computeAzimuthElevation ( const std::vector<double>& stateTOPO1,
                                                                   const std::vector<double>& stateTOPO2)
{
    // Sanity check on dimensions
    if ( stateTOPO1.size()!=6 )
        smartastro_throw("First state size differs from 6");
    if ( stateTOPO2.size()!=6 )
        smartastro_throw("Second state size differs from 6");

    // Compute azimuthElevation from relative state routine
    std::vector<double> relTOPOCartState (6);
    for (unsigned int i = 0 ; i < 6; i++)
        relTOPOCartState[i] = stateTOPO2[i]-stateTOPO1[i];

    return computeAzimuthElevation(relTOPOCartState);
}


/**
 * Compute azimuthElevation given relative state vector
 *
 * @return AzimuthElevation
 */
std::vector<double> smartastro::sensors::computeAzimuthElevation (const std::vector<double>& relTOPOCartState)
{
    // Sanity check on dimensions
    if ( relTOPOCartState.size()!=6 )
        smartastro_throw("Relative state size differs from 6");

    // Initialise quantities
    std::vector<double> relTOPOSpherState (6);
    std::vector<double> azimuthElevation (2);

    // Convert topocentric reference frame to spherical coordinated
    smartastro::astrocore::conversion_coordinates::car2spher(relTOPOCartState,relTOPOSpherState);

    // Azimuth [clockwise from north direction]
    azimuthElevation[0] = - relTOPOSpherState[1];
    if (azimuthElevation[0]<0.0)
        azimuthElevation[0] += 2.0*M_PI;
    if (azimuthElevation[0]>2.0*M_PI)
        azimuthElevation[0] -= 2.0*M_PI;

    // Elevation [positive if +z-axis of topocentric frame]
    azimuthElevation[1] = relTOPOSpherState[2];

    return azimuthElevation;
}
