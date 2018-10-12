/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------ Copyright (C) 2018 University of Strathclyde and Authors ------
-------------------- e-mail: c.greco@strath.ac.uk --------------------
----------------------- Author: Cristian Greco -----------------------
*/


#include "rangeSensor.h"

using namespace smartastro::sensors;


/**
 * Default constructor
 *
 */
rangeSensor::rangeSensor(const base_sensor::sensorParams* pParams) :
        base_sensor(pParams)
{
    // Nothing to initialise
}


/**
 * Default destructor
 *
 */
rangeSensor::~rangeSensor()
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
std::vector<double> rangeSensor::getPerfectMeasurement(const double &t)
{
    // Return value
    double range = -1.0;

    if (m_absoluteEphemeris) {

        // Compute state of observer sensor and target at time t
        std::vector<double> sensorPosition = m_pParams->getSensorEphemeris(t);
        std::vector<double> targetPosition = m_pParams->getTargetEphemeris(t);

        // Sanity check on dimensions
        if (sensorPosition.size() != 3 && sensorPosition.size() != 6)
            smartastro_throw("Sensor state size = " + std::to_string(sensorPosition.size()) + " differs from pure position or cartesian state");
        if (targetPosition.size() != 3 && targetPosition.size() != 6)
            smartastro_throw("Target state size = " + std::to_string(targetPosition.size()) + " differs from pure position or cartesian state");

        // Resize to be sure to keep only position
        sensorPosition.resize(3);
        targetPosition.resize(3);

        // Compute range
        range = smartastro::sensors::computeRange(sensorPosition, targetPosition);
    }
    else if (m_relativeEphemeris)
    {
        // Compute relative state of observer-target at time t
        std::vector<double> relativePosition = m_pParams->getSensorTargetRelativeEphemeris(t);

        // Sanity check on dimensions
        if (relativePosition.size() != 3 && relativePosition.size() != 6)
            smartastro_throw("Relative state size differs from pure position or cartesian state");

        // Resize to be sure to keep only position
        relativePosition.resize(3);

        // Compute range
        range = smartastro::sensors::computeRange(relativePosition);
    }
    else
        smartastro_throw("Range requested but function pointers not assigned");

    // Sanity check on sign
    if(range<0.0)
        smartastro_throw("Range should always be a nonnegative quantity");

    return std::vector<double>(1,range);
}



/**
 * Compute range given two position vectors in the same reference frame
 *
 * @return Range between two position vectors
 */
double smartastro::sensors::computeRange ( const std::vector<double>& pos1,
                                           const std::vector<double>& pos2)
{
    // Sanity check on dimensions
    if ( pos1.size()!=3 )
        smartastro_throw("First position size differs from 3");
    if ( pos2.size()!=3 )
        smartastro_throw("Second position size differs from 3");

    // Compute range
    double range = std::sqrt(
            std::pow(pos1[0]-pos2[0],2.0) + std::pow(pos1[1]-pos2[1],2.0) + std::pow(pos1[2]-pos2[2],2.0)
            );

    return range;
}

/**
 * Compute range given relative position vector
 *
 * @return Range of relative position vector
 */
double smartastro::sensors::computeRange ( const std::vector<double>& relPos )
{
    // Sanity check on dimensions
    if ( relPos.size()!=3 )
        smartastro_throw("Relative position size differs from 3");

    // Compute range
    double range = std::sqrt(
            std::pow(relPos[0],2.0) + std::pow(relPos[1],2.0) + std::pow(relPos[2],2.0)
    );

    return range;
}