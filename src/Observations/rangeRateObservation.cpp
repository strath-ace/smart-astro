/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------ Copyright (C) 2018 University of Strathclyde and Authors ------
-------------------- e-mail: c.greco@strath.ac.uk --------------------
----------------------- Author: Cristian Greco -----------------------
*/


#include "rangeRateObservation.h"

using namespace smartastro::observations;


/**
 * Default constructor
 *
 */
rangeRateObservation::rangeRateObservation(const base_observation::observationParams* pParams) :
        base_observation(pParams,RANGERATE)
{
    // Nothing to initialise
}


/**
 * Default destructor
 *
 */
rangeRateObservation::~rangeRateObservation()
{
    // Nothing to delete
}


/**
 * getObservation: Function that returns measurements at time t
 *
 * @param t: time at which the measurement(time system is defined in derived classes)
 * @return Measurement vector
 *
 */
std::vector<double> rangeRateObservation::getPerfectObservation(const double &t)
{
    // Return value
    double rangeRate = 0.0;

    if (m_absoluteEphemeris) {

        // Compute state of observer sensor and target at time t
        std::vector<double> sensorState = m_pParams->getSensorEphemeris(t);
        std::vector<double> targetState = m_pParams->getTargetEphemeris(t);

        // Compute rangeRate [sanity check on dimensions inside]
        rangeRate = smartastro::observations::computeRangeRate(sensorState, targetState);
    }
    else if (m_relativeEphemeris)
    {
        // Compute state of observer sensor and target at time t
        std::vector<double> relativeState = m_pParams->getSensorTargetRelativeEphemeris(t);

        // Compute rangeRate [sanity check on dimensions inside]
        rangeRate = smartastro::observations::computeRangeRate(relativeState);
    }
    else
        smartastro_throw("RangeRate requested but function pointers not assigned");

    return std::vector<double>(1,rangeRate);
}



/**
 * Compute rangeRate given two state vectors in the same reference frame
 *
 * @return RangeRate between two state vectors
 */
double smartastro::observations::computeRangeRate ( const std::vector<double>& state1,
                                               const std::vector<double>& state2)
{
    // Sanity check on dimensions
    if ( state1.size()!=6 )
        smartastro_throw("First state size differs from 6");
    if ( state2.size()!=6 )
        smartastro_throw("Second state size differs from 6");

    // Compute rangeRate from relative state routine
    std::vector<double> relState (6);
    for (unsigned int i = 0 ; i < 6; i++)
        relState[i] = state2[i]-state1[i];
    double rangeRate = computeRangeRate(relState);

    return rangeRate;
}

/**
 * Compute rangeRate given relative state vector
 *
 * @return RangeRate
 */
double smartastro::observations::computeRangeRate (const std::vector<double>& relState)
{
    // Sanity check on dimensions
    if ( relState.size()!=6 )
        smartastro_throw("Relative state size differs from 6");

    // Compute range
    double range = std::sqrt(
            std::pow(relState[0],2.0) + std::pow(relState[1],2.0) + std::pow(relState[2],2.0)
    );

    // Range rate as projection of velocity on scaled position vector
    double rangeRate = 0.0;
    for ( unsigned int i = 0 ; i < 3; i++ )
        rangeRate += relState[i+3]*relState[i];
    rangeRate /= range;

    return rangeRate;
}
