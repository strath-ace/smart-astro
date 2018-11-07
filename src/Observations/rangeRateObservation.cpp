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
 * getObservation: Function that returns measurements
 *
 * @return Measurement vector
 *
 */
std::vector<double> rangeRateObservation::getPerfectObservation( const std::vector<double>& targetState,
                                                                 const std::vector<double>& sensorState )
{
    // Return value
    double rangeRate ;

    // Compute rangeRate [sanity check on dimensions inside]
    rangeRate = smartastro::observations::computeRangeRate(targetState, sensorState);

    return std::vector<double>(1,rangeRate);
}



/**
 * Compute rangeRate given two state vectors in the same reference frame
 *
 * @return RangeRate between two state vectors
 */
double smartastro::observations::computeRangeRate ( const std::vector<double>& targetState,
                                                    const std::vector<double>& sensorState)
{
    // Sanity check on dimensions
    if ( targetState.size()!=6 )
        smartastro_throw("First state size differs from 6");
    if ( !sensorState.empty()&&sensorState.size()!=6 )
        smartastro_throw("Second state size differs from 6");

    // Compute rangeRate from relative state routine
    std::vector<double> relState (targetState);

    if (!sensorState.empty())
        for (unsigned int i = 0 ; i < 6; i++)
            relState[i] -= sensorState[i];

    return computeRangeRate(relState);
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
