/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------ Copyright (C) 2018 University of Strathclyde and Authors ------
-------------------- e-mail: c.greco@strath.ac.uk --------------------
----------------------- Author: Cristian Greco -----------------------
*/


#include "rangeObservation.h"

using namespace smartastro::observations;


/**
 * Default constructor
 *
 */
rangeObservation::rangeObservation(const base_observation::observationParams* pParams) :
        base_observation(pParams,RANGE)
{
    // Nothing to initialise
}


/**
 * Default destructor
 *
 */
rangeObservation::~rangeObservation()
{
    // Nothing to delete
}


/**
 * getObservation: Function that returns measurements
 *
 * @return Measurement vector
 *
 */
std::vector<double> rangeObservation::getPerfectObservation( const std::vector<double>& sensorState,
                                                             const std::vector<double>& targetState )
{
    // Return value
    double range ;

    // Compute range
    range = smartastro::observations::computeRange(sensorState, targetState);

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
double smartastro::observations::computeRange ( const std::vector<double>& state1,
                                                const std::vector<double>& state2)
{
    // Sanity check on dimensions
    if ( state1.size()!=3 && state1.size()!=6 )
        smartastro_throw("First state size differs from 3 or 6");
    if ( state2.size()!=3 && state2.size()!=6 )
        smartastro_throw("Second state size differs from 3 or 6");

    // Compute range
    double range = std::sqrt(
            std::pow(state1[0]-state2[0],2.0) + std::pow(state1[1]-state2[1],2.0) + std::pow(state1[2]-state2[2],2.0)
    );

    return range;
}

/**
 * Compute range given relative position vector
 *
 * @return Range of relative position vector
 */
double smartastro::observations::computeRange ( const std::vector<double>& relState )
{
    // Sanity check on dimensions
    if ( relState.size()!=3 && relState.size()!=6 )
        smartastro_throw("Relative position size differs from 3 or 6");

    // Compute range
    double range = std::sqrt(
            std::pow(relState[0],2.0) + std::pow(relState[1],2.0) + std::pow(relState[2],2.0)
    );

    return range;
}