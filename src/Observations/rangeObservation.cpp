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
std::vector<double> rangeObservation::getPerfectObservation( const std::vector<double>& targetState,
                                                             const std::vector<double>& sensorState )
{
    // Return value
    double range ;

    // Compute range
    range = smartastro::observations::computeRange(targetState, sensorState);

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
double smartastro::observations::computeRange ( const std::vector<double>& targetPos,
                                                const std::vector<double>& sensorPos)
{
    // Sanity check on dimensions
    if ( targetPos.size()!=3 && targetPos.size()!=6 )
        smartastro_throw("Target state size differs from 3 or 6");
    if ( !sensorPos.empty() && sensorPos.size()!=3 && sensorPos.size()!=6 )
        smartastro_throw("Target state size differs from 3 or 6");

    // Compute range
    std::vector<double> relPos (targetPos.begin(),targetPos.begin()+3);
    if (!sensorPos.empty())
        for (unsigned int i = 0 ; i < 3; i++)
            relPos[i] -= sensorPos[i];

    double range = computeRange(relPos);

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