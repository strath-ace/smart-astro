/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------ Copyright (C) 2018 University of Strathclyde and Authors ------
-------------------- e-mail: c.greco@strath.ac.uk --------------------
----------------------- Author: Cristian Greco -----------------------
*/


#include "Observations/azimuthElevationObservation.h"

using namespace smartastro::observations;


/**
 * Default constructor
 *
 */
azimuthElevationObservation::azimuthElevationObservation(const base_observation::observationParams* pParams) :
        base_observation(pParams,AZIMUTHELEVATION)
{
    // Nothing to initialise
}


/**
 * Default destructor
 *
 */
azimuthElevationObservation::~azimuthElevationObservation()
{
    // Nothing to delete
}



/**
 * getObservation: Function that returns measurements
 *
 * @return Measurement vector
 *
 */
std::vector<double> azimuthElevationObservation::getPerfectObservation( const std::vector<double>& targetStateTOPO,
                                                                        const std::vector<double>& sensorStateTOPO )
{
    // Return value
    std::vector<double> azimuthElevation (2);

    // Compute azimuthElevation [sanity check on dimensions inside]
    azimuthElevation = smartastro::observations::computeAzimuthElevation(targetStateTOPO, sensorStateTOPO);

    return azimuthElevation;
}



/**
 * Compute azimuthElevation given two state vectors in the same reference frame
 *
 * @return AzimuthElevation between two state vectors
 */
std::vector<double> smartastro::observations::computeAzimuthElevation ( const std::vector<double>& targetTOPO,
                                                                        const std::vector<double>& sensorTOPO)
{
    // Sanity check on dimensions
    if ( targetTOPO.size()!=6 )
        smartastro_throw("First state size differs from 6");
    if ( !sensorTOPO.empty()&&sensorTOPO.size()!=6 )
        smartastro_throw("Second state size differs from 6");

    // Compute azimuthElevation from relative state routine
    std::vector<double> relTOPOCartState (targetTOPO);
    if (!sensorTOPO.empty())
        for (unsigned int i = 0 ; i < 6; i++)
            relTOPOCartState[i] -= sensorTOPO[i];

    return computeAzimuthElevation(relTOPOCartState);
}


/**
 * Compute azimuthElevation given relative state vector
 *
 * @return AzimuthElevation
 */
std::vector<double> smartastro::observations::computeAzimuthElevation (const std::vector<double>& relTOPOCartState)
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
