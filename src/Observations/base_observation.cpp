/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------ Copyright (C) 2018 University of Strathclyde and Authors ------
-------------------- e-mail: c.greco@strath.ac.uk --------------------
----------------------- Author: Cristian Greco -----------------------
*/


#include "base_observation.h"


using namespace smartastro::observations;


/**
 * Default constructor
 *
 */
base_observation::base_observation(const observationParams* pParams,
                                   const observationTypes&  obsType) :
        m_pParams(pParams), m_observationType(obsType)
{
    // Check input
    checkInput();
}


/**
 * Default destructor
 *
 */
base_observation::~base_observation()
{
    // Nothing to delete
}


/**
 * getObservation: Function that returns noisy measurements
 *
 * @return Measurement vector
 *
 */
std::vector<double> base_observation::getNoisyObservation( const std::vector<double>& sensorState,
                                                           const std::vector<double>& targetState )
{
    // Get perfect measurement
    std::vector<double> measurement = getPerfectObservation(sensorState,targetState);

    // Add noisy if requested
    if (m_noisyObservation)
    {
        // Get noise sample
        std::vector<double> noise = getNoiseSample();

        // Check dimension
        if (noise.size()!=measurement.size())
            smartastro_throw("Measurement and noise sizes do not coincide");

        // Add noise
        for ( unsigned int i = 0 ; i < measurement.size(); i++ )
            measurement[i] += noise[i];
    }

    return measurement;
}


/**
 * Functions to set member attributes
 */


// Get observation type
observationTypes base_observation::getObservationType () const
{
    return m_observationType;
}


/**
 * Private routine called by constructor
 */

// Check whether relative or absolute positions shall be used
void base_observation::checkInput()
{
    if(m_pParams->generateMeasurementNoise)
        m_noisyObservation = true;
    else
        m_noisyObservation = false;
}



// Get noise value and save it
std::vector<double> base_observation::getNoiseSample()
{
    if (!m_noisyObservation)
        smartastro_throw("Noise sample requested but noise generator not set");

    // Compute noise sample
    std::vector<double> noiseSample = m_pParams->generateMeasurementNoise();

    // Save current value
    m_currentNoiseSample = noiseSample;

    return m_currentNoiseSample;
}