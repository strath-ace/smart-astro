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
base_observation::base_observation(const observationParams* pParams) :
    m_pParams(pParams), m_noisyObservation(false)
{
    // Initialise flags
    this->initialise();
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
 * getObservation: Function that returns noisy measurements at time t
 *
 * @param t: time at which the measurement(time system is defined in derived classes)
 * @return Measurement vector
 *
 */
std::vector<double> base_observation::getObservation( const double& t )
{
    // Get perfect measurement
    std::vector<double> measurement = getPerfectObservation(t);

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


// Set noisy measurements
void base_observation::setNoisyObservation( const bool noisyObservation )
{
    // Assign value
    m_noisyObservation = noisyObservation;

    // Check if generating function exists
    if (m_noisyObservation && !m_pParams->generateMeasurementNoise)
        smartastro_throw("Noisy measurements requested but function to generate noise not assigned");
}




/**
 * Private routine called by constructor
 */

// Check whether relative or absolute positions shall be used
void base_observation::initialise()
{
    // Absolute state flag
    if(m_pParams->getSensorEphemeris && m_pParams->getTargetEphemeris)
        m_absoluteEphemeris = true;
    else
        m_absoluteEphemeris = false;

    // Relative state flag
    if(m_pParams->getSensorTargetRelativeEphemeris)
        m_relativeEphemeris = true;
    else
        m_relativeEphemeris = false;
}



// Get noise value and save it
std::vector<double> base_observation::getNoiseSample()
{
    // Compute noise sample
    std::vector<double> noiseSample = m_pParams->generateMeasurementNoise();

    // Save current value
    m_currentNoiseSample = noiseSample;

    return m_currentNoiseSample;
}