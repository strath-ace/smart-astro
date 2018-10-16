/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------ Copyright (C) 2018 University of Strathclyde and Authors ------
-------------------- e-mail: c.greco@strath.ac.uk --------------------
----------------------- Author: Cristian Greco -----------------------
*/


#ifndef SMART_ASTRO_BASE_SENSOR_H
#define SMART_ASTRO_BASE_SENSOR_H

#include <vector>
#include <functional>
#include <cmath>

#include "../exception.h"

namespace smartastro
{
namespace observations
{

    // Types of observations
    enum observationTypes {

        RANGE,
        RANGERATE,
        AZIMUTHELEVATION,
        MAX_OBSERVATIONTYPES

    };


    class base_observation {


        /**
         * Parameter struct
         */

    public:

        struct observationParams {

            // Function to get noise on measurements
            std::function<std::vector<double>()>          generateMeasurementNoise = nullptr;

        }; // observationParams


        /**
         * List of class members
         */

    protected:

        // Input parameters
        const observationParams*     m_pParams;

        // Flag for noise [false by default]
        bool                         m_noisyObservation;

        // Last noise value
        std::vector<double>          m_currentNoiseSample;

        // Observation type
        observationTypes             m_observationType;

        /**
         * Class functions
         */

    public:


        /**
         * Default constructor
         *
         */
        base_observation( const observationParams* pParams,
                          const observationTypes&  obsType = MAX_OBSERVATIONTYPES );


        /**
         * Default destructor
         *
         */
        virtual ~base_observation();


        /**
         * getObservation: Function that returns noisy measurements
         *
         * @return Measurement vector
         *
         */
        virtual std::vector<double> getNoisyObservation( const std::vector<double>& sensorState,
                                                         const std::vector<double>& targetState ) ;


        /**
         * getObservation: Function that returns perfect measurement
         *
         * @return Measurement vector
         *
         */
        virtual std::vector<double> getPerfectObservation( const std::vector<double>& sensorState,
                                                           const std::vector<double>& targetState ) = 0 ;



        /**
         * Functions to get/set member attributes
         */

    public:

        // Get observation type
        observationTypes getObservationType () const;


        /**
         * Private routines
         */

    private:

        // Constructor call: check inputs
        void checkInput() ;

        // Get noise value and save it
        std::vector<double> getNoiseSample() ;


    }; // class base_observation

} // namespace observations
} // namespace smartastro


#endif //SMART_ASTRO_BASE_SENSOR_H
