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
namespace sensors
{

    class base_sensor {


        /**
         * Parameter struct
         */

    public:

        struct sensorParams {

            // Function to get noise on measurements
            std::function<std::vector<double>()>          generateMeasurementNoise          = nullptr;

            // Function to get sensor state
            std::function<std::vector<double>(double)>    getSensorEphemeris                = nullptr;

            // Function to get target state
            std::function<std::vector<double>(double)>    getTargetEphemeris                = nullptr;

            // Function to get sensor-target relative state
            std::function<std::vector<double>(double)>    getSensorTargetRelativeEphemeris  = nullptr;

        }; // sensorParams


        /**
         * List of class members
         */

    protected:

        // Input parameters
        const sensorParams*     m_pParams;

        // Flag for noise [false by default]
        bool                    m_noisyMeasurement;

        // Flag to check absolute ephemerides [priority]
        bool                    m_absoluteEphemeris;

        // Flag to check relative ephemerides
        bool                    m_relativeEphemeris;

        // Last noise value
        std::vector<double>     m_currentNoiseSample;


        /**
         * Class functions
         */

    public:


        /**
         * Default constructor
         *
         */
        base_sensor(const sensorParams* pParams);


        /**
         * Default destructor
         *
         */
        virtual ~base_sensor();


        /**
         * getMeasurement: Function that returns noisy measurements at time t [could be overloaded to include other effects]
         *
         * @param t: time at which the measurement(time system is defined in derived classes)
         * @return Measurement vector
         *
         */
        virtual std::vector<double> getMeasurement( const double& t ) ;


        /**
         * getMeasurement: Function that returns perfect measurements at time t
         *
         * @param t: time at which the measurement(time system is defined in derived classes)
         * @return Measurement vector
         *
         */
        virtual std::vector<double> getPerfectMeasurement( const double& t ) = 0 ;



        /**
         * Functions to set member attributes
         */

    public:

        // Set noisy measurements
        void setNoisyMeasurement( const bool noisyMeasurement );



        /**
         * Private routines
         */

    private:

        // Constructor call: check whether relative or absolute positions shall be used
        void initialise() ;

        // Get noise value and save it
        std::vector<double> getNoiseSample() ;


    }; // class base_sensor

} // namespace sensors
} // namespace smartastro


#endif //SMART_ASTRO_BASE_SENSOR_H
