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
#include <string>
#include <cmath>
#include <functional>

#include "../exception.h"


namespace smartastro
{

namespace sensors
{

    // Observation Value info
    struct OutputSensor {

        // Name of sensor
        std::string             sensor;

        // Number of measurement
        unsigned int            number;

        // Time of observation [jd]
        double                  time;

        // Value of observation
        std::vector<double>     value;

    };

    class base_sensor {

        /**
         * Parameter struct
         */

    public:

        struct baseSensorParams {

            // Sensor name
            std::string                                           name;

            // Function to compute observations for input sensor state and target state
            std::vector<std::function<
                std::vector<double>(
                    std::vector<double>,std::vector<double>)>>    vGetObservations;

            // Function to get sensor state
            std::function<std::vector<double>(double)>            getSensorEphemeris = nullptr;

        }; // ephemerisParams


        /**
         * List of class members
         */

    protected:

        // Input parameters
        baseSensorParams*               m_pParams;

        // Current time
        double                          m_currTime;

        // Sensor current state
        std::vector<double>             m_currSensorState;


        /**
         * Class functions
         */

    public:


        /**
         * Default constructor
         *
         * @param Parameters structure
         */
        base_sensor(baseSensorParams* pParams);


        /**
         * Default destructor
         *
         */
        virtual ~base_sensor();


        /**
         * getMeasurements: Function that returns measurement of target state targetState
         *
         */
        virtual std::vector<OutputSensor> getMeasurements( const std::vector<double>& targetState,
                                                           const double& jd = std::nan("NaN") ) ;


    public:


        /**
         * Update class member at time jd
         *
         * @param jd = time in Julian Dates
         */
        virtual void update(const double& jd) ;



        /**
         * Private routines
         */

    private:

        /**
         * Check routine that ensures good input
         */
        void checkInput() ;


    };

} // namespace sensors
} // namespace smartastro

#endif //SMART_ASTRO_BASE_SENSOR_H
