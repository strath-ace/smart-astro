/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------ Copyright (C) 2018 University of Strathclyde and Authors ------
-------------------- e-mail: c.greco@strath.ac.uk --------------------
----------------------- Author: Cristian Greco -----------------------
*/


#ifndef SMART_ASTRO_GROUNDSTATION_H
#define SMART_ASTRO_GROUNDSTATION_H

#include <vector>
#include <string>
#include <cmath>

#include "Observations/smartastro_observations.h"

namespace smartastro
{

namespace observations
{

    // Observation Value info
    struct OutputMeasurementValue {

        // Type of measurement
        observationTypes        type;

        // Name of ground station
        std::string             station;

        // Value of observation
        std::vector<double>     value;

        // Time of observation [jd]
        double                  time;

    };

}

namespace sensors
{

    class groundStation {


        /**
         * Parameter struct
         */

    public:

        struct groundStationParams {

            // Reference frame to which ephemerides are referred
            std::string                                     name = "";

            // Latitude [-pi/2,+pi/2]
            double                                          lat = -M_PI;

            // Longitude [-pi,+pi]
            double                                          lon = -2.0*M_PI;

            // Types of measurements
            std::vector<observations::observationTypes>     obsTypes;

            // Function to get target state
            std::function<std::vector<double>(double)>      getTargetEphemeris = nullptr;

        }; // ephemerisParams


        /**
         * List of class members
         */

    protected:

        // Input parameters
        const groundStationParams*                      m_pParams;

        // Pointer to observations
        std::vector<observations::base_observation>     m_observations;


        /**
         * Class functions
         */

    public:


        /**
         * Default constructor
         *
         * @param Parameters structure
         */
        groundStation(const groundStationParams* pParams);


        /**
         * Default destructor
         *
         */
        virtual ~groundStation();


        /**
         * getMeasurements: Function that returns measurement of target for time t
         *
         * @param t: time at which the measurements are desired
         * @return Measurements at Julian Date jd
         *
         */
        virtual std::vector<observations::OutputMeasurementValue> getMeasurements( const double& jd ) ;


    private:

        /**
         * Initialise routine that create observation classes in constructor
         */
        void initialise() ;


    };
}
}


#endif //SMART_ASTRO_GROUNDSTATION_H
