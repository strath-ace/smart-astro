/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------ Copyright (C) 2018 University of Strathclyde and Authors ------
-------------------- e-mail: c.greco@strath.ac.uk --------------------
----------------------- Author: Cristian Greco -----------------------
*/


#ifndef SMART_ASTRO_AZIMUTHELEVATIONSENSOR_H
#define SMART_ASTRO_AZIMUTHELEVATIONSENSOR_H

#include "Sensors/base_sensor.h"
#include "Astro-Core/conversion_coordinates.h"

namespace smartastro
{
    namespace sensors
    {

        /**
         * Compute azimuth and elevation angles
         * ASSUMPTION: State are in topocentric reference frame of first body defined as:
         * x-axis pointing local north, y-axis pointing local west, z-axis pointing local zenith
         *
         */
        class azimuthElevationSensor : public base_sensor {


            /**
             * Class functions
             *
             */

        public:


            /**
             * Default constructor
             *
             */
            azimuthElevationSensor(const sensorParams* pParams);


            /**
             * Default destructor
             *
             */
            virtual ~azimuthElevationSensor();


            /**
             * getMeasurement: Function that returns measurements at time t
             *
             * @param t: time at which the measurement(time system is defined in derived classes)
             * @return Measurement vector
             *
             */
            virtual std::vector<double> getPerfectMeasurement( const double& t ) ;


        }; // class azimuthElevationSensor



        /**
         * Compute range given two state vectors in the topocentric reference frame of state1
         *
         * @return Range between two state vectors
         */
        std::vector<double> computeAzimuthElevation ( const std::vector<double>& stateTOPO1,
                                                      const std::vector<double>& stateTOPO2 );

        /**
         * Compute azimuthElevation given relative state vector in topocentric reference frame
         *
         * @return azimuthElevation
         */
        std::vector<double> computeAzimuthElevation ( const std::vector<double>& relStateTOPO );


    } // namespace sensors
} // namespace smartastro


#endif //SMART_ASTRO_AZIMUTHELEVATIONSENSOR_H
