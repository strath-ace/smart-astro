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

#include "Observations/base_observation.h"
#include "Astro-Core/conversion_coordinates.h"

namespace smartastro
{
    namespace observations
    {

        /**
         * Compute azimuth and elevation angles
         * ASSUMPTION: State are in topocentric reference frame of first body defined as:
         * x-axis pointing local north, y-axis pointing local west, z-axis pointing local zenith
         *
         */
        class azimuthElevationObservation : public base_observation {


            /**
             * Class functions
             *
             */

        public:


            /**
             * Default constructor
             *
             */
            azimuthElevationObservation(const observationParams* pParams);


            /**
             * Default destructor
             *
             */
            virtual ~azimuthElevationObservation();


            /**
             * getObservation: Function that returns measurements
             *
             * @return Measurement vector
             *
             */
            virtual std::vector<double> getPerfectObservation( const std::vector<double>& targetStateTOPO,
                                                               const std::vector<double>& sensorStateTOPO = {}) ;


        }; // class azimuthElevationObservation



        /**
         * Compute range given two state vectors in the topocentric reference frame of state1
         *
         * @return Range between two state vectors
         */
        std::vector<double> computeAzimuthElevation ( const std::vector<double>& targetTOPO,
                                                      const std::vector<double>& sensorTOPO );

        /**
         * Compute azimuthElevation given relative state vector in topocentric reference frame
         *
         * @return azimuthElevation
         */
        std::vector<double> computeAzimuthElevation ( const std::vector<double>& relStateTOPO );


    } // namespace observations
} // namespace smartastro


#endif //SMART_ASTRO_AZIMUTHELEVATIONSENSOR_H
