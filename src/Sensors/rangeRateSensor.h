/* This Source Code Form is subject to the terms of the Mozilla Public
* License, v. 2.0. If a copy of the MPL was not distributed with this
* file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------ Copyright (C) 2018 University of Strathclyde and Authors ------
-------------------- e-mail: c.greco@strath.ac.uk --------------------
----------------------- Author: Cristian Greco -----------------------
*/


#ifndef SMART_ASTRO_RANGERATESENSOR_H
#define SMART_ASTRO_RANGERATESENSOR_H

#include "Sensors/base_sensor.h"

namespace smartastro
{
namespace sensors
{

    class rangeRateSensor : public base_sensor {


        /**
         * Class functions
         */

    public:


        /**
         * Default constructor
         *
         */
        rangeRateSensor(const sensorParams* pParams);


        /**
         * Default destructor
         *
         */
        virtual ~rangeRateSensor();


        /**
         * getMeasurement: Function that returns measurements at time t
         *
         * @param t: time at which the measurement(time system is defined in derived classes)
         * @return Measurement vector
         *
         */
        virtual std::vector<double> getPerfectMeasurement( const double& t ) ;


    }; // class rangeRateSensor



    /**
     * Compute range given two state vectors in the same reference frame
     *
     * @return Range between two state vectors
     */
    double computeRangeRate ( const std::vector<double>& state1,
                              const std::vector<double>& state2);

    /**
     * Compute rangeRate given relative state vector
     *
     * @return RangeRate
     */
    double computeRangeRate ( const std::vector<double>& relState);


} // namespace sensors
} // namespace smartastro


#endif //SMART_ASTRO_RANGERATESENSOR_H
