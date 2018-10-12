/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------ Copyright (C) 2018 University of Strathclyde and Authors ------
-------------------- e-mail: c.greco@strath.ac.uk --------------------
----------------------- Author: Cristian Greco -----------------------
*/


#ifndef SMART_ASTRO_RANGESENSOR_H
#define SMART_ASTRO_RANGESENSOR_H

#include "Sensors/base_sensor.h"

namespace smartastro
{
namespace sensors
{

    class rangeSensor : public base_sensor {


        /**
         * Class functions
         */

    public:


        /**
         * Default constructor
         *
         */
        rangeSensor(const sensorParams* pParams);


        /**
         * Default destructor
         *
         */
        virtual ~rangeSensor();


        /**
         * getMeasurement: Function that returns measurements at time t
         *
         * @param t: time at which the measurement(time system is defined in derived classes)
         * @return Measurement vector
         *
         */
        virtual std::vector<double> getPerfectMeasurement( const double& t ) ;






    }; // class rangeSensor



    /**
     * Compute range given two position vectors in the same reference frame
     *
     * @return Range between two position vectors
     */
    double computeRange ( const std::vector<double>& pos1,
                          const std::vector<double>& pos2);

    /**
     * Compute range given relative position vector
     *
     * @return Range of relative position vector
     */
    double computeRange ( const std::vector<double>& relPos );

} // namespace sensors
} // namespace smartastro


#endif //SMART_ASTRO_RANGESENSOR_H
