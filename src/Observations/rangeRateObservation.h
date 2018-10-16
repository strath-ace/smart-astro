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

#include "Observations/base_observation.h"

namespace smartastro
{
namespace observations
{

    class rangeRateObservation : public base_observation {


        /**
         * Class functions
         */

    public:


        /**
         * Default constructor
         *
         */
        rangeRateObservation(const observationParams* pParams);


        /**
         * Default destructor
         *
         */
        virtual ~rangeRateObservation();


        /**
         * getObservation: Function that returns measurements
         *
         * @return Measurement vector
         *
         */
        virtual std::vector<double> getPerfectObservation( const std::vector<double>& sensorState,
                                                           const std::vector<double>& targetState ) ;


    }; // class rangeRateObservation



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


} // namespace observations
} // namespace smartastro


#endif //SMART_ASTRO_RANGERATESENSOR_H
