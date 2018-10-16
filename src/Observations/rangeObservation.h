/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------ Copyright (C) 2018 University of Strathclyde and Authors ------
-------------------- e-mail: c.greco@strath.ac.uk --------------------
----------------------- Author: Cristian Greco -----------------------
*/


#ifndef SMART_ASTRO_rangeObservation_H
#define SMART_ASTRO_rangeObservation_H

#include "Observations/base_observation.h"

namespace smartastro
{
namespace observations
{


    class rangeObservation : public base_observation {


        /**
         * Class functions
         */

    public:


        /**
         * Default constructor
         *
         */
        rangeObservation(const observationParams* pParams);


        /**
         * Default destructor
         *
         */
        virtual ~rangeObservation();


        /**
         * getObservation: Function that returns measurements
         *
         * @return Measurement vector
         *
         */
        virtual std::vector<double> getPerfectObservation( const std::vector<double>& sensorState,
                                                           const std::vector<double>& targetState ) ;


    }; // class rangeObservation



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

} // namespace observations
} // namespace smartastro


#endif //SMART_ASTRO_rangeObservation_H
