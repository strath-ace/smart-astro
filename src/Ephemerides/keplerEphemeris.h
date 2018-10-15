/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------ Copyright (C) 2018 University of Strathclyde and Authors ------
-------------------- e-mail: c.greco@strath.ac.uk --------------------
----------------------- Author: Cristian Greco -----------------------
*/


#ifndef SMART_ASTRO_KEPLEREPHEMERIS_H
#define SMART_ASTRO_KEPLEREPHEMERIS_H


#include "Ephemerides/base_ephemeris.h"
#include "Astro-Core/conversion_coordinates.h"
#include "../exception.h"


namespace smartastro {
namespace ephemerides {


    class keplerEphemeris : public base_ephemeris {


        /**
         * Parameter struct
         */

    public:

        struct keplerEphemerisParams : public ephemerisParams {

            // Central body gravitational paremater
            double*        centralBodyGravitationalParameter = nullptr;

        }; // ephemerisParams


        /**
         * List of class members
         */

    protected:

        // Input parameters
        const keplerEphemerisParams* m_pParams;



        /**
         * Class functions
         */

    public:


        /**
         * Default constructor
         *
         * @param pParams: Input parameters
         * !! Be careful to mu units !!
         */
        keplerEphemeris(const ephemerisParams* pParams);


        /**
         * Default destructor
         *
         */
        virtual ~keplerEphemeris();


        /**
         * getCartesianState: Function that returns object position-velocity for time t
         *
         * Virtual function to be implemented in derived class
         *
         * @param t: time at which the Cartesian State is desired (time system is defined in derived classes)
         * @return Cartesian state at time t
         *
         */
        virtual std::vector<double> getCartesianState( const double& t ) ;


        /**
         * getCartesianState: Function that returns object position-velocity for time t
         *
         * Virtual function to be implemented in derived class
         *
         * @param t: time at which the Cartesian State is desired (time system is defined in derived classes)
         * @return Cartesian state at time t
         *
         */
        virtual std::vector<double> getKeplerianState( const double& t ) = 0;


    }; // class base_ephemeris


} // namespace smartastro
} // namespace ephemerides


#endif //SMART_ASTRO_KEPLEREPHEMERIS_H
