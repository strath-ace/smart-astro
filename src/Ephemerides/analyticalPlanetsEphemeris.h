/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------ Copyright (C) 2018 University of Strathclyde and Authors ------
-------------------- e-mail: c.greco@strath.ac.uk --------------------
----------------------- Author: Cristian Greco -----------------------
*/


#ifndef SMART_ASTRO_ANALYTICALPLANETSEPHEMERIS_H
#define SMART_ASTRO_ANALYTICALPLANETSEPHEMERIS_H

#include <algorithm>
#include <functional>

#include "Ephemerides/keplerEphemeris.h"
#include "Ephemerides/analytical_planets.h"
#include "AstroUtils/smartastro_utils.h"


namespace smartastro {
namespace ephemerides {


    class analyticalPlanetsEphemeris : public keplerEphemeris {



        /**
         * Parameter struct
         */

    public:

        struct analyticalPlanetsParams : public keplerEphemerisParams {

            // Planet of which ephemeris are desired
            std::string             planet;

        }; // ephemerisParams



        /**
         * List of class members
         */

    protected:

        // Input parameters
        const analyticalPlanetsParams* m_pParams;

        // Planet internal ID
        unsigned int m_ID;




        /**
         * Class functions
         */

    public:


        /**
         * Default constructor
         *
         * @param Input struct of parameters:
         *
         * Planet possible choiches:
         * MERCURY, VENUS, EARTH, MARS, JUPITER, SATURN, URANUS, NEPTUNE
         *
         * Mu shall be in [km^3/s^2] !
         *
         */
        analyticalPlanetsEphemeris(const ephemerisParams* pParams);


        /**
         * Default destructor
         *
         */
        virtual ~analyticalPlanetsEphemeris();



        /**
         * getKeplerianState: Function that returns object keplerian State for time t mjd2000
         *
         * @param mjd: modified Julian Date 2000 at which the Cartesian State is desired
         * @return Keplerian state at mjd2000
         *         a [km], e [-], i [rad], RAAN [rad], argument of pericenter [rad], true anomaly [rad]
         */
        std::vector<double> getKeplerianState( const double& mjd2000 ) const ;



        /**
         * List of private routines
         */

    private:

        /**
         * Check input planet requested
         */
        void initialise() ;


    };


} // namespace smartastro
} // namespace ephemerides

#endif //SMART_ASTRO_ANALYTICALPLANETSEPHEMERIS_H
