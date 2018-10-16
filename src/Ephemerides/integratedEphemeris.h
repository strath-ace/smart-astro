/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------ Copyright (C) 2018 University of Strathclyde and Authors ------
-------------------- e-mail: c.greco@strath.ac.uk --------------------
----------------------- Author: Cristian Greco -----------------------
*/


#ifndef SMART_ASTRO_INTEGRATEDEPHEMERIS_H
#define SMART_ASTRO_INTEGRATEDEPHEMERIS_H

#include <functional>

#include "Ephemerides/base_ephemeris.h"
#include "Integrators/euler.h"


namespace smartastro {
    namespace ephemerides {


        class integratedEphemeris : public base_ephemeris {


            /**
             * Parameter struct
             */

        public:

            struct integratedEphemerisParams : public ephemerisParams {

                // Integration function -> int (const double &ti, const double &tend, const int &nsteps, const std::vector<T> &x0, std::vector<T> &xfinal)
                std::function<int(double,double,int,std::vector<double>,std::vector<double>)>       integrate;

                // Initial time
                double                                                                              ti;

                // Initial state
                std::vector<double>                                                                 xi;

                // Number of integration steps
                int                                                                                 nsteps;

            }; // ephemerisParams


            /**
             * List of class members
             */

        protected:

            // Input parameters
            const integratedEphemerisParams* m_pParams;

            // Current time
            double                           m_currT;

            // Current state
            std::vector<double>              m_currX;



            /**
             * Class functions
             */

        public:


            /**
             * Default constructor
             *
             * @param Input parameters as defined in integratedEphemerisParams
             */
            integratedEphemeris(const ephemerisParams* pParams);


            /**
             * Default destructor
             *
             */
            virtual ~integratedEphemeris();


            /**
             * getCartesianState: Function that returns object position-velocity for time t by numerical integration
             *
             * @param t: time at which the Cartesian State is desired
             * @return Cartesian state at time t
             *
             */
            virtual std::vector<double> getCartesianState( const double& t ) ;




        public:

            /**
             * Reset current time and state
             */
             void resetCurrentTimeState() ;


        }; // class integratedEphemeris



    } // namespace smartastro
} // namespace ephemerides


#endif //SMART_ASTRO_INTEGRATEDEPHEMERIS_H