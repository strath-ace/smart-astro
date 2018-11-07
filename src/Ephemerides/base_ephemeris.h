/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------ Copyright (C) 2018 University of Strathclyde and Authors ------
-------------------- e-mail: c.greco@strath.ac.uk --------------------
----------------------- Author: Cristian Greco -----------------------
*/

#ifndef SMART_ASTRO_BASE_EPHEMERIDES_H
#define SMART_ASTRO_BASE_EPHEMERIDES_H

#include <string>
#include <vector>

#ifdef __USE_CSPICE
    #include <cspice/SpiceUsr.h>
#endif // __USE_CSPICE

#include "../exception.h"

namespace smartastro
{
namespace ephemerides
{

    class base_ephemeris {


        /**
         * Parameter struct
         */

    public:

        struct ephemerisParams {

            // Reference frame to which ephemerides are referred
            std::string          referenceFrame = "";

            // Center of Reference frame to which ephemerides are referred
            std::vector<double>  referenceFrameCenter;

        }; // ephemerisParams


        /**
         * List of class members
         */

    protected:

        // Input parameters
        const ephemerisParams* m_pParams;



        /**
         * Class functions
         */

    public:


        /**
         * Default constructor
         *
         * @param referenceFrame: name of reference frame to which ephemerides are referred
         */
        base_ephemeris(const ephemerisParams* pParams);


        /**
         * Default destructor
         *
         */
        virtual ~base_ephemeris();


        /**
         * getCartesianState: Function that returns object position-velocity for time t
         *
         * Virtual function to be implemented in derived class
         *
         * @param t: time at which the Cartesian State is desired (time system is defined in derived classes)
         * @return Cartesian state at time t
         *
         */
        virtual std::vector<double> getCartesianState( const double& t ) = 0;

//
//        /**
//         * getCartesianState: Function that returns object position-velocity for time t
//         *
//         * @param jd: time at which the Cartesian State is desired in jd
//         * @param jd: Reference frame name of output state
//         * @param jd: Reference frame center of output state
//         * @return Cartesian state at time jd
//         *
//         */
//        virtual std::vector<double> getCartesianState( const double& jd,
//                                                       const std::string& outputReferenceFrame,
//                                                       const std::vector<double>& outputReferenceFrameCenter ={} ) ;
//



    }; // class base_ephemeris


} // namespace ephemerides
} // namespace smartastro




#endif //SMART_ASTRO_BASE_EPHEMERIDES_H
