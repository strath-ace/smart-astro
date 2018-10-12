/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------ Copyright (C) 2018 University of Strathclyde and Authors ------
-------------------- e-mail: c.greco@strath.ac.uk --------------------
----------------------- Author: Cristian Greco -----------------------
*/


#ifndef SMART_ASTRO_SPICEEPHEMERIS_H
#define SMART_ASTRO_SPICEEPHEMERIS_H

#include <iostream>
#include <iomanip>

#include "Ephemerides/base_ephemeris.h"
#include "Astro-Core/conversion_time.h"
#include "AstroData/SpiceKernels/spiceKernelNames.h"


#ifdef __USE_CSPICE
    #include <cspice/SpiceUsr.h>
#endif // __USE_CSPICE

namespace smartastro {
namespace ephemerides {


    class spiceEphemeris : public base_ephemeris {


        /**
         * Parameter struct
         */

    public:

        struct spiceEphemerisParams : public ephemerisParams {

            // Kernel to load body data from
            std::vector<std::string>        kernelToLoad ;

            // Name of object to get ephemerides
            std::string                     target;

            // Name of observer
            std::string                     observer;

            // Aberration correction
            std::string                     abberrationCorrection;

        }; // ephemerisParams


        /**
         * List of class members
         */

    protected:

        // Input parameters
        const spiceEphemerisParams* m_pParams;



        /**
         * Class functions
         */

    public:


        /**
         * Default constructor
         *
         * @param Input parameters as defined in spiceEphemerisParams
         */
        spiceEphemeris(const ephemerisParams* pParams);


        /**
         * Default destructor
         *
         */
        virtual ~spiceEphemeris();


        /**
         * getCartesianState: Function that returns object position-velocity for time t
         *
         * @param jd: time at which the Cartesian State is desired in jd
         * @return Cartesian state at time jd
         *
         */
        virtual std::vector<double> getCartesianState( const double& jd ) const ;


    }; // class spiceEphemeris



} // namespace smartastro
} // namespace ephemerides


#endif //SMART_ASTRO_SPICEEPHEMERIS_H
