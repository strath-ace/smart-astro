/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
-------Copyright (C) 2018 University of Strathclyde and Authors-------
--------- Author: Romain Serra (GitHub: Serrof)-----------------------
-------- e-mail: serra.romain@gmail.com ------------------------------
*/

#ifndef SMARTASTRO_ANALYTICAL_PLANETS_H
#define SMARTASTRO_ANALYTICAL_PLANETS_H

#include "../exception.h"
#include "../constants.h"
#include "../Astro-Core/conversion_time.h"
#include <math.h>
#include <vector>

namespace smartastro{
    namespace ephemerides{

        class analytical_planets{
        public:
            /**
             * @brief analytically computes orbital elements in a Sun-centered ecliptic system for all planets in the solar system, from P. Dysli - 1977
             * @param[in] modified Julian date
             * @param[in] index of planet 
             *                          1:   Mercury
             *                          2:   Venus
             *                          3:   Earth
             *                          4:   Mars
             *                          5:   Jupiter
             *                          6:   Saturn
             *                          7:   Uranus
             *                          8:   Neptune
             * @param[out] kep orbital elements (semi-major axis is in kilometers and angles in radians, angular variable is true anomaly)
             * @return exit flag (0=success)
             */
            static bool get_orbel(const double &mjd2000, const int &index_planet, std::vector<double> &kep);

        };
    }
}

#endif // SMARTASTRO_ANALYTICAL_PLANETS_H
