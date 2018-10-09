/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
-------Copyright (C) 2018 University of Strathclyde and Authors-------
--------- Author: Romain Serra (GitHub: Serrof)-----------------------
-------- e-mail: serra.romain@gmail.com ------------------------------
*/

#ifndef SMARTASTRO_CONVERSION_FRAMES_H
#define SMARTASTRO_CONVERSION_FRAMES_H

#include <vector>
#include <cmath>
#include "../exception.h"
#include "../constants.h"

namespace smartastro{
    namespace astrocore{

        class conversion_frames{
        public:
            /**
             * @brief computes columns of rotation matrix to go from Earth-centered inertial frame to true of date
             * @param[in] Julian date
             * @param[out] rotation matrix in a vector form
             * @return exit flag (0=success)
             */
            static bool inertial_to_tod(const double &jd, std::vector<double> &PN);

            /**
             * @brief computes columns of rotation matrix to go from Earth-centered true of date to body-fixed frame
             * @param[in] Julian date
             * @param[out] rotation matrix in a vector form
             * @return exit flag (0=success)
             */
            static bool tod_to_bf(const double &jd, std::vector<double> &R);

                        /**
             * @brief computes columns of rotation matrix to go from Earth-centered inertial to body-fixed frame
             * @param[in] Julian date
             * @param[out] rotation matrix in a vector form
             * @return exit flag (0=success)
             */
            static bool inertial_to_bf(const double &jd, std::vector<double> &M);

	   };
    }
}


#endif /* SMARTASTRO_CONVERSION_FRAMES_H */

