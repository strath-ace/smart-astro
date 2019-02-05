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

            /**
             * @brief converts rth vector vRth to cartesian vector vCar given the spacecraft state {pos,vel} scCar in cartesian coordinates
             * @param[in]  scCar: spacecraft state {pos,vel} in cartesian coordinates
             * @param[in]  vRth : vector in rth coordinates to be converted
             * @param[out] vCar : vector in cartesian coordinates converted
             * @return exit flag (0=success)
             * 
             *
             * EXAMPLE:
             *
             *   Given a spacecraft in orbit:
             *       - we have the thrust vector in {r,t,h};
             *       - we want the thrust vector in {x,y,z}.
             *   In this case:
             *       scCar = [position, velocity] of the spacecraft in {x,y,z};
             *       vRth  = Thrust vector in {r,t,h};
             *       vCar  = Thrust vector, transformed in {x,y,z}.
             *
             * FUNCTIONS CALLED: none
             *
             * C++ conversion of Matlab function rth_carT written by:
             * - Camilla Colombo  - 03/03/2006
             * - Matteo Ceriotti  - 10/01/2007 : Revision
             * - Matteo Ceriotti  - 11/02/2008 : Help improved.
             * - Cristian Greco   - 05/02/2019 : C++ conversion
             * 
             */
            static int rth2car(const std::vector<double> &scCar, const std::vector<double> &vRth, std::vector<double> &vCar);

            /**
             * @brief converts cartesian vector vCar to rth vector vRth given the spacecraft state {pos,vel} scCar in cartesian coordinates
             * @param[in]  scCar: spacecraft state {pos,vel} in cartesian coordinates
             * @param[in]  vCar : vector in cartesian coordinates to be converted
             * @param[out] vRth : vector in rth coordinates converted
             * @return exit flag (0=success)
             *
             * car reference frame: {x,y,z}
             *   inertial reference frame
             * rth reference frame: {r,t,h}
             *   r-axis: direction of the orbit radius
             *   h-axis: direction of angular momentum
             *   t-axis: in the orbit plane, completes the reference frame (inward)
             *
             * EXAMPLE:
             *
             *   Given a spacecraft in orbit:
             *       - we have the thrust vector in {x,y,z}.
             *       - we want the thrust vector in {r,t,h};
             *   In this case:
             *       scCar = [position, velocity] of the spacecraft in {x,y,z};
             *       vCar  = Thrust vector in {x,y,z}.
             *       vRth  = Thrust vector, transformed in {r,t,h};
             *
             * FUNCTIONS CALLED: none
             *
             * C++ conversion of Matlab function rth_carT written by:
             * - Camilla Colombo                   - 03/03/2006
             * - Matteo Ceriotti                   - 10/01/2007 : Revision
             * - Matteo Ceriotti, Nicolas Croisard - 24/01/2008 : Help improved.
             * - Cristian Greco                    - 05/02/2019 : C++ conversion
             *
             */
            static int car2rth(const std::vector<double> &scCar, const std::vector<double> &vCar, std::vector<double> &vRth);


	   };
    }
}


#endif /* SMARTASTRO_CONVERSION_FRAMES_H */

