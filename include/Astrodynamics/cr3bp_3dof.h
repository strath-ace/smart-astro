/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
-------Copyright (C) 2018 University of Strathclyde and Authors-------
--------- Author: SMART team----------------------
-------- e-mail: smart@strath.ac.uk ------------------------------
*/


#ifndef SMARTASTRO_CR3BP_3DOF_TEST_H
#define SMARTASTRO_CR3BP_3DOF_TEST_H

#include "base_astrodynamics.h"

namespace smartastro
{
    namespace astrodynamics{
        /**
         * @brief The cr3bp_3dof class models the dynamics of the Circular Restricted Three Body Problem in cartesian coordinates. . All units must be consistents
         *
         * Models the dynamics of a system composed by three gravitationally interacting point masses, m1, m2, and m3.\\
         * It is supposed, however, that the third mass, m3, is much smaller than the other two so its gravitational effects will be neglected.
         * Moreover the first two masses, m1 and m2, execute circular orbits about their common center of mass.\\
         * The dynamics is given in a reference frame with origin in the center of mass of m1 and m2, x-axis directed toward m2, y-axis lies
         * in the orbital plane and z-axis is perpendicular.
         * \f$\ddot{x} = 2\Omega \dot{y}+\Omega^2 x - \frac{\mu_1}{r_1^3}(x+\pi_2 r_{12}) - \frac{\mu_2}{r_2^3}(x-\pi_1 r_{12})\f$\n
         * \f$\ddot{y} = -2\Omega\dot{x}+\Omega^2 y - \frac{\mu_1}{r_1^3} y - \frac{\mu_2}{r_2^3} y\f$\n
         * \f$\ddot{z} = -\frac{\mu_1}{r_1^3} z - \frac{\mu_2}{r_2^3} z\f$\n
         * where \f$r_{12}\f$ is the radius of the circular orbit, \f$\Omega\f$ the inertial angular velocity, \f$r_1\f$ and \f$r_2\f$
         * the distance between the s/c and the first and second mass respectively. \f$\pi_1\f$ and \f$\pi_2\f$ are the dimensionless mass ratio of the two bodies.
         */
	template < class T >
        class cr3bp_3dof_test: public base_astrodynamics<T>
        {
	    private:
            using base_astrodynamics<T>::m_name;

	    public:
            /**
             * @brief cr3bp_3dof Circular Restricted Three Body Problem constructor
             * @param body1 first body
             */
            cr3bp_3dof_test(const double &mu1) ;
            ~cr3bp_3dof_test();

            /**
             * @brief Evaluate the dynamics in a given state and time.
             * @param[in] t time istant
             * @param[in] state state values
             * @param[in] dstate derivative of the state
             * @return exit flag (0=success)
             */
            int evaluate(const double &t, const std::vector<T> &state, std::vector<T> &dstate) const;

            int gravity(const std::vector<T> &state, std::vector<T> &dstate) const;

        private:
            double m_mu1;

        };
    }
}


#endif // SMARTASTRO_CR3BP_3DOF_TEST_H
