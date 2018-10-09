/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
-------Copyright (C) 2018 University of Strathclyde and Authors-------
--------- Author: SMART team----------------------
-------- e-mail: smart@strath.ac.uk ------------------------------
*/


#ifndef SMARTASTRO_D2BP_3DOF_H
#define SMARTASTRO_D2BP_3DOF_H

#include "base_astrodynamics.h"

namespace smartastro
{
    namespace astrodynamics{
        /**
         * @brief The d2bp_3dof class models the dynamics of the Two Body Problem in cartesian coordinates. All units must be consistents
         *
         * Models the dynamics of a small object, of negligible mass, subject to the gravitational field of a secondary body of mass m1.\n
         * The equations of motion are\n
         * \f$\ddot{x} = -\frac{\mu}{r^3} \cdot x\f$\n
         */
	    template < class T >
        class d2bp_3dof: public base_astrodynamics<T>
        {
	    private:
            T m_mu;
            // astrocore::propagatable_object<T>* m_body;

        public:
            /**
             * @brief Contructor for the two body dynamical system.
             * @param mass mass of the central body
             * @param L_scale length scaling factor. Default value = 1, meters (no scaling).
             * @param T_scale time scaling factor. Default value = 1, seconds (no scaling)
             */
            d2bp_3dof(const T &mass,
                      const double &L_scale = 1.0,
                      const double &T_scale = 1.0);
            
            /**
             * @brief Contructor for the two body dynamical system.
             * @param mass mass of the central body
             * @param L_scale length scaling factor. Default value = 1, meters (no scaling).
             * @param T_scale time scaling factor. Default value = 1, seconds (no scaling)
             */
            // d2bp_3dof(astrocore::propagatable_object<T>* body,
            //           const double &L_scale = 1.0,
            //           const double &T_scale = 1.0,
            //           base_noise<T>* noise_gen = NULL,
            //           controllers::base_controller<T>* controller = NULL);
            
            ~d2bp_3dof();

            /**
             * @brief Evaluate the dynamics in a given state and time.
             * @param[in] t time istant
             * @param[in] state state values
             * @param[in] dstate derivative of the state
             * @return exit flag (0=success)
             */
            int evaluate(const double &t, const std::vector<T> &state, std::vector<T> &dstate) const;            
        };
    }
}


#endif // SMARTASTRO_D2BP_3DOF_H
