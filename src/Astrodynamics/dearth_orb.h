/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
-------Copyright (C) 2018 University of Strathclyde and Authors-------
--------- Author: Romain Serra (GitHub: Serrof)-----------------------
-------- e-mail: serra.romain@gmail.com ------------------------------
*/

#ifndef SMARTASTRO_DEARTH_ORB_H
#define SMARTASTRO_DEARTH_ORB_H

#include "base_dearth.h"
#include "dearth_3dof.h"
#include "../config.h"

#ifdef ENABLE_SMARTUQ
#include "smartuq.h"
#endif

namespace smartastro
{
    namespace astrodynamics{
        /**
         * @brief dearth_orb is a 3 degrees of freedom propagator for Earth orbits using orbital elements for integration
         *
         * 
         * This 3dof dynamics simulates Keplerian motion as well as a range of gravitational and non-graviational perturbations around the Earth. 
         * It uses the Gauss equations with the non-Kepelrian accelerations computed from dearth_3dof
         * The state vector consists of the orbital elements, nore precisely the modified equinoctial coordinates.
         * The reference frame is the Earth centered inertial.
         */
        template < class T >
        class dearth_orb: public base_dearth<T>
        {

        private:     
            double m_mu;
            using base_dearth<T>::m_L_scale;
            using base_dearth<T>::m_T_scale; 
            using base_dearth<T>::m_max_degree_Earth_gravity;
            using base_dearth<T>::m_max_degree_Earth_magnetic;
            using base_dearth<T>::m_flags;
            using base_dearth<T>::m_drag;
            using base_dearth<T>::m_other;
            using base_dearth<T>::m_C;
            using base_dearth<T>::m_S;    
            /**
             * @brief m_mass (constant) mass of the spacecraft
             */              
            T m_mass; 
            /**
             * @brief m_dyn pointer to corresponding dearth_3dof 
             */                  
            smartastro::astrodynamics::dearth_3dof<T> *m_dyn; 

        public:

            /**
             * @brief dearth_orb constructor
             * @param scaling factor in meters for distance
             * @param scaling factor in seconds for time
             * @param expansion order for geopotential
             */
            dearth_orb(const double &L_scale = 1.0, const double &t_scale = 1.0, const int &n = 0, const T &mass = 1.0, const std::vector<int> &flags=std::vector<int>(4, 0), const std::vector<T> &p = std::vector<T>(1, 0.0), const std::vector<T> &p2 = std::vector<T>(2, 0.0), const std::vector< std::vector<double> > &F10dot7 = std::vector< std::vector<double> >(2));

            ~dearth_orb();

            /**
             * @brief evaluate high-fidelity Earth-orbiting dynamics in Earth-centered inertial frame
             * @param[in] time in scaled units
             * @param[in] state vector in scaled units
             * @param[out] state derivative in scaled units
             * @return exit flag (0=success)
             */
            int evaluate(const double &t, const std::vector<T> &state, std::vector<T> &dstate) const;

            /**
             * @brief gauss_modeq Gauss equations for modified equinoctial elements
             * @param[in] x state vector in scaled units
             * @param[in] mu planetary gravitational constant 
             * @param[in] perturb vector of accelerations from orbital perturbations
             * @param[out] dx time derivative of state vector in scaled units
             * @return exit flag (0=success)
             */
            static int gauss_modeq(const std::vector<T> &x, const double &mu, const std::vector<T> &perturb, std::vector<T> &dx);

            /**
             * @brief charging method to compute the charging parameters of the S/C
             * @param[in] jd Julian date             
             * @param[in] x state vector in scaled units
             * @param[out] C electrostatic capacitance in standard unit
             * @param[out] V surface potential in standard unit
             * @return exit flag (0=success)
             */
            int charging(const double &jd, const std::vector<T> &x, T &C, T &V) const;
        };
    }
}


#endif // SMARTASTRO_DEARTH_ORB_H
