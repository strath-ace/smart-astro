/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
-------Copyright (C) 2018 University of Strathclyde and Authors-------
--------- Author: Romain Serra (GitHub: Serrof)-----------------------
-------- e-mail: serra.romain@gmail.com ------------------------------
*/

#ifndef SMARTASTRO_DEARTH_3DOF_H
#define SMARTASTRO_DEARTH_3DOF_H

#include "base_dearth.h"

namespace smartastro
{
    namespace astrodynamics{
        /**
         * @brief dearth_3dof is a 3 degrees of freedom propagator for Earth orbits
         *
         * 
         * This 3dof dynamics simulates Keplerian motion as well as a range of gravitational and non-graviational perturbations around the Earth. 
         * The state vector consists of the Cartesian coordinates and the mass for a total of 7 variables.
         * The reference frame is the Earth centered inertial.
         */
        template < class T >
        class dearth_3dof: public base_dearth<T>
        {

        private:
            using base_dearth<T>::m_L_scale;
            using base_dearth<T>::m_T_scale;
            using base_dearth<T>::m_max_degree_Earth_gravity;
            using base_dearth<T>::m_max_degree_Earth_magnetic;
            using base_dearth<T>::m_flags;
            using base_dearth<T>::m_drag;
            using base_dearth<T>::m_other;
            using base_dearth<T>::m_C;
            using base_dearth<T>::m_S;                

        public:

            using base_dearth<T>::harmonics;
            using base_dearth<T>::Earth_gravity;
            using base_dearth<T>::geodetic;
            using base_dearth<T>::density;
            using base_dearth<T>::density_log;
            using base_dearth<T>::He_number_density_log;
            using base_dearth<T>::density_log_sigmoid_Tinf;
            using base_dearth<T>::Sun;
            using base_dearth<T>::lunisolar;
            using base_dearth<T>::shadow;
            using base_dearth<T>::Gaussian_coeff;
            using base_dearth<T>::IRI90;
            using base_dearth<T>::surface_potential_LEO;
            using base_dearth<T>::signed_square;

            /**
             * @brief dearth_3dof constructor
             * @param scaling factor in meters for distance
             * @param scaling factor in seconds for time
             * @param expansion order for geopotential
             * @param flags for perturbations: first one is drag, second is lunisolar, third is radiation pressure, fourth is magnetic force
             * @param drag parameters: drag area in square meters, drag coefficient (optional), mean solar flux (optional), geomagnetic index (optional) 
             * @param other parameters: SRP area in square meters, radiation pressure coefficient, constant surface charge in Coulomb (optional)
             * @param polynomial coefficients of empirical acceleration
             */
            dearth_3dof(const double &L_scale = 1.0, const double &t_scale = 1.0, const int &n=0, const std::vector<int> &flags=std::vector<int>(4,0), const std::vector<T> &p=std::vector<T>(1,0.0), const std::vector<T> &p2=std::vector<T>(2,0.0), const std::vector< std::vector<double> > &F10dot7 = std::vector< std::vector<double> >(2));

            ~dearth_3dof();

            /**
             * @brief evaluate high-fidelity Earth-orbiting dynamics in Earth-centered inertial frame
             * @param[in] time in scaled units
             * @param[in] state (position-velocity-mass) vector in scaled units
             * @param[out] state derivative in scaled units
             * @return exit flag (0=success)
             */
            int evaluate(const double &t, const std::vector<T> &state, std::vector<T> &dstate) const;

            /**
            * @brief computes acceleration due to atmospheric drag
            * @param[in] julian date
            * @param[in] state (position-velocity-mass) vector in scaled units
            * @param[out] acceleration in scaled units due to atmospheric drag
            * @return exit flag (0=success)
            */
            int aero(const double &d, const std::vector<T> &x, std::vector<T> &acc) const;

            /**
            * @brief computes acceleration due to solar radiation pressure
            * @param[in] julian date            
            * @param[in] state (position-velocity-mass) vector in scaled units
            * @param[out] acceleration in scaled units due to solar radiation pressure
            * @return exit flag (0=success)
            */
            int SRP(const double &jd,const std::vector<T> &x, std::vector<T> &acc) const;  

            /**
            * @brief computes acceleration due to Earth's magnetic field
            * @param[in] julian date            
            * @param[in] state (position-velocity-mass) vector in scaled units
            * @param[out] acceleration in scaled units due to magnetic force
            * @return exit flag (0=success)
            */
            int magnetism(const double &jd,const std::vector<T> &x, std::vector<T> &acc) const;          

            /**
            * @brief computes electrostatic surface charge due to ambient plasma
            * @param[in] julian date     
            * @param[in] state (position-velocity-mass) vector in scaled units
            * @param[out] electrostatic capacitance in SI units
            * @param[out] electrostatic charge in SI units
            * @return exit flag (0=success)
            */
            int charging(const double &jd, const std::vector<T> &x, T &C, T &V) const;

        };
    }
}


#endif // SMARTASTRO_DEARTH_3DOF_H
