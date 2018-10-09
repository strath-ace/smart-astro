/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
-------Copyright (C) 2018 University of Strathclyde and Authors-------
--------- Author: Romain Serra (GitHub: Serrof)-----------------------
-------- e-mail: serra.romain@gmail.com ------------------------------
*/

#ifndef SMARTASTRO_CONSTANT_THRUST_H
#define SMARTASTRO_CONSTANT_THRUST_H

#include "base_astrodynamics.h"

namespace smartastro
{
    namespace astrodynamics{

        template < class T >
        class constant_thrust: public base_astrodynamics<T>
        {

        protected:
            using base_astrodynamics<T>::m_name;
            using base_astrodynamics<T>::m_L_scale;
            using base_astrodynamics<T>::m_T_scale;
            /**
             * @brief m_mu scaled gravitational constant
             */            
            double m_mu;
            /**
             * @brief m_thrust scaled constant thrust 
             */               
            std::vector<double> m_thrust;
            /**
             * @brief m_Isp scaled? ISP
             */                
            double m_Isp;
            /**
             * @brief m_inertial true if thrust is constant in inertial frame, false if constant in local orbital frame
             */                  
            bool m_inertial;
            /**
             * @brief m_rate scaled rate for mass consumption 
             */               
            T m_rate;
            
        public:
            /**
             * @brief constant_thrust constructor
             * @param L_scale scaling parameter for distances
             * @param T_scale scaling parameter for time
             * @param mu scaled gravitational constant for central body
             * @param thrust propulsive force constant in local orbital frame (radial/transversal/out-of-plane) or inertial frame, in Newtons     
             * @param Isp specific impulse in s^-1 
             * @param intertial boolean, set to true if thrust is passed in inertial coordinates 
             */
            constant_thrust(const double &L_scale, const double &t_scale, const double &mu, const std::vector<double> &thrust, const double &Isp, bool inertial=false);

            ~constant_thrust();

            /**
             * @brief evaluates Keplerian dynamics with constant thrust in local orbital or inertial frame
             * @param[in] t time
             * @param[in] x state vector in scaled units (position-velocity-mass)
             * @param[out] dx state derivative in scaled units
             * @return exit flag (0=success)
             */
            using smartmath::dynamics::base_dynamics<T>::evaluate;

            int evaluate(const double &t, const std::vector<T> &x, std::vector<T> &dx) const;
            /**
             * @brief sets attribute constant thrust in local orbital or inertial frame
             * @param[in] thrust propulsive force constant in local orbital frame (radial/transversal/out-of-plane) or inertial frame, in Newtons 
             */
            void set_thrust(const std::vector<double> &thrust);
        };
    }
}


#endif // SMARTASTRO_CONSTANT_THRUST_H
